#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

repo_dir <- Sys.getenv("EP_REPO_DIR", unset = normalizePath(getwd(), mustWork = FALSE))
setwd(repo_dir)

default_out_dir <- file.path(
  "results",
  "eposter",
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_interpretability_support")
)
out_dir <- Sys.getenv("EP_SUPPORT_OUT_DIR", unset = default_out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_eposter <- function(base_size = 15) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 5),
      plot.subtitle = element_text(size = base_size - 1),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "grey20"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold")
    )
}

safe_entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

normalized_entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) <= 1) return(0)
  safe_entropy(p) / log(length(p))
}

calc_umap_dispersion <- function(obj, label) {
  emb <- as.data.frame(Embeddings(obj, "umap"))
  colnames(emb)[seq_len(2)] <- c("UMAP_1", "UMAP_2")
  df <- cbind(obj@meta.data, emb)

  centroids <- df |>
    group_by(day) |>
    summarise(
      centroid_umap_1 = mean(UMAP_1, na.rm = TRUE),
      centroid_umap_2 = mean(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )

  df_dist <- df |>
    left_join(centroids, by = "day") |>
    mutate(
      umap_distance_to_day_centroid = sqrt(
        (UMAP_1 - centroid_umap_1)^2 + (UMAP_2 - centroid_umap_2)^2
      )
    )

  summary <- df_dist |>
    group_by(day) |>
    summarise(
      dataset_scope = label,
      n_cells = n(),
      mean_umap_distance = mean(umap_distance_to_day_centroid, na.rm = TRUE),
      median_umap_distance = median(umap_distance_to_day_centroid, na.rm = TRUE),
      q25_umap_distance = quantile(umap_distance_to_day_centroid, 0.25, na.rm = TRUE),
      q75_umap_distance = quantile(umap_distance_to_day_centroid, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  list(cell_level = df_dist, summary = summary)
}

# 1. Day vs pseudotime-bin spread -------------------------------------------

obj_traj <- readRDS("data/intermediate/05_trajectory_pseudotime.rds")

traj_df <- obj_traj@meta.data |>
  mutate(
    day = factor(day, levels = c("D8", "D14", "D21", "D28", "D35"))
  ) |>
  filter(!is.na(pseudotime)) |>
  mutate(
    pseudotime_rank = percent_rank(pseudotime),
    pseudotime_bin = ntile(pseudotime_rank, 10)
  )

day_bin_fraction <- traj_df |>
  count(day, pseudotime_bin, name = "n") |>
  complete(day, pseudotime_bin = 1:10, fill = list(n = 0)) |>
  group_by(day) |>
  mutate(
    day_total = sum(n),
    fraction_within_day = n / day_total
  ) |>
  ungroup()

write_csv(day_bin_fraction, file.path(out_dir, "11_day_by_pseudotime_bin_fraction.csv"))

day_pseudotime_spread <- day_bin_fraction |>
  group_by(day) |>
  summarise(
    n_cells = first(day_total),
    n_bins_above_5pct = sum(fraction_within_day >= 0.05),
    max_bin_fraction = max(fraction_within_day),
    entropy_pseudotime_bins = safe_entropy(fraction_within_day),
    normalized_entropy_pseudotime_bins = normalized_entropy(fraction_within_day),
    .groups = "drop"
  ) |>
  left_join(
    traj_df |>
      group_by(day) |>
      summarise(
        median_pseudotime = median(pseudotime),
        q25_pseudotime = quantile(pseudotime, 0.25),
        q75_pseudotime = quantile(pseudotime, 0.75),
        iqr_pseudotime = IQR(pseudotime),
        min_pseudotime = min(pseudotime),
        max_pseudotime = max(pseudotime),
        .groups = "drop"
      ),
    by = "day"
  )

write_csv(day_pseudotime_spread, file.path(out_dir, "11_day_pseudotime_spread_summary.csv"))

p_day_bin <- ggplot(day_bin_fraction, aes(x = pseudotime_bin, y = day, fill = fraction_within_day)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(fraction_within_day >= 0.05, percent(fraction_within_day, accuracy = 1), "")),
            size = 3.8, color = "grey10") +
  scale_fill_gradient(low = "#f7fbff", high = "#2166ac", labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Experimental day distribution across pseudotime bins",
    subtitle = "Fractions are calculated within each day; spread across bins indicates day-state mismatch",
    x = "Rank-based pseudotime bin",
    y = "Experimental day",
    fill = "Within-day\nfraction"
  ) +
  theme_eposter(15)

ggsave(file.path(out_dir, "11_day_by_pseudotime_bin_heatmap.png"), p_day_bin, width = 8.8, height = 4.8, dpi = 320)

day_bin_regions <- tibble::tribble(
  ~xmin, ~xmax, ~region_label, ~region_color,
  0.5, 2.5, "D8-enriched", "#66c2a5",
  2.5, 5.5, "D14-enriched", "#8da0cb",
  5.5, 7.5, "D21-enriched", "#fc8d62",
  7.5, 10.5, "D28/D35-enriched", "#e78ac3"
)

day_bin_fraction_annotated <- day_bin_fraction |>
  mutate(day_y = as.numeric(day))

p_day_bin_annotated <- ggplot(
  day_bin_fraction_annotated,
  aes(x = pseudotime_bin, y = day_y, fill = fraction_within_day)
) +
  geom_tile(color = "white", linewidth = 0.55, height = 0.95) +
  geom_text(
    aes(label = ifelse(fraction_within_day >= 0.05, percent(fraction_within_day, accuracy = 1), "")),
    size = 3.8,
    color = "grey10"
  ) +
  geom_rect(
    data = day_bin_regions,
    aes(xmin = xmin, xmax = xmax, ymin = 5.55, ymax = 5.82, fill = NULL),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.55,
    fill = day_bin_regions$region_color
  ) +
  geom_text(
    data = day_bin_regions,
    aes(x = (xmin + xmax) / 2, y = 5.95, label = region_label),
    inherit.aes = FALSE,
    size = 3.4,
    fontface = "bold",
    color = "grey15"
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(day_bin_fraction$day)),
    labels = levels(day_bin_fraction$day),
    limits = c(0.5, 6.15),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_x_continuous(breaks = 1:10) +
  scale_fill_gradient(low = "#f7fbff", high = "#2166ac", labels = percent_format(accuracy = 1)) +
  labs(
    title = "Experimental day distribution across pseudotime bins",
    subtitle = "Colored bars indicate day-enriched pseudotime regions, not predefined chronological intervals",
    x = "Rank-based pseudotime bin",
    y = "Experimental day",
    fill = "Within-day\nfraction"
  ) +
  coord_cartesian(clip = "off") +
  theme_eposter(15) +
  theme(plot.margin = margin(t = 14, r = 12, b = 8, l = 8))

ggsave(
  file.path(out_dir, "11_day_by_pseudotime_bin_heatmap_annotated.png"),
  p_day_bin_annotated,
  width = 9.4,
  height = 5.2,
  dpi = 320
)

# 2. State-region and UMAP distribution breadth -----------------------------

obj_all <- readRDS("data/intermediate/03_marker_scored_all_cells.rds")

state_region_fraction_all <- obj_all@meta.data |>
  mutate(day = factor(day, levels = c("D8", "D14", "D21", "D28", "D35"))) |>
  count(day, state_region, name = "n") |>
  group_by(day) |>
  mutate(
    day_total = sum(n),
    fraction_within_day = n / day_total
  ) |>
  ungroup()

write_csv(state_region_fraction_all, file.path(out_dir, "12_all_cells_state_region_fraction_by_day.csv"))

state_region_display_levels <- c(
  "D8-enriched early/progenitor-like",
  "D14-enriched progenitor-pattern",
  "D21-enriched intermediate-like",
  "D21-D35 mixed intermediate-like",
  "D35-associated progenitor/regional-like",
  "SOX2/NES-low DCX-high\nneurogenic-like island",
  "Broadly mixed"
)

state_region_fraction_all_plot <- state_region_fraction_all |>
  mutate(
    state_region_display = case_when(
      state_region == "SOX2/NES-low DCX-high neurogenic island" ~
        "SOX2/NES-low DCX-high\nneurogenic-like island",
      state_region == "broadly mixed" ~ "Broadly mixed",
      TRUE ~ state_region
    ),
    state_region_display = factor(state_region_display, levels = state_region_display_levels)
  )

state_region_diversity_all <- state_region_fraction_all |>
  group_by(day) |>
  summarise(
    dataset_scope = "all_QC_cells",
    n_cells = first(day_total),
    n_state_regions_present = sum(n > 0),
    max_state_region_fraction = max(fraction_within_day),
    entropy_state_regions = safe_entropy(fraction_within_day),
    normalized_entropy_state_regions = normalized_entropy(fraction_within_day),
    .groups = "drop"
  )

cluster_fraction_all <- obj_all@meta.data |>
  mutate(day = factor(day, levels = c("D8", "D14", "D21", "D28", "D35"))) |>
  count(day, seurat_clusters, name = "n") |>
  group_by(day) |>
  mutate(
    day_total = sum(n),
    fraction_within_day = n / day_total
  ) |>
  ungroup()

cluster_diversity_all <- cluster_fraction_all |>
  group_by(day) |>
  summarise(
    dataset_scope = "all_QC_cells",
    n_clusters_present = sum(n > 0),
    max_cluster_fraction = max(fraction_within_day),
    entropy_clusters = safe_entropy(fraction_within_day),
    normalized_entropy_clusters = normalized_entropy(fraction_within_day),
    .groups = "drop"
  )

write_csv(cluster_fraction_all, file.path(out_dir, "12_all_cells_cluster_fraction_by_day.csv"))

obj_traj_meta <- obj_traj@meta.data |>
  mutate(day = factor(day, levels = c("D8", "D14", "D21", "D28", "D35")))

state_region_fraction_traj <- obj_traj_meta |>
  count(day, state_region, name = "n") |>
  group_by(day) |>
  mutate(
    day_total = sum(n),
    fraction_within_day = n / day_total
  ) |>
  ungroup()

write_csv(state_region_fraction_traj, file.path(out_dir, "12_trajectory_subset_state_region_fraction_by_day.csv"))

state_region_diversity_traj <- state_region_fraction_traj |>
  group_by(day) |>
  summarise(
    dataset_scope = "trajectory_subset",
    n_cells = first(day_total),
    n_state_regions_present = sum(n > 0),
    max_state_region_fraction = max(fraction_within_day),
    entropy_state_regions = safe_entropy(fraction_within_day),
    normalized_entropy_state_regions = normalized_entropy(fraction_within_day),
    .groups = "drop"
  )

cluster_fraction_traj <- obj_traj_meta |>
  count(day, seurat_clusters, name = "n") |>
  group_by(day) |>
  mutate(
    day_total = sum(n),
    fraction_within_day = n / day_total
  ) |>
  ungroup()

cluster_diversity_traj <- cluster_fraction_traj |>
  group_by(day) |>
  summarise(
    dataset_scope = "trajectory_subset",
    n_clusters_present = sum(n > 0),
    max_cluster_fraction = max(fraction_within_day),
    entropy_clusters = safe_entropy(fraction_within_day),
    normalized_entropy_clusters = normalized_entropy(fraction_within_day),
    .groups = "drop"
  )

diversity_summary <- bind_rows(state_region_diversity_all, state_region_diversity_traj) |>
  left_join(bind_rows(cluster_diversity_all, cluster_diversity_traj), by = c("day", "dataset_scope"))

write_csv(diversity_summary, file.path(out_dir, "12_day_state_region_cluster_diversity_summary.csv"))

p_state_heatmap <- ggplot(state_region_fraction_all_plot, aes(x = day, y = state_region_display, fill = fraction_within_day)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(fraction_within_day >= 0.05, percent(fraction_within_day, accuracy = 1), "")),
    size = 3.4,
    color = "grey12"
  ) +
  scale_fill_gradient(low = "#f7fbff", high = "#2166ac", labels = percent_format(accuracy = 1)) +
  labs(
    title = "State-region composition by experimental day",
    subtitle = "Fractions are calculated within each day across all QC cells",
    x = "Experimental day",
    y = "Provisional differentiation region",
    fill = "Within-day\nfraction"
  ) +
  theme_eposter(14) +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "12_state_region_composition_by_day_heatmap.png"), p_state_heatmap, width = 9.5, height = 5.8, dpi = 320)

disp_all <- calc_umap_dispersion(obj_all, "all_QC_cells")
disp_traj <- calc_umap_dispersion(obj_traj, "trajectory_subset")
umap_dispersion_summary <- bind_rows(disp_all$summary, disp_traj$summary)

write_csv(umap_dispersion_summary, file.path(out_dir, "12_umap_dispersion_by_day_summary.csv"))

p_disp <- ggplot(umap_dispersion_summary, aes(x = day, y = median_umap_distance, group = dataset_scope, color = dataset_scope)) +
  geom_errorbar(aes(ymin = q25_umap_distance, ymax = q75_umap_distance), width = 0.15, linewidth = 0.85, alpha = 0.75) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 3) +
  scale_color_manual(values = c(all_QC_cells = "#b2182b", trajectory_subset = "#2166ac")) +
  labs(
    title = "UMAP dispersion by experimental day",
    subtitle = "Distance from each cell to the centroid of its day group",
    x = "Experimental day",
    y = "Median UMAP distance",
    color = "Dataset scope"
  ) +
  theme_eposter(15)

ggsave(file.path(out_dir, "12_umap_dispersion_by_day.png"), p_disp, width = 7.2, height = 4.8, dpi = 320)

# 3. Marker-control interpretability summary --------------------------------

control_scores <- read_csv("data/intermediate/09_random_module_control_scores.csv", show_col_types = FALSE) |>
  filter(!is.na(pseudotime)) |>
  mutate(
    day = factor(day, levels = c("D8", "D14", "D21", "D28", "D35")),
    marker_control_delta = real_marker_module - random_mean,
    abs_marker_control_delta = abs(marker_control_delta),
    pseudotime_rank = percent_rank(pseudotime),
    pseudotime_bin = ntile(pseudotime_rank, 10),
    pseudotime_region = case_when(
      pseudotime_bin <= 3 ~ "early_bins_1_3",
      pseudotime_bin >= 8 ~ "late_bins_8_10",
      TRUE ~ "middle_bins_4_7"
    )
  )

marker_control_bin_summary <- control_scores |>
  group_by(pseudotime_bin) |>
  summarise(
    n_cells = n(),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    median_observed_marker_module = median(real_marker_module, na.rm = TRUE),
    median_random_module = median(random_mean, na.rm = TRUE),
    median_signed_delta = median(marker_control_delta, na.rm = TRUE),
    q25_signed_delta = quantile(marker_control_delta, 0.25, na.rm = TRUE),
    q75_signed_delta = quantile(marker_control_delta, 0.75, na.rm = TRUE),
    median_abs_separation = median(abs_marker_control_delta, na.rm = TRUE),
    q25_abs_separation = quantile(abs_marker_control_delta, 0.25, na.rm = TRUE),
    q75_abs_separation = quantile(abs_marker_control_delta, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(marker_control_bin_summary, file.path(out_dir, "13_marker_control_delta_by_pseudotime_bin.csv"))

marker_control_region_summary <- control_scores |>
  group_by(pseudotime_region) |>
  summarise(
    n_cells = n(),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    median_signed_delta = median(marker_control_delta, na.rm = TRUE),
    q25_signed_delta = quantile(marker_control_delta, 0.25, na.rm = TRUE),
    q75_signed_delta = quantile(marker_control_delta, 0.75, na.rm = TRUE),
    median_abs_separation = median(abs_marker_control_delta, na.rm = TRUE),
    q25_abs_separation = quantile(abs_marker_control_delta, 0.25, na.rm = TRUE),
    q75_abs_separation = quantile(abs_marker_control_delta, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(factor(pseudotime_region, levels = c("early_bins_1_3", "middle_bins_4_7", "late_bins_8_10")))

write_csv(marker_control_region_summary, file.path(out_dir, "13_marker_control_early_middle_late_summary.csv"))

p_abs <- ggplot(marker_control_bin_summary, aes(x = pseudotime_bin, y = median_abs_separation)) +
  geom_errorbar(aes(ymin = q25_abs_separation, ymax = q75_abs_separation), width = 0.12, color = "#d98c99", linewidth = 0.9) +
  geom_line(color = "#b2182b", linewidth = 1.1) +
  geom_point(color = "#b2182b", size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Marker-control separation across pseudotime bins",
    subtitle = "Lower values indicate weaker marker-specific separation from expression-matched random modules",
    x = "Rank-based pseudotime bin",
    y = "Median absolute difference"
  ) +
  theme_eposter(15)

ggsave(file.path(out_dir, "13_marker_control_abs_separation_by_pseudotime_bin.png"), p_abs, width = 8.5, height = 5.0, dpi = 320)

p_signed <- ggplot(marker_control_bin_summary, aes(x = pseudotime_bin, y = median_signed_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_signed_delta, ymax = q75_signed_delta), width = 0.12, color = "#9bb8dc", linewidth = 0.9) +
  geom_line(color = "#2166ac", linewidth = 1.1) +
  geom_point(color = "#2166ac", size = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Observed marker module minus random module",
    subtitle = "Values closer to zero indicate more similar median scores",
    x = "Rank-based pseudotime bin",
    y = "Median signed difference"
  ) +
  theme_eposter(15)

ggsave(file.path(out_dir, "13_marker_control_signed_delta_by_pseudotime_bin.png"), p_signed, width = 8.5, height = 5.0, dpi = 320)

interpretation_notes <- c(
  "Supporting interpretability analyses",
  "",
  "1. day_by_pseudotime_bin_heatmap:",
  "   Quantifies whether cells from each experimental day are concentrated in one pseudotime region or spread across multiple inferred states.",
  "",
  "2. state_region and UMAP dispersion summaries:",
  "   Quantify whether later timepoints show broader cell-state distributions. These are descriptive summaries, not causal tests.",
  "",
  "3. marker-control summaries:",
  "   Compare observed marker module scores with expression-matched random gene modules. The control is gene-level, not a biological control sample.",
  "",
  "Claim boundary:",
  "   These analyses support region-specific interpretability of pseudotime. They do not prove exact chronological time, same-cell tracking, lineage tracing, or disease progression."
)

writeLines(interpretation_notes, file.path(out_dir, "00_interpretation_support_notes.txt"))

cat("Additional interpretability-support analyses completed.\n")
cat("Output directory:", normalizePath(out_dir), "\n")

cat("\n[CHECK] Day-pseudotime spread summary\n")
print(day_pseudotime_spread, n = Inf)

cat("\n[CHECK] State-region / cluster diversity summary\n")
print(diversity_summary, n = Inf)

cat("\n[CHECK] UMAP dispersion summary\n")
print(umap_dispersion_summary, n = Inf)

cat("\n[CHECK] Marker-control early/middle/late summary\n")
print(marker_control_region_summary, n = Inf)
