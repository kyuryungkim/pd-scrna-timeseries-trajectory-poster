#!/usr/bin/env Rscript

# Export a minimal, stable set of poster-facing outputs from one timestamped run.
#
# Usage from repository root:
#   Rscript scripts/export_public_poster_results.R
#
# Optional:
#   EP_FINAL_RUN_DIR=results/eposter/20260606_174928 Rscript scripts/export_public_poster_results.R

repo_dir <- Sys.getenv("EP_REPO_DIR", unset = normalizePath(getwd(), mustWork = FALSE))
setwd(repo_dir)

run_dir <- Sys.getenv("EP_FINAL_RUN_DIR", unset = "")
if (!nzchar(run_dir)) {
  runs <- list.dirs(file.path("results", "eposter"), recursive = FALSE, full.names = TRUE)
  runs <- runs[dir.exists(runs)]
  if (length(runs) == 0) {
    stop("No results/eposter run directories found.")
  }
  run_dir <- runs[which.max(file.info(runs)$mtime)]
}

if (!dir.exists(run_dir)) {
  stop("Run directory not found: ", run_dir)
}

fig_dir <- file.path("results", "poster_figures")
tab_dir <- file.path("results", "poster_tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

copy_named <- function(from, to) {
  src <- file.path(run_dir, from)
  if (!file.exists(src)) {
    stop("Expected output is missing: ", src)
  }
  ok <- file.copy(src, to, overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy ", src, " to ", to)
  }
}

figure_map <- c(
  "fig1_umap_by_day.png" = "03_all_cell_umap_by_day_clean.png",
  "fig2_cluster_umap.png" = "03_all_cell_umap_by_cluster_labeled.png",
  "fig2_marker_dotplot.png" = "04_marker_dotplot_by_cluster.png",
  "fig2_focused_subset.png" = "05_focused_trajectory_candidate_subset.png",
  "fig3_pseudotime_umap.png" = "06_trajectory_subset_umap_pseudotime.png",
  "fig3_pseudotime_vs_day.png" = "07_pseudotime_vs_day_violin.png",
  "fig3_day_pseudotime_heatmap.png" = "11_day_by_pseudotime_bin_heatmap_annotated.png",
  "fig4_module_trends.png" = "08_module_ranked_pseudotime_bins.png",
  "fig4_marker_control_separation.png" = "13_marker_control_abs_separation_by_pseudotime_bin.png"
)

table_map <- c(
  "samples_used.csv" = "01_metadata_raw_cell_counts.csv",
  "trajectory_subset_definition.txt" = "05_trajectory_subset_definition_note.txt",
  "pseudotime_path_day_alignment.csv" = "06_pseudotime_path_day_alignment.csv",
  "pseudotime_by_day_summary.csv" = "07_pseudotime_by_day_summary.csv",
  "pseudotime_day_alignment_stats.txt" = "07_pseudotime_day_alignment_stats.txt",
  "day_pseudotime_bin_fraction.csv" = "11_day_by_pseudotime_bin_fraction.csv",
  "marker_control_delta_by_pseudotime_bin.csv" = "13_marker_control_delta_by_pseudotime_bin.csv",
  "marker_control_early_middle_late_summary.csv" = "13_marker_control_early_middle_late_summary.csv"
)

for (target in names(figure_map)) {
  copy_named(figure_map[[target]], file.path(fig_dir, target))
}

for (target in names(table_map)) {
  copy_named(table_map[[target]], file.path(tab_dir, target))
}

manifest <- c(
  "# Public Poster Results Manifest",
  "",
  paste0("- Source run: `", run_dir, "`"),
  paste0("- Exported at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Figures",
  paste0("- `results/poster_figures/", names(figure_map), "`"),
  "",
  "## Tables",
  paste0("- `results/poster_tables/", names(table_map), "`"),
  "",
  "These files are a minimal poster-facing subset. Timestamped full run folders,",
  "raw data, and intermediate R objects are not required for public viewing."
)

writeLines(manifest, file.path("results", "poster_manifest.md"))

message("Exported public poster figures to: ", normalizePath(fig_dir, mustWork = FALSE))
message("Exported public poster tables to: ", normalizePath(tab_dir, mustWork = FALSE))
message("Manifest: ", normalizePath(file.path("results", "poster_manifest.md"), mustWork = FALSE))
