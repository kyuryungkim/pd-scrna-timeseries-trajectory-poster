# Reconstructing Progression-like Trajectories from Static scRNA-seq Snapshots

**Interpreting what pseudotime can and cannot tell us**

This repository contains the public code and poster-facing results for an e-poster project on static scRNA-seq pseudotime interpretation. The analysis uses a public known-day dopaminergic differentiation time-course as a calibration case, not as direct evidence of Parkinson's disease progression.

## Abstract

Neurodegenerative disorders, such as Parkinson's disease, progress over time. However, single-cell RNA sequencing (scRNA-seq) captures each cell at a single time point. Consequently, even time-course scRNA-seq data provide static molecular snapshots rather than direct trajectories of the same cells. This presents a significant challenge for trajectory analysis: reconstructing progression-like cellular structures from snapshot-based data while determining how the resulting pseudotime should be interpreted.

We reanalyzed the public D8-D35 dopaminergic differentiation scRNA-seq time-course from Xu et al. (2022, JCI), which was generated in the context of Parkinson's disease cell therapy research. Through trajectory inference, we organized static single-cell snapshots along a progression-like cell-state axis and compared the inferred ordering with known experimental days and marker-associated expression patterns.

The inferred pseudotime showed a progression-like ordering that broadly corresponded to experimental days, supporting the use of static snapshots to recover cell-state structures across differentiation. However, this relationship was not strictly one-to-one: cells sampled on the same day occupied different inferred states, and later timepoints exhibited broader state distributions. Comparisons with expression-matched random gene modules further indicated reduced marker-specific interpretability in later pseudotime regions.

Collectively, these findings show that public time-course scRNA-seq snapshots can be used not only to reconstruct progression-like cell-state trajectories but also to assess where these trajectories are most biologically interpretable. This reanalysis underscores the importance of distinguishing between sampling time and inferred cell-state progression when applying trajectory methods to neurodegenerative and rare brain disease research.

## Core Question

What can static snapshot-derived pseudotime reflect from a known-day differentiation time-course, and where should its interpretation remain limited?

## Study Design

This is a reanalysis of a public scRNA-seq differentiation time-course. The dataset has known experimental day labels, but cells were sampled destructively at each timepoint. Therefore, experimental day is used as comparison metadata, not as same-cell longitudinal ground truth.

| Item | Value |
|---|---|
| Public dataset | GSE204796 |
| Source study | Xu et al., JCI, 2022 |
| System | hPSC-derived dopaminergic neuron differentiation |
| Experimental days | D8, D14, D21, D28, D35 |
| Assay | scRNA-seq |
| QC-filtered cells | 37,163 |
| Main pseudotime method | Slingshot |

## Analysis Overview

The public workflow follows the same logic as the e-poster:

1. Load public D8-D35 scRNA-seq data.
2. Perform QC, normalization, clustering, and UMAP embedding.
3. Inspect cluster-level marker expression and marker/module scores.
4. Define a focused trajectory candidate subset from cluster and marker structure.
5. Infer pseudotime within the focused subset using Slingshot.
6. Compare inferred pseudotime with known experimental days.
7. Bin cells by pseudotime rank and evaluate day-enriched pseudotime regions.
8. Compare marker/module trends with expression-matched random gene modules.
9. State what the analysis supports and what it does not support.

The detailed public workflow is summarized in [docs/public_workflow.md](docs/public_workflow.md).

## Focused Trajectory Subset

The focused trajectory subset was defined before pseudotime inference from the cluster map and marker patterns.

| Group | Clusters | Use |
|---|---:|---|
| Included trajectory candidates | 0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14 | included |
| Separate SOX2/NES-low, DCX-high island | 10, 11, 15, 16, 17 | excluded |
| Broadly mixed cluster | 5 | excluded |

This subset is a candidate ordering space. It should not be interpreted as lineage tracing.

## Main Results

### 1. Pseudotime aligned with known experimental day

Slingshot returned multiple computational paths. The poster axis used `Lineage3`, which showed the strongest alignment with experimental day among the returned paths.

Focused path day alignment:

| Metric | Value |
|---|---:|
| Spearman rho with experimental day | 0.8067 |
| Interpretation | temporal consistency, not validation of true elapsed time |

Median inferred pseudotime increased across experimental days:

| Day | Cells on focused path | Median pseudotime |
|---|---:|---:|
| D8 | 2,724 | 24.91 |
| D14 | 6,248 | 30.53 |
| D21 | 4,217 | 42.03 |
| D28 | 2,872 | 46.79 |
| D35 | 2,801 | 47.59 |

### 2. Experimental day was not equivalent to exact cell state

The day-by-pseudotime heatmap showed day-enriched pseudotime regions, but cells from the same experimental day were distributed across multiple inferred pseudotime bins. This supports calibration against known day labels, but it does not make pseudotime a direct clock.

### 3. Marker/module trends supported interpretation, with limits

Marker/module scores changed across rank-based pseudotime bins and helped interpret early, intermediate, and later regions of the focused axis. However, marker-control comparison showed weaker marker-specific separation in later pseudotime bins. This supports a cautious interpretation of the late focused axis.

## Poster-facing Figures

The stable public figure set is exported to `results/poster_figures/`.

| Figure | File | Purpose |
|---|---|---|
| Fig. 1 | [fig1_umap_by_day.png](results/poster_figures/fig1_umap_by_day.png) | all cells colored by experimental day |
| Fig. 2A | [fig2_cluster_umap.png](results/poster_figures/fig2_cluster_umap.png) | all cells colored by cluster |
| Fig. 2B | [fig2_marker_dotplot.png](results/poster_figures/fig2_marker_dotplot.png) | marker expression by cluster |
| Fig. 2C | [fig2_focused_subset.png](results/poster_figures/fig2_focused_subset.png) | included/excluded trajectory subset |
| Fig. 3A | [fig3_pseudotime_umap.png](results/poster_figures/fig3_pseudotime_umap.png) | focused subset colored by pseudotime |
| Fig. 3B | [fig3_pseudotime_vs_day.png](results/poster_figures/fig3_pseudotime_vs_day.png) | pseudotime compared with experimental day |
| Fig. 3C | [fig3_day_pseudotime_heatmap.png](results/poster_figures/fig3_day_pseudotime_heatmap.png) | day distribution across pseudotime bins |
| Fig. 4A | [fig4_module_trends.png](results/poster_figures/fig4_module_trends.png) | marker/module dynamics |
| Fig. 4B | [fig4_marker_control_separation.png](results/poster_figures/fig4_marker_control_separation.png) | marker-control separation by pseudotime bin |

## Poster-facing Tables

The stable public table set is exported to `results/poster_tables/`.

| Table | File |
|---|---|
| Sample metadata | [samples_used.csv](results/poster_tables/samples_used.csv) |
| Focused subset definition | [trajectory_subset_definition.txt](results/poster_tables/trajectory_subset_definition.txt) |
| Pseudotime path alignment | [pseudotime_path_day_alignment.csv](results/poster_tables/pseudotime_path_day_alignment.csv) |
| Pseudotime by day summary | [pseudotime_by_day_summary.csv](results/poster_tables/pseudotime_by_day_summary.csv) |
| Spearman alignment statistic | [pseudotime_day_alignment_stats.txt](results/poster_tables/pseudotime_day_alignment_stats.txt) |
| Day-by-pseudotime bin fraction | [day_pseudotime_bin_fraction.csv](results/poster_tables/day_pseudotime_bin_fraction.csv) |
| Marker-control delta by bin | [marker_control_delta_by_pseudotime_bin.csv](results/poster_tables/marker_control_delta_by_pseudotime_bin.csv) |
| Marker-control region summary | [marker_control_early_middle_late_summary.csv](results/poster_tables/marker_control_early_middle_late_summary.csv) |

## Interpretation Boundary

This analysis supports:

- transcriptomic ordering within a focused candidate differentiation-state space
- temporal consistency with known experimental days
- marker/module-based interpretation of pseudotime regions
- use of a public known-day time-course as a calibration case

This analysis does not support:

- same-cell tracking
- exact chronological time
- true lineage tracing
- direct Parkinson's disease progression
- validation of the original study's biological claims
- new marker discovery beyond the context of the original dataset

## Technical Notes

- Experimental day labels are sampling times, not direct measurements of continuous cell-state progression.
- Day-associated batch effects, cell-cycle effects, and scRNA-seq sparsity may contribute to the inferred structure and reduced late marker-control separation.
- RNA velocity was not used here; the analysis is intentionally limited to static scRNA-seq snapshots and marker/module-based checks.
- Pseudotime bins are rank-based regions of inferred pseudotime, not chronological intervals.

## Reproducibility

Raw GEO data and intermediate R objects are not stored in this repository. Public users should obtain source data from GEO and run the scripts locally.

Core scripts:

```text
scripts/eposter_code_md_pipeline.R
scripts/eposter_supporting_interpretability.R
scripts/export_public_poster_results.R
```

Typical run order from the repository root:

```bash
Rscript --vanilla scripts/eposter_code_md_pipeline.R
Rscript --vanilla scripts/eposter_supporting_interpretability.R
Rscript --vanilla scripts/export_public_poster_results.R
```

The final command copies the selected poster-facing figures and tables from the latest timestamped run into stable public paths under `results/poster_figures/` and `results/poster_tables/`.

Expected R packages include Seurat, SingleCellExperiment, slingshot, Matrix, ggplot2, dplyr, tidyr, readr, tibble, and scales.

## Repository Scope

This branch is poster-focused. It is intended to share enough code and result files to understand and reproduce the e-poster analysis without exposing raw data, intermediate RDS objects, conference template files, or internal drafting notes.
