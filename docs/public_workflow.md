# Public e-Poster Workflow

This document describes the public, poster-facing analysis workflow. It is a
condensed English version of the internal analysis notes and includes only the
steps needed to understand or rerun the e-poster figures.

## Scope

This workflow uses the public GSE204796 D8-D35 in vitro differentiation
single-cell RNA-seq samples as a known-day calibration case for pseudotime
interpretation.

The analysis does not claim:

- direct Parkinson's disease progression
- same-cell tracking
- exact chronological time
- true lineage tracing
- new marker discovery from the original study

## Data

Raw data are not stored in this repository. They should be downloaded from GEO:

- GSE204796
- D8, D14, D21, D28, and D35 10x-format matrices

Place the downloaded files under:

```text
data/raw/
```

The repository `.gitignore` excludes raw data and intermediate R objects.

## Primary Analysis

Run from the repository root:

```bash
Rscript scripts/eposter_code_md_pipeline.R
```

This script performs:

1. Load D8-D35 scRNA-seq matrices.
2. Apply QC thresholds.
3. Normalize, cluster, and embed cells by UMAP.
4. Inspect marker modules and provisional differentiation regions.
5. Define the focused trajectory candidate subset.
6. Infer pseudotime using Slingshot.
7. Compare focused pseudotime with known experimental day.
8. Summarize marker/module trends and expression-matched random controls.

The script writes a timestamped folder under:

```text
results/eposter/
```

## Supporting Interpretability Analyses

Run after the primary analysis has generated intermediate objects:

```bash
Rscript scripts/eposter_supporting_interpretability.R
```

This script generates:

- experimental day distribution across rank-based pseudotime bins
- state-region breadth summaries
- UMAP dispersion summaries
- marker-control separation summaries

These analyses support cautious interpretation of pseudotime. They do not prove
exact time, lineage tracing, or disease progression.

## Public Poster Export

After choosing the final timestamped run, export only the minimal public-facing
figure and table set:

```bash
EP_FINAL_RUN_DIR=results/eposter/20260606_174928 \
Rscript scripts/export_public_poster_results.R
```

This creates stable public filenames:

```text
results/poster_figures/
results/poster_tables/
results/poster_manifest.md
```

The timestamped full run folders are useful locally, but the public repository
should expose only the stable poster subset unless additional output is needed
for interpretation.

## Main Interpretation Boundary

Static snapshot-derived pseudotime can support a progression-like
differentiation-state ordering when checked against known experimental days and
marker/module patterns. It should not be interpreted as direct time, same-cell
tracking, true lineage tracing, or Parkinson's disease progression.
