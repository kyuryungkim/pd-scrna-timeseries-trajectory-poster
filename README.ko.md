# Static scRNA-seq snapshot에서 progression-like trajectory 재구성하기

**Pseudotime이 말할 수 있는 것과 말할 수 없는 것**

이 문서는 전자포스터 내용을 더 자세히 이해하기 위한 한국어 README입니다. 공개 여부는 나중에 결정하더라도, 현재 포스터의 논리, 분석 흐름, 결과, 해석 경계를 한 번에 확인할 수 있도록 정리했습니다.

## Abstract

파킨슨병과 같은 신경퇴행성 질환은 시간에 따라 진행됩니다. 그러나 single-cell RNA sequencing(scRNA-seq)은 각 세포를 하나의 시점에서만 측정합니다. 따라서 time-course scRNA-seq 데이터라 하더라도 동일한 세포를 시간에 따라 직접 추적한 것이 아니라, 서로 다른 시점에서 얻은 정적인 분자 snapshot의 집합입니다. 이 점은 trajectory 분석에서 중요한 문제를 만듭니다. 즉, snapshot 기반 데이터에서 progression-like cellular structure를 재구성할 수 있는지, 그리고 그 결과로 얻은 pseudotime을 어디까지 해석할 수 있는지를 구분해야 합니다.

본 분석에서는 Xu et al. (2022, JCI)의 공개 D8-D35 dopaminergic differentiation scRNA-seq time-course 데이터를 재분석했습니다. 이 데이터는 파킨슨병 cell therapy 연구 맥락에서 생성된 hPSC-derived dopaminergic neuron differentiation 데이터입니다. Trajectory inference를 통해 static single-cell snapshot을 progression-like cell-state axis 위에 배열하고, 이 inferred ordering을 알려진 experimental day 및 marker-associated expression pattern과 비교했습니다.

분석 결과, inferred pseudotime은 experimental day와 전반적으로 대응되는 progression-like ordering을 보였습니다. 이는 static snapshot 데이터에서도 differentiation 과정의 cell-state structure를 어느 정도 회복할 수 있음을 뒷받침합니다. 그러나 이 관계는 1:1이 아니었습니다. 같은 day에서 샘플링된 세포들도 서로 다른 inferred state에 위치했고, later timepoint에서는 state distribution이 더 넓게 나타났습니다. 또한 expression-matched random gene module과의 비교에서는 later pseudotime region에서 marker-specific interpretability가 낮아지는 양상이 관찰되었습니다.

종합하면, public time-course scRNA-seq snapshot은 progression-like cell-state trajectory를 재구성하는 데 사용할 수 있을 뿐 아니라, 그 trajectory가 어느 구간에서 더 생물학적으로 해석 가능한지 평가하는 데에도 사용할 수 있습니다. 이 재분석은 neurodegenerative 및 rare brain disease 연구에서 trajectory method를 적용할 때 sampling time과 inferred cell-state progression을 구분하는 것이 중요함을 보여줍니다.

## 핵심 질문

Known-day differentiation time-course에서 static snapshot-derived pseudotime은 무엇을 반영할 수 있고, 어디서부터 해석을 제한해야 하는가?

## 연구 설계

이 분석은 공개 scRNA-seq differentiation time-course 데이터의 재분석입니다. 데이터에는 experimental day label이 있지만, 같은 세포를 반복 측정한 longitudinal 데이터는 아닙니다. 따라서 experimental day는 pseudotime을 비교하고 보정하는 metadata이지, 같은 세포의 실제 시간 경과를 의미하지 않습니다.

| 항목 | 값 |
|---|---|
| 공개 데이터셋 | GSE204796 |
| 원 논문 | Xu et al., JCI, 2022 |
| 시스템 | hPSC-derived dopaminergic neuron differentiation |
| Experimental day | D8, D14, D21, D28, D35 |
| Assay | scRNA-seq |
| QC 후 세포 수 | 37,163 |
| 주요 pseudotime 방법 | Slingshot |

## 분석 흐름

전자포스터의 흐름과 동일하게 다음 순서로 분석했습니다.

1. Public D8-D35 scRNA-seq 데이터를 불러옵니다.
2. QC, normalization, clustering, UMAP embedding을 수행합니다.
3. Cluster-level marker expression과 marker/module score를 확인합니다.
4. Cluster 구조와 marker pattern을 바탕으로 focused trajectory candidate subset을 정의합니다.
5. Focused subset 안에서 Slingshot으로 pseudotime을 추론합니다.
6. Inferred pseudotime을 known experimental day와 비교합니다.
7. 세포를 pseudotime rank 기반 bin으로 나누고 day-enriched pseudotime region을 확인합니다.
8. Marker/module trend를 expression-matched random gene module과 비교합니다.
9. 이 분석이 지지하는 주장과 지지하지 않는 주장을 분리합니다.

공개용 workflow 설명은 [docs/public_workflow.md](docs/public_workflow.md)에 정리되어 있습니다.

## Focused Trajectory Subset

Pseudotime을 추론하기 전에 cluster map과 marker pattern을 보고 focused trajectory subset을 정의했습니다.

| Group | Cluster | 사용 여부 |
|---|---:|---|
| Included trajectory candidates | 0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14 | 포함 |
| Separate SOX2/NES-low, DCX-high island | 10, 11, 15, 16, 17 | 제외 |
| Broadly mixed cluster | 5 | 제외 |

이 subset은 pseudotime ordering을 위한 후보 공간입니다. Lineage tracing으로 해석하면 안 됩니다.

## 주요 결과

### 1. Pseudotime은 known experimental day와 전반적으로 정렬되었습니다

Slingshot은 여러 computational path를 반환했습니다. 포스터에서는 experimental day와 가장 강하게 정렬된 `Lineage3`를 focused pseudotime axis로 사용했습니다.

Focused path와 day의 정렬:

| 지표 | 값 |
|---|---:|
| Experimental day와 Spearman rho | 0.8067 |
| 해석 | temporal consistency, not validation of true elapsed time |

Experimental day가 늦어질수록 median inferred pseudotime도 증가했습니다.

| Day | Focused path 세포 수 | Median pseudotime |
|---|---:|---:|
| D8 | 2,724 | 24.91 |
| D14 | 6,248 | 30.53 |
| D21 | 4,217 | 42.03 |
| D28 | 2,872 | 46.79 |
| D35 | 2,801 | 47.59 |

### 2. Experimental day는 exact cell state와 같지 않았습니다

Day-by-pseudotime heatmap에서는 day-enriched pseudotime region이 보였습니다. 그러나 같은 experimental day의 세포들도 여러 pseudotime bin에 분포했습니다. 따라서 known day는 pseudotime 해석을 보정하는 기준이 될 수 있지만, pseudotime이 실제 시간과 동일하다는 증거는 아닙니다.

### 3. Marker/module trend는 해석을 도왔지만, 한계도 있었습니다

Marker/module score는 rank-based pseudotime bin에 따라 변했고, early/intermediate/later region을 해석하는 데 도움이 되었습니다. 그러나 marker-control comparison에서는 later pseudotime bin에서 marker-specific separation이 약해졌습니다. 따라서 late focused axis는 더 조심스럽게 해석해야 합니다.

## 포스터용 공개 Figure

공개용 figure는 `results/poster_figures/`에 stable filename으로 export했습니다.

| Figure | 파일 | 목적 |
|---|---|---|
| Fig. 1 | [fig1_umap_by_day.png](results/poster_figures/fig1_umap_by_day.png) | experimental day별 전체 세포 UMAP |
| Fig. 2A | [fig2_cluster_umap.png](results/poster_figures/fig2_cluster_umap.png) | cluster별 전체 세포 UMAP |
| Fig. 2B | [fig2_marker_dotplot.png](results/poster_figures/fig2_marker_dotplot.png) | cluster별 marker expression |
| Fig. 2C | [fig2_focused_subset.png](results/poster_figures/fig2_focused_subset.png) | trajectory subset 포함/제외 표시 |
| Fig. 3A | [fig3_pseudotime_umap.png](results/poster_figures/fig3_pseudotime_umap.png) | pseudotime으로 색칠한 focused subset |
| Fig. 3B | [fig3_pseudotime_vs_day.png](results/poster_figures/fig3_pseudotime_vs_day.png) | pseudotime과 experimental day 비교 |
| Fig. 3C | [fig3_day_pseudotime_heatmap.png](results/poster_figures/fig3_day_pseudotime_heatmap.png) | pseudotime bin별 day 분포 |
| Fig. 4A | [fig4_module_trends.png](results/poster_figures/fig4_module_trends.png) | marker/module dynamics |
| Fig. 4B | [fig4_marker_control_separation.png](results/poster_figures/fig4_marker_control_separation.png) | marker-control separation |

## 포스터용 공개 Table

공개용 table은 `results/poster_tables/`에 stable filename으로 export했습니다.

| Table | 파일 |
|---|---|
| Sample metadata | [samples_used.csv](results/poster_tables/samples_used.csv) |
| Focused subset definition | [trajectory_subset_definition.txt](results/poster_tables/trajectory_subset_definition.txt) |
| Pseudotime path alignment | [pseudotime_path_day_alignment.csv](results/poster_tables/pseudotime_path_day_alignment.csv) |
| Pseudotime by day summary | [pseudotime_by_day_summary.csv](results/poster_tables/pseudotime_by_day_summary.csv) |
| Spearman alignment statistic | [pseudotime_day_alignment_stats.txt](results/poster_tables/pseudotime_day_alignment_stats.txt) |
| Day-by-pseudotime bin fraction | [day_pseudotime_bin_fraction.csv](results/poster_tables/day_pseudotime_bin_fraction.csv) |
| Marker-control delta by bin | [marker_control_delta_by_pseudotime_bin.csv](results/poster_tables/marker_control_delta_by_pseudotime_bin.csv) |
| Marker-control region summary | [marker_control_early_middle_late_summary.csv](results/poster_tables/marker_control_early_middle_late_summary.csv) |

## 해석 경계

이 분석이 지지하는 것:

- focused candidate differentiation-state space 안에서의 transcriptomic ordering
- known experimental day와의 temporal consistency
- marker/module pattern 기반 pseudotime region 해석
- public known-day time-course를 calibration case로 사용하는 접근

이 분석이 지지하지 않는 것:

- same-cell tracking
- exact chronological time
- true lineage tracing
- direct Parkinson's disease progression
- 원 논문 biological claim의 검증
- 원 논문 맥락을 넘어서는 새로운 marker discovery

## Technical Notes

- Experimental day label은 sampling time이지 continuous cell-state progression의 직접 측정값이 아닙니다.
- Day-associated batch effect, cell-cycle effect, scRNA-seq sparsity가 inferred structure와 late marker-control separation 감소에 영향을 줄 수 있습니다.
- RNA velocity는 사용하지 않았습니다. 이 분석은 static scRNA-seq snapshot과 marker/module 기반 check에 의도적으로 제한되어 있습니다.
- Pseudotime bin은 inferred pseudotime rank 기반 region이며 chronological interval이 아닙니다.

## 재현 방법

Raw GEO data와 intermediate R object는 이 repo에 저장하지 않습니다. 공개 사용자는 GEO에서 source data를 받은 뒤 script를 실행해야 합니다.

주요 script:

```text
scripts/eposter_code_md_pipeline.R
scripts/eposter_supporting_interpretability.R
scripts/export_public_poster_results.R
```

Repository root에서 실행 순서:

```bash
Rscript --vanilla scripts/eposter_code_md_pipeline.R
Rscript --vanilla scripts/eposter_supporting_interpretability.R
Rscript --vanilla scripts/export_public_poster_results.R
```

마지막 명령은 timestamped run folder에서 포스터에 필요한 figure/table만 골라 `results/poster_figures/`와 `results/poster_tables/`에 stable filename으로 복사합니다.

필요한 R package는 Seurat, SingleCellExperiment, slingshot, Matrix, ggplot2, dplyr, tidyr, readr, tibble, scales 등이 포함됩니다.

## Repository Scope

이 branch는 포스터 전용입니다. Raw data, intermediate RDS object, conference template, 내부 draft note를 공개하지 않으면서도, 포스터 분석을 이해하고 재현하는 데 필요한 code와 결과 파일만 공유하는 것을 목표로 합니다.
