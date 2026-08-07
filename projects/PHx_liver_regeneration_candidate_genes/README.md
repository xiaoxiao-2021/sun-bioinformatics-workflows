# PHx liver regeneration candidate genes

## 分析目标

Compare the expression trajectories of candidate genes during mouse liver regeneration after partial hepatectomy (PHx), distinguish PHx-associated changes from sham-surgery responses, and validate promising patterns across independent datasets.

比较小鼠部分肝切除术后肝再生过程中候选基因的表达轨迹，区分 PHx 相关变化与假手术反应，并在独立数据集中验证有潜力的表达模式。

## Current candidates

- **Adcy7** - matrix symbol `Adcy7`, ENSMUSG00000031659
- **Piezo1** - current symbol `Piezo1`; the mm9-era matrix uses the former symbol `Fam38a`, ENSMUSG00000014444
- **Cxcl2** - matrix symbol `Cxcl2`, ENSMUSG00000058427
- **Ifrd1** - matrix symbol `Ifrd1`, ENSMUSG00000001627; the paper's Figure 1A candidate

## 当前状态

进行中。GSE95135 的数据清单、R/Python 分析入口和候选基因分析框架已经建立；跨数据集整合与最终候选基因结论仍待后续验证。

## Project organization

This is a topic-level project. Each public dataset is stored under `datasets/<accession>/`; scripts and outputs use matching accession folders. Cross-dataset validation belongs in `scripts/integration/`, with its numeric outputs under the appropriate `results/` category.

- `metadata/dataset_manifest.tsv`: project-level registry of discovery and validation datasets.
- `datasets/<accession>/raw/`: original downloaded data; ignored by Git by default.
- `datasets/<accession>/processed/`: cleaned expression tables that remain inputs to later analysis.
- `datasets/<accession>/metadata/`: sample identifiers, groups, time points and source annotations.
- `scripts/<accession>/`: dataset-specific analysis code.
- `results/summary/<accession>/`: descriptive time-course summaries.
- `results/statistics/<accession>/`: statistical comparisons, P values, FDR and model results.
- `figures/<accession>/`: all PDF, PNG and SVG plots.
- `objects/`: reusable RDS, RData or h5ad objects when needed.
- `docs/`: project notes, methods and result interpretation when needed.

## 数据来源

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95135
- Organism: *Mus musculus*
- Assay: bulk liver RNA-seq
- Expression unit: author-provided log2 RPKM
- Design: untreated control, sham surgery, and post-PHx time-course samples

## 使用方法

The R scripts are the primary, human-readable workflow. They contain Chinese comments and separate data preparation, plotting, matched-time comparison, and trajectory testing.

```bash
# Run once to prepare the candidate-gene table
Rscript scripts/GSE95135/01_prepare_GSE95135.R

# Plot any one extracted gene; change Piezo1 to Ifrd1, Adcy7, or Cxcl2
Rscript scripts/GSE95135/02_plot_single_gene.R Piezo1

# Compare Post-PH and Sham only at genuinely matched time points
Rscript scripts/GSE95135/03_compare_PHx_vs_sham.R

# Exploratory group-by-time interaction test
Rscript scripts/GSE95135/04_time_group_interaction.R
```

Only `ggplot2` is required beyond base R. Install it once with `install.packages("ggplot2")`.

The R plots deliberately do not connect measured points separated by more than 12 hours. In particular, the missing Sham measurements at 28, 36, 44, and 72 h are not visually interpolated.

## Previous Python workflow

```bash
python3 scripts/GSE95135/01_plot_candidate_gene_timecourses.py
```

The Python script remains as a previously validated reference implementation. The R workflow is intended for routine use and modification.

## Design note

Control and Sham are not pooled. Untreated Control 0 h is used as a shared baseline; subsequent reference points come from Sham samples. Missing Sham time points are not imputed.
