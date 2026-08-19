# Reusable bulk RNA-seq workflow

## 1. 适用范围

适用于已经获得 gene-level expression matrix 的两组 bulk RNA-seq 差异表达分析。当前面向 human、Ensembl gene ID、case vs control，包含 QC、limma、TREAT、基因注释、volcano、Top DEG heatmap，以及 GO-BP、KEGG 和 Hallmark pathway analysis。

流程同时输出 standard limma FDR DEG、nominal P-value exploratory DEG 和 strict TREAT DEG。`downstream_deg_method` 可选 `"limma"`、`"limma_pvalue"` 或 `"treat"`，分别用于 heatmap 和 ORA；模板默认使用正式的 `"limma"`，LM3 示例使用探索性的 `"limma_pvalue"`。

Volcano 会同时输出 FDR 版和 nominal P-value 探索版，文件名明确标注所用显著性指标。

Volcano label filtering: `volcano_label_filter` supports `all`, `annotated`, and `protein_coding` (default when omitted). It changes labels only; every tested gene remains in the volcano plot and keeps the original significance classification.
Heatmap gene filtering: `heatmap_gene_filter` supports `all`, `annotated`, and `protein_coding` (default). Filtering affects heatmap visualization only; all genes remain in differential-expression results and enrichment analysis.

## 2. 输入方式

- `config.yml`：项目路径、组名和分析阈值；路径相对于仓库根目录，`input_file` 和 `sample_metadata` 相对于 `project_dir`。
- sample metadata：制表符分隔，必须含唯一的 `sample` 和对应的 `group` 两列。
- expression file：Excel 第一张表；包含配置指定的 gene ID 列以及 metadata 中的全部样本列。

可从 `config_template.yml` 和 `sample_metadata_template.tsv` 复制后修改。

Enrichment outputs include all results, nominal P exploratory results, and FDR-significant results. At pathway level, `P < 0.05` is exploratory and `FDR < 0.05` is formal.

`enrichment_gene_logFC_cutoff` controls only the gene-level `|logFC|` threshold used to define the ORA foreground; `0` means no effect-size filtering. It does not change the ORA background universe.

Steps 10-11 add GO Biological Process, KEGG, and Hallmark GSEA using all limma-tested genes ranked by the moderated t statistic. Positive NES indicates enrichment in the configured case group; negative NES indicates the control group. Nominal P results are exploratory and FDR results are formal.

## 4. Hallmark pathway analysis

Hallmark is the MSigDB H gene-set collection, not a new enrichment statistic. It is analyzed with the existing methods: Hallmark ORA uses thresholded foreground genes against the tested-gene background and reports `GeneRatio`, `Count`, `pvalue`, and `p.adjust`; Hallmark GSEA uses the complete ranking and reports `ES`, `NES`, `pvalue`, and `p.adjust`.

The authoritative bulk RNA-seq GSEA ranking snapshot is `<project_dir>/results/<dataset_id>/GSEA/<dataset_id>_<case>_vs_<control>_GSEA_ranked_genes.tsv`; its `t` column is the limma moderated t statistic used in `geneList`. Hallmark GSEA writes the matching `<prefix>_GSEA_Hallmark_all.tsv` and `<prefix>_GSEA_Hallmark_TERM2GENE.tsv`. Hallmark ORA writes directional `<prefix>_ORA_Hallmark_UP_all.tsv` and `<prefix>_ORA_Hallmark_DOWN_all.tsv`, a combined `<prefix>_ORA_Hallmark_all.tsv`, and `<prefix>_ORA_Hallmark_TERM2GENE.tsv`. The TERM2GENE files are snapshots of the exact objects passed to `GSEA()` or `enricher()`.

Specialized publication-style Hallmark figures should be generated from `docs/visualization_catalog/` (for example `gsea_hallmark_nes_summary` and `gsea_enrichment_curve`), rather than hard-coded inside this bulk RNA-seq workflow.

## 3. Linux 终端使用方法

```bash
cd ~/sun-bioinformatics-workflows

cp workflows/bulk-rnaseq/config_template.yml \
  projects/U251_LM3_RNA-seq/config/LM3.yml
nano projects/U251_LM3_RNA-seq/config/LM3.yml

Rscript workflows/bulk-rnaseq/run_pipeline.R \projects/U251_LM3_RNA-seq/config/LM3.yml
```

运行日志位于 `<project_dir>/logs/<dataset_id>/pipeline.log`：

```bash
tail -f projects/U251_LM3_RNA-seq/logs/LM3/pipeline.log
```
