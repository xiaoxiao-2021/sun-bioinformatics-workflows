# GSEA Enrichment Curve

## Purpose

This lightweight module is a post-GSEA visualization tool. It reconstructs a selected Hallmark pathway's running enrichment trajectory from three snapshots produced by the same GSEA analysis:

1. GSEA result TSV
2. ranked gene list TSV
3. Hallmark TERM2GENE snapshot TSV

The module does not rerun GSEA, recalculate formal NES/P/FDR values, or perform DEG filtering.

## Module responsibilities

- `config/selected_pathways.yml` — decides `WHAT to plot`; every listed pathway must be matched in both the GSEA result and TERM2GENE snapshot.
- `config/plot_config.yml` — provides input paths, gene/rank columns, group labels, exponent, output directory, and figure dimensions.
- `code_template.R` — rebuilds the running ES trajectory, extracts the original GSEA statistics, creates the figure, and writes curve data.
- `run_plot.R` — CLI entry point that reads the YAML files and runs one curve per selected pathway.

## Reproducibility principle

The three input snapshots should come from the same GSEA analysis:

- GSEA result → authoritative NES, nominal P, FDR, and set size.
- Ranked genes → ranked axis and ranking statistic.
- TERM2GENE → Hallmark gene-set membership.

Do not download Hallmark gene sets again during plotting. Prefer the workflow output named `*_GSEA_Hallmark_TERM2GENE.tsv`, because it records the membership used for the actual GSEA calculation.

## Data workflow

GSEA result + ranked genes + TERM2GENE + selected pathways  
↓  
Exact pathway matching  
↓  
Ranked-gene matching  
↓  
Running ES reconstruction  
↓  
Hit positions and ranked-metric strip  
↓  
Original NES/P annotation  
↓  
PNG / PDF / curve data TSV

Pathway selection and sorting are separate concepts: this module does not apply an automatic Top-N selection or remove selected pathways according to NES or P value.

## Configuration

Example `config/plot_config.yml`:

```yaml
gsea_result_file: "projects/ADCY7_proteomics/results/GSEA/ADCY7_proteomics_GSEA_Hallmark_all.tsv"
ranked_gene_file: "projects/ADCY7_proteomics/results/GSEA/ADCY7_proteomics_GSEA_ranked_genes.tsv"
term2gene_file: "projects/ADCY7_proteomics/results/GSEA/ADCY7_proteomics_GSEA_Hallmark_TERM2GENE.tsv"
pathway_file: "docs/visualization_catalog/gsea_enrichment_curve/config/selected_pathways.yml"
gene_id_column: "ENTREZID"
rank_column: "limma_t"
case_group: "ADCY7"
control_group: "NC"
gsea_exponent: 1
output_dir: "projects/ADCY7_proteomics/figures/GSEA/Hallmark_curves"
width: 7
height: 5
dpi: 300
```

Project and group names are configuration values; they are not hard-coded in the R implementation.

`selected_pathways.yml` contains only the pathway set:

```yaml
pathway_set:
  - HALLMARK_TNFA_SIGNALING_VIA_NFKB
  - HALLMARK_INFLAMMATORY_RESPONSE
  - HALLMARK_INTERFERON_ALPHA_RESPONSE
```

## Running ES reconstruction

For each selected pathway, the module intersects TERM2GENE genes with the ranked gene IDs and sorts the full ranking by the configured rank column in decreasing order.

For a hit at rank `i`, the increment is:

`abs(r_i)^p / sum(abs(r_hit)^p)`

where `p = gsea_exponent`. Each miss contributes `-1 / (N - Nh)`. The resulting running sum is used only to draw the visualization. The trajectory is expected to return to zero within floating-point error.

If the original result contains `enrichmentScore`, the reconstructed visual ES is compared with it for QC. A discrepancy produces a warning; the original GSEA statistics are never replaced.

## Recommended CLI usage

Run from the repository root:

```bash
Rscript docs/visualization_catalog/gsea_enrichment_curve/run_plot.R \docs/visualization_catalog/gsea_enrichment_curve/config/plot_config.yml
```

Defaults are `gsea_exponent = 1`, `width = 7`, `height = 5`, and `dpi = 300`. `case_group` and `control_group` must be explicitly supplied.

## Outputs

Each successfully plotted pathway is saved independently; this V0.1 module does not use patchwork, cowplot, or another multi-plot composition framework.

For a pathway such as `HALLMARK_TNFA_SIGNALING_VIA_NFKB`, the output directory receives:

- `TNFA_SIGNALING_VIA_NFKB.png`
- `TNFA_SIGNALING_VIA_NFKB.pdf`
- `TNFA_SIGNALING_VIA_NFKB_curve_data.tsv`

Each curve data TSV contains `rank`, `gene_id`, `rank_metric`, `is_hit`, and `running_ES`.

The module also writes `gsea_enrichment_curve_summary.tsv`, with one row per successfully plotted pathway and fields including `ID`, `Description`, `NES`, `pvalue`, `p.adjust`, `qvalue`, `setSize`, `matched_gene_count`, `visual_ES`, and `visual_ES_rank`.

The PNG/PDF annotation uses NES and nominal P directly from the original GSEA result. No biological interpretation is added.
