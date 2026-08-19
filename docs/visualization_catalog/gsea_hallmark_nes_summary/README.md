# GSEA Hallmark NES Summary

## Purpose

This is a lightweight post-GSEA visualization module. It visualizes pathway-level results from Hallmark GSEA, using an already completed `*_GSEA_Hallmark_all.tsv` result to plot NES values for a user-selected pathway set. Core statistics are `NES`, `pvalue`, `p.adjust`, and `setSize`; NES is GSEA-specific.

This module is specific to Hallmark GSEA results and should not be used for Hallmark ORA results, because ORA does not produce NES. A future ORA visualization can use a separate module such as `ora_hallmark_summary`.

The module does not recalculate GSEA, ES, NES, p values, adjusted p values, or gene set statistics.

## Module responsibilities

The module keeps pathway selection, plotting configuration, implementation, and execution separate:

- `config/hallmark_pathways.yml` — `WHAT to plot`. Its `hallmark_set` is the only pathway selection. Every listed pathway is retained when it matches the input; no second truncation is applied.
- `config/plot_config.yml` — `HOW / WHERE`. It defines the input TSV, pathway YAML path, sorting, color and size mappings, output prefix, and dimensions.
- `code_template.R` — implementation. It receives a completed GSEA `data.frame`, filters the configured pathways, checks matches, removes invalid NES values, sorts, reorders the display factor, plots, and writes outputs.
- `run_plot.R` — execution. It reads the plot configuration and input TSV, sources the function template, and runs one configured plot.

## Workflow

Existing Hallmark GSEA result  
↓  
`plot_config.yml` provides input location  
↓  
`hallmark_pathways.yml` provides pathway selection  
↓  
Filter selected pathways  
↓  
Check requested, matched, and missing pathways  
↓  
Remove invalid NES  
↓  
Sort  
↓  
Factor reorder  
↓  
Plot  
↓  
PNG / PDF / `plot_data.tsv`

`gsea_hallmark_nes_summary` provides a pathway-level Hallmark GSEA overview. The separate `gsea_enrichment_curve` module provides a detailed running enrichment trajectory for selected pathways.

## Input requirements

The input TSV is read as a data frame and must contain:

- `ID`
- `Description`
- `NES`

Optional columns are:

- `pvalue`
- `p.adjust`
- `qvalue`
- `setSize`

Matching uses `ID` first with exact membership in `hallmark_set`. If `ID` is absent, exact `Description` membership is used as an explicit compatibility fallback. No fuzzy or partial matching is used.

## Pathway selection versus sorting

Pathway selection is controlled only by `config/hallmark_pathways.yml`; it determines `WHAT to plot`.

Sorting is controlled by `sort_by` in `config/plot_config.yml`; it determines the `ORDER` of the selected pathways and does not remove any selected pathway.

Supported sorting:

- `sort_by: "NES"` — NES descending.
- `sort_by: "p.adjust"` — p.adjust ascending. This option stops if `p.adjust` is absent.

The sorted first pathway is shown at the top by reversing the factor levels used on the y-axis.

## Visual mapping

- X-axis: `NES`, labeled `Normalized enrichment score (NES)`.
- Y-axis: pathway display label.
- `color_by: "NES"` — diverging continuous color centered at zero.
- `color_by: "p.adjust"` — color uses `-log10(p.adjust)`; zero values are clamped to a small positive value before the log transform.
- `color_by: "none"` — fixed point color.
- `size_by: "setSize"` — point size uses `setSize`; if absent, a warning is issued and fixed point size is used.
- `size_by: "none"` — fixed point size.

Display labels are separate from pathway identity. The original `ID` is not changed. If `Description` is already readable it is preferred; otherwise an ID-like label is displayed without the `HALLMARK_` prefix and with underscores replaced by spaces.

## Configuration template

`config/plot_config.yml`:

```yaml
input_file: "path/to/Hallmark_GSEA_result.tsv"
pathway_file: "docs/visualization_catalog/gsea_hallmark_nes_summary/config/hallmark_pathways.yml"
sort_by: "NES"
color_by: "NES"
size_by: "setSize"
output_prefix: "figures/Hallmark_NES_summary"
width: 7
height: 5
dpi: 300
```

The pathway list belongs in `config/hallmark_pathways.yml`, not in this plot configuration.

## Recommended CLI usage

Run from the repository root:

```bash
Rscript docs/visualization_catalog/gsea_hallmark_nes_summary/run_plot.R \
  docs/visualization_catalog/gsea_hallmark_nes_summary/config/plot_config.yml
```

The CLI requires one plot configuration YAML. `sort_by`, `color_by`, and `size_by` default to `NES`, `NES`, and `setSize`; `width`, `height`, and `dpi` default to `7`, `5`, and `300` when omitted.

## Direct function usage

Advanced users may source the implementation and call the function directly:

```r
source(
  "docs/visualization_catalog/gsea_hallmark_nes_summary/code_template.R"
)

gsea_result <- read.delim(
  "path/to/Hallmark_GSEA_result.tsv",
  check.names = FALSE
)

plot_hallmark_nes_summary(
  gsea_result = gsea_result,
  pathway_file =
    "docs/visualization_catalog/gsea_hallmark_nes_summary/config/hallmark_pathways.yml",
  sort_by = "NES",
  color_by = "NES",
  size_by = "setSize",
  output_prefix = "figures/Hallmark_NES_summary",
  width = 7,
  height = 5,
  dpi = 300
)
```

## Outputs

For `output_prefix = "figures/Hallmark_NES_summary"`, the module writes:

- `figures/Hallmark_NES_summary.png` at the configured DPI.
- `figures/Hallmark_NES_summary.pdf` as a vector graphic.
- `figures/Hallmark_NES_summary_plot_data.tsv` containing the finite, YAML-filtered, sorted rows actually used in the plot. It retains `ID`, `Description`, `NES`, and any available optional result columns.

No biological interpretation is added by this module.
