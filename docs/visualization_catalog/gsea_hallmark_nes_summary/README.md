# GSEA Hallmark NES Summary

## Purpose

This module visualizes pathway-level results from an already completed Hallmark
GSEA analysis. It does not recalculate GSEA, ES, NES, nominal p-values, adjusted
p-values, q-values, or gene-set statistics.

## Visual mapping

- Y-axis: Hallmark pathway.
- X-axis: `NES`.
- Point size: absolute NES, `abs(NES)`.
- Point color: nominal p-value transformed as `-log10(pvalue)`.

The continuous color scale runs from `#e07e65` for lower transformed values to
`#9b1f2a` for higher transformed values.

The legend titles are `|NES|` and `-log10(p-value)`. Pathways are sorted only by
NES (descending, with the highest NES at the top), never by p-value.

## Input requirements

The input is a completed Hallmark GSEA result table and must contain:

- `Description`
- `NES`
- `pvalue`

`ID` is preferred for exact pathway matching. If `ID` is absent, exact
`Description` matching is used as a compatibility fallback. The module stops
with `Required column pvalue missing` when nominal p-values are absent; it does
not substitute `p.adjust` or `qvalue`.

## Configuration

`config/hallmark_pathways.yml` defines the pathways to display.
`config/plot_config.yml` defines paths, output dimensions, and the fixed mapping:

```yaml
input_file: "path/to/Hallmark_GSEA_result.tsv"
pathway_file: "docs/visualization_catalog/gsea_hallmark_nes_summary/config/hallmark_pathways.yml"
sort_by: "NES"
x_variable: "NES"
size_variable: "NES"
size_transform: "abs"
color_variable: "pvalue"
color_transform: "-log10"
color_low: "#e07e65"
color_high: "#9b1f2a"
label_wrap_width: 0
output_prefix: "figures/Hallmark_NES_summary"
width: 6
height: 5
dpi: 300
```

The mapping fields are validated by the module so the chart cannot silently
fall back to `setSize`, NES color, adjusted p-values, or q-values.
Set `label_wrap_width` to a positive character width for long pathway labels;
the exported `plot_data.tsv` keeps the unwrapped description.

## Run

Run from the repository root:

```bash
Rscript docs/visualization_catalog/gsea_hallmark_nes_summary/run_plot.R \docs/visualization_catalog/gsea_hallmark_nes_summary/config/plot_config.yml
```

## Outputs

For an `output_prefix` of `figures/Hallmark_NES_summary`, the module writes:

- `figures/Hallmark_NES_summary.png`
- `figures/Hallmark_NES_summary.pdf`
- `figures/Hallmark_NES_summary_plot_data.tsv`

The TSV contains the exact filtered and NES-sorted plotting rows, including
`Description`, `NES`, `pvalue`, `abs_NES`, and `neg_log10_pvalue`. Zero nominal
p-values are clamped to the smallest positive double only for the finite color
transform; the GSEA result itself is not modified.
