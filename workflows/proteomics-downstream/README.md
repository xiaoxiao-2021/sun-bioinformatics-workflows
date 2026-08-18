# Proteomics downstream workflow

This workflow starts from an **unscreened protein-level differential-expression
CSV or Excel sheet** produced upstream. It does not perform raw-spectrum processing,
protein identification, missing-value imputation, or vendor result replacement.

## Analysis branches

- Volcano plots are emitted in two separate versions using the vendor-provided
  statistics: a q-value/FDR version and a nominal p-value version. The vendor
  q-value is retained as supplied and is never recalculated. Formal q-value and
  exploratory nominal-p-value DE tables are saved separately.
- ORA uses unique valid Entrez IDs. Its universe is the set of mapped proteins
  in the unscreened input, not the whole genome. Formal foregrounds use vendor
  q-value plus absolute log2 fold change; exploratory foregrounds use nominal
  p-value plus absolute log2 fold change.
- If the vendor sheet has no GeneID column, `01_import_and_standardize.R`
  maps the standardized `SYMBOL` to `ENTREZID` with `org.Hs.eg.db` without
  changing the original Gene Name column or dropping unmapped proteins.
- GSEA is the only branch that refits protein abundance values with limma. The
  signed moderated t statistic is collapsed to one representative protein per
  Entrez ID and used as the unfiltered ranked list.
- Hallmark provides a compact biological overview. GO Biological Process and
  KEGG provide more detailed functional results. Vendor GO/pathway annotation
  columns are retained for reference but are not used as the primary enrichment
  databases.

Sample membership comes only from the metadata TSV. The workflow verifies the
metadata against the Excel abundance columns and checks whether abundance-scale
case-minus-control differences agree with the vendor log2 fold changes. It does
not silently log-transform or otherwise change a mismatched abundance scale.

For a simplified vendor-statistics-only sheet, keep the complete configuration
below but set `sample_metadata: null`, `case_group: null`,
`control_group: null`, `gene_id_col: null`, and `run_gsea: false`. Volcano and ORA remain available;
the heatmap is skipped because there are no replicate-level abundance columns,
and the runner skips both GSEA scripts with the normal message
`GSEA skipped by config: run_gsea = false`.

## Run

From the repository root:

```bash
Rscript workflows/proteomics-downstream/run_pipeline.R \
  projects/ADCY7_proteomics/config/ADCY7_proteomics.yml
```

The seven scripts run in order. Output is written below the configured project,
and the complete console record is written to `logs/pipeline.log`. On failure,
the runner reports the failed step and preserves the original R error message.
DE plot and audit filenames include the active P/q-value, absolute log2FC, and
label/top-N settings, so rerunning a different threshold combination does not
overwrite earlier parameter sweeps.

Set `gsea_draw_curves: false` to skip individual enrichment-curve plots while
retaining the GSEA NES dot plots and lollipop plots. `gsea_curve_n` is used only
when curve drawing is enabled.

CSV inputs are read with base R; Excel inputs use `readxl`. Required R packages are `yaml`, `readxl`, `ggplot2`, `pheatmap`, `limma`,
`clusterProfiler`, and `org.Hs.eg.db`. `AnnotationDbi` is required when
SYMBOL-to-ENTREZID fallback mapping is used; `msigdbr` is required for GSEA,
and `enrichplot` is required only when `gsea_draw_curves: true`.
