## 2114INF

This repository contains **processed, analysis-ready tables** and **R scripts** to reproduce the figures in the 2114INF paper.

### Repository structure

- `data/tables/`: **base analysis-ready tables** (metadata, CyTOF frequencies, Olink NPX, …)
- `scripts/figures/`: figure scripts (one script per panel)
- `scripts/lib/`: shared helpers (loading tables, common functions)
- `scripts/export/`: scripts that rebuild `data/tables/` from the current workspace sources
- `output/`: locally generated figures

### Re-generate base tables

From repository root:

```bash
Rscript scripts/export/export_base_tables.R
```

### Reproduce figures

Examples:

```bash
Rscript scripts/figures/fig2/Fig2A_growth_zscores.R
Rscript scripts/figures/fig5/Fig5B_mds_panelB_timepoint.R
Rscript scripts/figures/fig6/Fig6D_cytof_topClusters_PCA_timepoint.R
```

By default, scripts write PDFs into `output/figures/manuscript/`.

### QC vs. manuscript figures (convention added 2026-09-03)

`output/figures/` is split in two:

- `output/figures/manuscript/` -- figures that are actual results, candidates
  for the paper or its supplement.
- `output/figures/qc/` -- diagnostic/sanity-check figures (batch-effect
  checks, before/after comparisons, outlier hunting, "did this correction
  actually work" plots). These document how an analysis choice was arrived
  at, but are not themselves meant for the manuscript.

When a script produces both (e.g. a corrected result plot plus a sanity-check
plot confirming the correction worked), split the two `save_pdf()` calls
across both folders in the same script rather than writing everything to one
place.
