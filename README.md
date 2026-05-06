## AMPLIO – reproducible figures

This repository contains **processed, analysis-ready tables** and **R scripts** to reproduce the figures in the AMPLIO paper.

### Repository structure

- `data/tables/`: **base analysis-ready tables** (metadata, CyTOF frequencies, Olink NPX, …)
- `scripts/figures/`: figure scripts (one script per panel)
- `scripts/lib/`: shared helpers (loading tables, common functions)
- `scripts/export/`: scripts that rebuild `data/tables/` from the current workspace sources
- `output/`: locally generated figures (not intended for GitHub)
- `scripts/amplio_merge.Rmd`, `scripts/amplio_FCT.Rmd`: legacy analysis reports (reference only)

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
