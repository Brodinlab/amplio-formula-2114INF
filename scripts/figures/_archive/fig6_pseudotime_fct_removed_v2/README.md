# Archived: pseudotime-by-FCT-trajectory analysis (removed from manuscript in Sep2026 v2)

These 4 scripts implement the "immune pseudotime varies by fecal-community-type (FCT)
trajectory" analysis (EF-associated trajectories T2/T3 showing lower pseudotime at 2 and
4 months). This claim was in the Aug2026 manuscript draft but was **removed** in the
Sep2026 v2 revision — the manuscript now only claims FCT trajectory associates with fecal
markers (pH, secretory IgA, antitrypsin-related measures), not with systemic immune
pseudotime. This also brought the manuscript in line with the lab's independent
WGCNA re-analysis, which found no FCT association for the flagged CD8 module
(Kruskal-Wallis p=0.83).

Archived here (not deleted) so the analysis history is preserved. See
`projects/amplio.md` in the Brodin wiki for the full narrative.

- `Fig6D_cytof_topClusters_PCA_pseudotime_proxy.R`
- `Fig6E_olink_topProteins_PCA_pseudotime_proxy.R`
- `Fig6F_cytof_delta_pseudotime_density_faceted_FCT.R`
- `Fig6G_olink_delta_pseudotime_density_faceted_FCT.R`

The feeding-group (not FCT) version of the pseudotime analysis is still current and
now lives in `scripts/figures/fig5/` (`fig5_pseudotime_*_by_feeding.R`), since it
supports the manuscript's current Fig. 5 pseudotime panel.
