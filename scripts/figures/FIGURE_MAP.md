# Figure map — current manuscript (`Amplio, manuscript_Sep2026.docx`, Nature Medicine submission)

Maps each manuscript figure/panel to the script(s) that generate it. Reorganized 2026-09-08
(this repo, not the upstream `nehciq/2114INF` original, is now the canonical source for
figure reproduction going forward).

## Figure 1 — Study design and cohort
No script — graphical trial-design drawing (CONSORT diagram + demographics table + sampling
schematic), not code-generated.

## Figure 2 — Growth and GI tolerance
- a–d (growth z-scores): `fig2/Fig2_growth_zscores_panels_abcd.R`
- e (IGSQ-13 GI tolerance): `fig2/Fig2E_igsq13_gi_tolerance.R`

Complete.

## Figure 3 — Markers of intestinal health
No script in this repo. **Pending: Petter to request from Nestlé.**

## Figure 4 — Microbiome and metabolite profiles
No script in this repo. **Pending: Petter to request from Nestlé.**

## Figure 5 — Immune development in relation to feeding
**Open issue:** the written legend only defines panels a–e, but Results/Methods text
cites panels up to g with conflicting letter assignments for pseudotime and the Olink
volcano (see wiki `projects/amplio.md`, 2026-09-08 entries, for detail). Scripts below
are grouped by content, not by letter, until the manuscript text is reconciled.

- CyTOF PCA colored by timepoint: `fig5/fig5_cytof_pca_timepoint.R`
- CyTOF PCA colored by feeding group (not described in the current 5-panel legend —
  confirm still wanted): `fig5/Fig5C_mds_panelA_feeding.R`
- CyTOF PCA biplot with loadings — **3 outlier-handling variants exist, canonical one
  not identified:** `Fig5_manualgating_pca_biplot_dist10_excluded.R`,
  `Fig5_manualgating_pca_biplot_filtered.R`, `Fig5_manualgating_pca_biplot_no_exclusion.R`
- Change-from-baseline effect-size ranking (best match to current legend panel c):
  `fig5/Fig5_manualgating_effect_rankings.R`
- Pseudotime by feeding group (current legend panel d) — moved here from `fig6/` on
  2026-09-08 since they compute pseudotime by `group_feeding`, not by FCT trajectory:
  `fig5/fig5_pseudotime_manualgating_by_feeding.R`,
  `fig5/fig5_pseudotime_olink_by_feeding.R`,
  `fig5/fig5_pseudotime_olink_cytokine_by_feeding.R`
  — 3 variants (CyTOF-derived, Olink-derived, Olink-cytokine-only); unclear which is the
  submission version vs. a supporting/methods-comparison variant.
- Olink volcano by feeding group: `fig5/fig5_olink_volcano_by_feeding.R` — **flagged
  mismatch:** this script plots raw log2FC vs. per-timepoint Wilcoxon p (3 panels:
  V1/V3/V5), but the current legend describes Cohen's d change-from-baseline vs.
  −log10(p) in only 2 panels (2mo/4mo). Not reconciled — see in-script note.
- Olink PCA by timepoint / by feeding group (not described in the current 5-panel
  legend — confirm still wanted): `Fig5D_olink_pca_timepoint.R`, `Fig5E_olink_pca_feeding.R`
- WGCNA six-module detection + Module 3 (cited in Results/Methods as "Fig. 5e/f/g" but
  absent from the written legend entirely):
  `scripts/analysis/wgcna_cytof_olink_modules_finegrained.R` (module detection),
  `fig5/Fig5_wgcna_module3_table.R` (kME/feature table),
  `scripts/analysis/wgcna_finegrained_brown_splitviolin.R` /
  `scripts/analysis/wgcna_module_feeding_interaction.R` (Module 3 eigengene by group)
- Remaining `Fig5_manualgating_*` / `Fig5_olink_*` variants (mds_age ×4, outlier_check ×2,
  batch_diagnostic, nlr, spaghetti_by_population, synf_ctrlf_4months_only,
  synf_ctrlf_stats ×2, cytokine_only_biplot, top_age_biplot, pca_age_clean) — not mapped
  to a specific legend panel; presumed QC/exploratory, left untouched.

## Figure 6 — Pace of immune-microbe development (FCT trajectories)
Legend narrowed in the Sep2026 v2 revision to 3 panels (transition graph / trajectory
fraction by feeding group / trajectory-vs-fecal-biomarker heatmap), removing the
pseudotime-by-FCT-trajectory claim.

- b (trajectory fraction by feeding group): `fig6/Fig6B_fct_distribution.R`
- a (transition graph), c (biomarker heatmap): **no script — not yet built. Unclear
  whether this is a Nestlé deliverable (like Fig 3/4/7) or an in-house analysis to do.**
- `fig6/Fig6D_cytof_topClusters_PCA_timepoint.R`, `fig6/Fig6E_olink_topProteins_PCA_timepoint.R`
  — no FCT dependency; likely QC/diagnostic supporting the pseudotime derivation
  methodology rather than a numbered manuscript panel. Left in place, unclassified.

**Archived** (superseded by the v2 revision, moved 2026-09-08, not deleted):
`_archive/fig6_pseudotime_fct_removed_v2/` — see that folder's README.

## Figure 7 — Upper respiratory tract infections
No script in this repo. **Pending: Petter to request from Nestlé.**

## Supplementary Figure 3 — CyTOF gating scheme
Not audited in this pass.

## Supplementary Figure 4 — Cell populations increasing/decreasing after birth
Not audited in this pass.

## Supplementary Figure 5 — Vaccine antibody responses by feeding group
`supplementary/SupplFig5_vaccine_synf_ctrlf_v5.R` — match.

## Table 1 — Adverse events
No script identified. Not audited in this pass.

## Table 2 — URTI adverse events
No script identified. Contains the internal count inconsistency flagged 2026-09-08
(19 vs 41, p=0.009 in this table vs. 25 vs 46, p=0.023 in Fig. 7/Discussion) — Petter
is removing this table from the manuscript rather than reconciling it.
