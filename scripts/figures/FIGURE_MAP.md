# Figure map — current manuscript (`Amplio, manuscript_Sep2026.docx`, Nature Medicine submission)

Maps each manuscript figure/panel to the script that generates it. This repo (not the
upstream `nehciq/2114INF` original) is the canonical source for figure reproduction.

**Verification method (2026-09-08):** for each panel below marked ✅ VERIFIED, the
candidate script was actually executed and its printed statistics (Cohen's d, p-values,
% variance explained, named populations/proteins, sample counts) were checked digit-for-
digit against the manuscript's own Results/Methods/Figure-legend text. A superficial
filename/content match was NOT treated as sufficient — several scripts that looked like
plausible matches on inspection turned out to use the wrong input table, wrong method, or
a filtered/exploratory subset once actually run. Those are documented under "Ruled out"
in each section and moved to `_archive/superseded/`.

Also fixed 2026-09-08: 17 scripts repo-wide had a stale `source("code/lib/...")` path
(the actual directory is `scripts/lib/`) and would not execute at all. Corrected.

## Figure 1 — Study design and cohort
No script — graphical trial-design drawing, not code-generated. (Confirmed with Petter.)

## Figure 2 — Growth and GI tolerance
- a–d (growth z-scores): `fig2/Fig2_growth_zscores_panels_abcd.R`
- e (IGSQ-13 GI tolerance): `fig2/Fig2E_igsq13_gi_tolerance.R`

**Structurally matches the legend exactly** (same 4 growth metrics + IGSQ-13 with the
same two clinical-threshold annotations), but **cannot be numerically verified**: both
scripts require `data/figure_data/fig2/*.csv`, which does not exist anywhere in this
repo. Either export that table here, or confirm these figures were finalized from data
that lives outside this repo.

## Figure 3 — Markers of intestinal health
No script in this repo. Pending: request from Nestlé.

## Figure 4 — Microbiome and metabolite profiles
No script in this repo. Pending: request from Nestlé.

## Figure 5 — Immune development in relation to feeding
### a, b — CyTOF PCA + biplot ✅ VERIFIED
`fig5/Fig5_manualgating_pca_biplot_dist10_excluded.R`
- PC1/PC2 variance explained: script gives 38.6%/10.8% vs. legend's "38.7%/10.8%" (n=248,
  3 samples excluded by the dist>10 outlier rule — matches the manuscript's n for the
  immune sub-study exactly).
- Top PC1 loadings: CD8T.TEM, CD8T.TEMRA, CD8T.TCM — matches legend's "CD8 memory/effector
  subsets — central memory, effector memory, and terminally-differentiated effector memory
  — load on PC1" exactly.
- Top PC2 loadings: Neutrophils, Nonclassical.monocytes — matches legend's "neutrophils and
  nonclassical monocytes load on PC2" exactly.
- **Open discrepancy:** this script runs on 32 populations (`base$cytof_manual`, "Kanth's
  re-gated v1.1" raw table), not the "27 manually-gated" populations stated in the legend
  and Methods. The statistics match too precisely to be the wrong script, but the
  population-count wording in the manuscript needs reconciling (32, not 27, for this
  specific panel — the CyTOF Δ-baseline and WGCNA panels genuinely do use 27, see below).
- Ruled out: `fig5_cytof_pca_timepoint.R` (uses the 16-population lineage table + MDS on
  Aitchison distance, not PCA-with-variance-explained — cannot produce the legend's
  numbers); `Fig5_manualgating_pca_age_clean.R` (same family, 1-sample exclusion instead
  of 3, PC1=38.8% — close but the biplot script's loadings are the exact textual match);
  `Fig5C_mds_panelA_feeding.R`, `Fig5D/E_olink_pca_*.R` (not described in the legend at
  all); `Fig5_manualgating_pca_biplot_no_exclusion.R`, `_filtered.R` (different
  sample/population filtering, variance numbers don't match).

### c — Change-from-baseline effect-size ranking ✅ VERIFIED
`fig5/Fig5_manualgating_synf_ctrlf_stats.R` (the `p_delta` plot object; saved as
`Fig5_manualgating_change_from_baseline.pdf`)
- Structure matches exactly: faceted by follow-up timepoint (2mo/4mo), each facet
  independently ranked by its own Cohen's d, nominal Wilcoxon p-value labeled per point,
  95% CI error bars.
- Confirmed to run on **27 populations** (own script header: "as used for the 27-population
  CyTOF figure" — cross-referenced from `Fig5_olink_synf_ctrlf_stats.R`'s header).
- Ruled out: `Fig5_manualgating_effect_rankings.R` — similar in spirit (Cohen's d
  change-from-baseline) but produces a single combined ranking colored by which timepoint
  won, not two independently-ranked faceted panels with p-values; also its header
  comment's population count ("32") is stale relative to the table it actually loads (28).

### d — Pseudotime (immune maturation) by feeding group ✅ VERIFIED
`fig5/fig5_pseudotime_manualgating_by_feeding.R` (moved here from `fig6/` earlier
2026-09-08; formerly `Fig6_manualgating_pseudotime.R`)
- Cohen's d / p-value at each timepoint: **V1 d=0.116→0.12, p=0.943 · V3 d=−0.084→−0.08,
  p=0.628 · V5 d=−0.460→−0.46, p=0.057** — matches the Results text's "d = 0.12, p = 0.943
  ... d = −0.08, p = 0.628 ... d = −0.46, p = 0.057" exactly, digit for digit.
- The parallel Olink-derived pseudotime scripts (`fig5_pseudotime_olink_by_feeding.R`,
  `fig5_pseudotime_olink_cytokine_by_feeding.R`) were not needed for this panel — the
  manuscript's Fig. 5d pseudotime is the CyTOF-derived one only.

### e — Olink volcano (change from baseline) ✅ VERIFIED (correction from earlier pass)
`fig5/Fig5_olink_synf_ctrlf_stats.R` (the `p_volcano_delta` plot object)
- Faceted by follow-up timepoint (2mo/4mo, matching the legend's "left"/"right"), x-axis
  Cohen's d change-from-baseline, y-axis −log10(nominal p), top 15 per panel labeled,
  colored by which feeding group had the higher change.
- All 16 proteins named in the Results text (KLRB1, CD160, TNFSF10, FASLG, GZMA, IL12RB1,
  CD70, CD6, LY9, LAMP3, CTSO, CCL3, CCL4, JUN, PLAUR, SCGN) appear among the top-15-labeled
  proteins across the two panels, all with `higher_in = "CtrlF"` — matching "increased less
  from baseline in SynF than in CtrlF."
- **Correction:** an earlier pass in this repo (2026-09-08, before this verification) had
  renamed `Fig5F_olink_volcano_top_hits.R` → `fig5_olink_volcano_by_feeding.R` as the
  presumed Fig. 5e source. That was wrong — it computes a different metric (raw log2FC,
  3 panels, no change-from-baseline). Archived to `_archive/superseded/`.

### WGCNA six-module detection, Module 3 eigengene, kME table
Cited in Results/Methods as "Fig. 5e/f/g" but **absent from the written legend entirely**
(the legend only describes a–e as above) — this is a genuine manuscript-legend gap to fix,
not a script problem.
- Module 3 ("brown", 40 features) eigengene by feeding group ✅ VERIFIED:
  `scripts/analysis/wgcna_module_feeding_interaction.R` — V1 d=0.249→0.25, p=0.368; V3
  d=0.047→0.05, p=0.879; V5 d=−0.436→−0.44, p=0.027 — matches the Results text exactly.
- Six-module detection: `scripts/analysis/wgcna_cytof_olink_modules_finegrained.R` (not
  re-verified numerically this pass; module/feature counts were independently corroborated
  earlier — see `projects/amplio.md` in the Brodin wiki).
- kME feature-ranking table: `fig5/Fig5_wgcna_module3_table.R` (not re-verified this pass).

## Figure 6 — Pace of immune-microbe development (FCT trajectories)
- b (trajectory fraction by feeding group) ✅ VERIFIED: `fig6/Fig6B_fct_distribution.R`
  — T1: CtrlF=18/SynF=5; T2: CtrlF=5/SynF=24; T3: CtrlF=5/SynF=9; T4: CtrlF=14/SynF=3 —
  matches the legend's "CtrlF-fed infants predominantly T1 and T4... SynF-fed infants
  predominantly T2 and T3" exactly. (This script had the stale `code/lib/` path bug,
  now fixed — it would not have run at all before today.)
- a (transition graph), c (biomarker heatmap): **no script — not yet built.** Unclear
  whether this is a Nestlé deliverable (like Fig 3/4/7) or an in-house analysis to do.
- `fig6/Fig6D_cytof_topClusters_PCA_timepoint.R`, `fig6/Fig6E_olink_topProteins_PCA_timepoint.R`
  — no FCT dependency; likely QC/diagnostic for the pseudotime derivation methodology
  rather than a numbered manuscript panel. Left in place, unclassified.

**Archived** (superseded by the v2 revision, not deleted):
`_archive/fig6_pseudotime_fct_removed_v2/` — pseudotime-by-FCT-trajectory scripts, removed
from the manuscript in the Sep2026 v2 revision.

## Figure 7 — Upper respiratory tract infections
No script in this repo. Pending: request from Nestlé.

## Supplementary Figure 3 — CyTOF gating scheme
Not audited.

## Supplementary Figure 4 — Cell populations increasing/decreasing after birth
Not audited.

## Supplementary Figure 5 — Vaccine antibody responses by feeding group
`supplementary/SupplFig5_vaccine_synf_ctrlf_v5.R` — not re-verified numerically this pass.

## Table 1 — Adverse events
No script identified. Not audited.

## Table 2 — URTI adverse events
No script identified. Contains the internal count inconsistency flagged 2026-09-08 (19 vs
41, p=0.009 in this table vs. 25 vs 46, p=0.023 in Fig. 7/Discussion) — Petter is removing
this table from the manuscript rather than reconciling it.
