# Archived: scripts identified as NOT the source of a manuscript panel

Unlike `_archive/fig6_pseudotime_fct_removed_v2/` (analysis removed from the
manuscript), these scripts were mis-identified in an earlier pass as the source of
a given panel, then ruled out once the real source was found and numerically
verified against the manuscript text (2026-09-08). Kept for reference, not deleted.

- `fig5_olink_volcano_by_feeding.R` (originally `Fig5F_olink_volcano_top_hits.R`)
  — computes raw log2FC + per-timepoint Wilcoxon p across 3 panels (V1/V3/V5).
  The manuscript's Fig. 5 Olink volcano panel instead uses Cohen's d
  change-from-baseline vs. -log10(p) across 2 panels (2mo/4mo). The correct
  script is `fig5/Fig5_olink_synf_ctrlf_stats.R` (`p_volcano_delta` object) —
  its top-15-per-panel labeled proteins reproduce all 16 proteins named in the
  manuscript's Results paragraph on the Olink attenuation pattern.
