# Builds the analysis-ready, QC'd manual-gating table from the raw
# cytof_manual_gating_frequency.csv (itself a 1:1 mirror of the flowjo-gating-
# freqs skill's flattened output -- never modified in place):
#
#   1. Excludes outlier samples per MANUAL_GATING_OUTLIER_IDS (common.R) --
#      final set as of 2026-09-04: 3 samples with distance-from-centroid > 10
#      in plate-corrected CLR/Aitchison space, computed on Kanth's re-gated
#      v1.1 base table (see Fig5_manualgating_outlier_check.R). Two of the
#      three (453612193, 453610960) were originally found as extreme points
#      in the plate-corrected Aitchison-distance MDS on the pre-v1.1 gating
#      (2026-09-03); the third (453611359) is new to v1.1.
#   2. Drops intermediate/complement gating-hierarchy populations per
#      MANUAL_GATING_EXCLUDED_POPULATIONS (common.R) -- NonNK, NonTB,
#      CD14neg.CD16neg, NonEosinophils (Petter's decision, 2026-09-04): not
#      informative end-populations, just parent-gate splits.
#   3. Regresses out cytof_plate per population -- plate does NOT confound
#      timepoint or (much) feeding group (checked 2026-09-03: balanced
#      V1/V3/V5 per plate; feeding group only mildly skewed on plate_4, 30
#      CtrlF/42 SynF), but plate_4 forms a completely separate cluster in
#      Aitchison-distance MDS, so it is a real source of extraneous variance
#      worth removing before any group or age comparison.
#
# CORRECTED 2026-09-03 (same day, later): the first version of this script
# regressed plate out on the RAW PERCENTAGE SCALE (residuals + global mean).
# That produced biologically impossible NEGATIVE percentages for 16 of the 32
# populations (worst: Neutrophils min -2.97%, CD8T.TEM 44 negative values) --
# an additive correction has no way to respect the [0, 100] boundary, and
# several low-abundance populations have values close enough to 0 that the
# correction routinely overshot past it. Found via a plain sanity check on
# the resulting neutrophil-to-lymphocyte ratio (Petter: "the neutrophil
# levels are low ... no?").
#
# Fixed by regressing on log(value + pseudocount) instead: residuals are
# added back on the log scale, then exponentiated, which is positive by
# construction regardless of how large the correction is. This is the same
# principle as the CLR-space correction used for the MDS scripts (Aitchison
# geometry only being valid on log-ratios), simplified to a single log
# rather than a full log-ratio, since here each population is corrected
# independently anyway (as before -- this does not attempt to keep the
# corrected populations summing to a coherent whole composition).
#
# Output: data/tables/cytof_manual_gating_frequency_clean.csv

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
})

load_required_packages(c("dplyr", "readr", "tidyr"))

root <- get_repo_root()
tables_dir <- file.path(root, "data", "tables")

raw <- readr::read_csv(file.path(tables_dir, "cytof_manual_gating_frequency.csv"), show_col_types = FALSE) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))
metadata <- readr::read_csv(file.path(tables_dir, "metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))

pop_cols <- setdiff(colnames(raw), "cytof_id")
pop_cols <- setdiff(pop_cols, MANUAL_GATING_EXCLUDED_POPULATIONS) # intermediate/complement gates, not analysed

cleaned <- raw |>
  dplyr::filter(!cytof_id %in% MANUAL_GATING_OUTLIER_IDS) |>
  dplyr::left_join(metadata |> dplyr::select(cytof_id, cytof_plate), by = "cytof_id")

stopifnot(!any(is.na(cleaned$cytof_plate))) # every remaining sample must have a plate to regress out

cleaned <- plate_correct_log_scale(cleaned, pop_cols, cleaned$cytof_plate) # common.R -- log-scale correction, see header

cleaned <- cleaned |> dplyr::select(cytof_id, dplyr::all_of(pop_cols))

readr::write_csv(cleaned, file.path(tables_dir, "cytof_manual_gating_frequency_clean.csv"))
cat("Wrote", nrow(cleaned), "rows (excluded", length(MANUAL_GATING_OUTLIER_IDS), "outliers from", nrow(raw), "), plate-corrected across", length(pop_cols), "populations\n")
