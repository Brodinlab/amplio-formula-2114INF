# Builds the analysis-ready, QC'd manual-gating table from the raw
# cytof_manual_gating_frequency.csv (itself a 1:1 mirror of the flowjo-gating-
# freqs skill's flattened output -- never modified in place):
#
#   1. Excludes the outlier sample identified 2026-09-03 (subject PHL001-0073,
#      V1, CtrlF, cytof_id 453612193) -- extreme outlier in the batch-corrected
#      MDS, confirmed not a small-plate artifact (its plate, plate_3, has 72
#      samples). Per Petter's instruction; Kanth separately double-checking.
#   2. Regresses out cytof_plate per population (linear residuals + global
#      mean, on the raw percentage scale) -- plate does NOT confound timepoint
#      or (much) feeding group (checked 2026-09-03: balanced V1/V3/V5 per
#      plate; feeding group only mildly skewed on plate_4, 30 CtrlF/42 SynF),
#      but plate_4 forms a completely separate cluster in Aitchison-distance
#      MDS, so it is a real source of extraneous variance worth removing
#      before any group or age comparison.
#
# This percentage-scale correction (not the CLR/Aitchison-space correction
# used for the MDS scripts, which is the mathematically correct approach for
# distance-based multivariate analysis) is the appropriate one for per-
# population univariate tests -- simple, directly interpretable in percentage-
# point units, standard practice for feature-by-feature batch correction.
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

OUTLIER_CYTOF_ID <- "453612193" # subject PHL001-0073, V1, CtrlF -- see header

pop_cols <- setdiff(colnames(raw), "cytof_id")

cleaned <- raw |>
  dplyr::filter(cytof_id != OUTLIER_CYTOF_ID) |>
  dplyr::left_join(metadata |> dplyr::select(cytof_id, cytof_plate), by = "cytof_id")

stopifnot(!any(is.na(cleaned$cytof_plate))) # every remaining sample must have a plate to regress out

for (pop in pop_cols) {
  fit <- stats::lm(cleaned[[pop]] ~ factor(cleaned$cytof_plate))
  cleaned[[pop]] <- stats::residuals(fit) + mean(cleaned[[pop]], na.rm = TRUE)
}

cleaned <- cleaned |> dplyr::select(cytof_id, dplyr::all_of(pop_cols))

readr::write_csv(cleaned, file.path(tables_dir, "cytof_manual_gating_frequency_clean.csv"))
cat("Wrote", nrow(cleaned), "rows (excluded 1 outlier from", nrow(raw), "), plate-corrected across", length(pop_cols), "populations\n")
