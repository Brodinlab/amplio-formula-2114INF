# Builds the analysis-ready, QC'd Olink NPX table from the raw
# olink_npx_wide.csv (itself already filtered to QC_Warning != "EXCLUDED" at
# export time, see export_base_tables.R -- never modified in place):
#
#   1. Excludes outlier samples per OLINK_OUTLIER_IDS (common.R) -- currently
#      empty; Petter reviewed the outlier-check diagnostic
#      (Fig5_olink_outlier_check.R) and decided to keep all samples, since
#      distances declined smoothly with no clear separation (2026-09-06).
#   2. Drops proteins missing in >20% of samples (a handful of assays are
#      near-entirely NA), then mean-imputes any remaining missing values per
#      protein -- same rule used in Fig5_olink_outlier_check.R, applied here
#      as the single source of truth so every downstream script sees the
#      same protein set and imputed values.
#
# No plate/batch correction is applied: olink_metadata does not track a
# batch variable (unlike cytof_plate for CyTOF), so there is nothing to
# regress out.
#
# Output: data/tables/olink_npx_clean.csv

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
})

load_required_packages(c("dplyr", "readr"))

root <- get_repo_root()
tables_dir <- file.path(root, "data", "tables")

raw <- readr::read_csv(file.path(tables_dir, "olink_npx_wide.csv"), show_col_types = FALSE) |>
  dplyr::mutate(SampleID = as.character(SampleID)) |>
  dplyr::filter(!SampleID %in% OLINK_OUTLIER_IDS)

mat <- as.data.frame(raw)
rownames(mat) <- mat$SampleID
mat$SampleID <- NULL

missing_frac <- colMeans(is.na(mat))
mat <- mat[, missing_frac <= 0.20, drop = FALSE]
mat <- apply(mat, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

cleaned <- as.data.frame(mat)
cleaned$SampleID <- rownames(mat)
cleaned <- cleaned |> dplyr::select(SampleID, dplyr::everything())

readr::write_csv(cleaned, file.path(tables_dir, "olink_npx_clean.csv"))
cat(
  "Wrote", nrow(cleaned), "rows (excluded", length(OLINK_OUTLIER_IDS), "outliers from", nrow(raw), "),",
  ncol(cleaned) - 1, "proteins retained (of", ncol(raw) - 1, "raw)\n"
)
