# Outlier check on the Olink NPX table (olink_npx_wide/olink_metadata), as
# the Olink counterpart to Fig5_manualgating_outlier_check.R. Diagnostic only
# -- no exclusion is applied here; Petter reviews the flagged samples and
# decides on a fixed cutoff (same two-step process used for the CyTOF
# manual-gating outlier check), recorded afterward as OLINK_OUTLIER_IDS in
# common.R and applied in scripts/export/export_olink_clean.R.
#
# Method: unlike the manually-gated CyTOF populations (compositional
# percentages of CD45+, analysed in CLR/Aitchison space), Olink NPX values
# are already log2-scale, non-compositional per-protein measurements -- so
# CLR does not apply here. Each protein is instead z-scored (center + unit
# variance) across samples, and per-sample distance-from-centroid in that
# scaled space is the outlier statistic, flagged via the same robust Tukey
# rule used throughout this repo (median distance + 3*MAD). There is no
# plate/batch variable tracked in olink_metadata (unlike cytof_plate for
# CyTOF), and Olink's own QC_Warning flag was already applied when the base
# table was exported (scripts/export/export_base_tables.R) -- so no
# additional batch correction is performed before computing distances.
#
# Output:
#   output/tables/olink_outlier_check.csv
#   output/figures/qc/Fig5_olink_outlier_check.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "ggplot2", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

mat <- as.data.frame(base$olink_npx_wide)
rownames(mat) <- as.character(mat$SampleID)
mat$SampleID <- NULL
mat <- as.matrix(mat)

# Proteins missing in >20% of samples are dropped before scaling/distance
# (a handful of assays are near-entirely NA; z-scoring on them would be
# dominated by imputation rather than signal). Remaining NAs are mean-imputed
# per protein so every sample gets a distance (consistent with how the
# downstream PCA scripts already zero/impute missing NPX values).
missing_frac <- colMeans(is.na(mat))
mat <- mat[, missing_frac <= 0.20, drop = FALSE]
mat <- apply(mat, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

scaled <- scale(mat) # per-protein z-score: center=0, sd=1

centroid <- colMeans(scaled)
dist_from_centroid <- sqrt(rowSums(sweep(scaled, 2, centroid, "-")^2))

med_d <- stats::median(dist_from_centroid)
mad_d <- stats::mad(dist_from_centroid)
mad_threshold <- med_d + 3 * mad_d

df <- tibble::tibble(
  SampleID = names(dist_from_centroid),
  distance = as.numeric(dist_from_centroid)
) |>
  dplyr::left_join(base$olink_metadata, by = "SampleID") |>
  dplyr::arrange(dplyr::desc(distance)) |>
  dplyr::mutate(
    flagged = distance > mad_threshold,
    status = ifelse(flagged, "Flagged (MAD rule)", "Not flagged"),
    rank = dplyr::row_number()
  )

readr::write_csv(df, file.path(root, "output", "tables", "olink_outlier_check.csv"))

status_colors <- c("Flagged (MAD rule)" = "#E8A33D", "Not flagged" = "grey40")

p <- ggplot2::ggplot(df, ggplot2::aes(x = rank, y = distance)) +
  ggplot2::geom_point(ggplot2::aes(color = status), size = 2) +
  ggplot2::geom_hline(yintercept = mad_threshold, linetype = "dashed", color = "#E8A33D") +
  ggplot2::scale_color_manual(values = status_colors, name = "Status") +
  ggplot2::labs(
    title = "Olink NPX: distance from centroid, per-protein z-scored space",
    subtitle = paste0("orange line = median+3*MAD (", round(mad_threshold, 2), "); no exclusion cutoff set yet"),
    x = "Rank (by distance, descending)", y = "Euclidean distance from centroid"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "qc", "Fig5_olink_outlier_check.pdf"), width = 6, height = 4.5)

cat("n samples:", nrow(df), "\n")
cat("n proteins retained (<=20% missing):", ncol(mat), "of", ncol(base$olink_npx_wide) - 1, "\n")
cat("MAD rule threshold:", round(mad_threshold, 3), "\n")
cat("\nFlagged by MAD rule (median+3*MAD) -- review before deciding on an exclusion cutoff:\n")
print(df |> dplyr::filter(flagged) |> dplyr::select(rank, SampleID, subject_id, timepoint, group_feeding, distance))
cat("\nTop 10 by distance (for context):\n")
print(df |> dplyr::slice_head(n = 10) |> dplyr::select(rank, SampleID, subject_id, timepoint, group_feeding, distance, status))
