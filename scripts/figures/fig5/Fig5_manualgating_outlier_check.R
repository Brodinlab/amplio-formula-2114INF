# Outlier check on Kanth's re-gated v1.1 manual-gating base table.
#
# The previous outlier exclusions (subject PHL001-0073 / cytof_id 453612193,
# subject PHL022-0012 / cytof_id 453610960; see
# scripts/export/export_cytof_manual_clean.R) were found ad hoc on the OLD
# gating scheme, before the v1.1 structural re-gate (2026-09-04). Per
# Petter's instruction, re-checking from scratch on the new data rather than
# assuming those two still apply -- no exclusion is applied here.
#
# Method: same plate-corrected Aitchison (CLR) space used throughout this
# repo's ordination scripts. Per-sample distance-from-centroid in that space
# is the outlier statistic; flagged via Tukey's rule (median distance +
# 3*MAD -- MAD instead of SD to stay robust to the outliers being sought).
#
# Per Petter's decision (2026-09-04): the actual exclusion applied downstream
# (Fig5_manualgating_pca_biplot_dist10_excluded.R) uses a fixed distance>10
# cutoff, not the MAD rule -- 3 of the 5 MAD-flagged samples cross it
# (453612193, 453611359, 453610960); the other 2 (453612170, 453611485) are
# flagged but kept. Both are shown here for context.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

mat <- as.data.frame(base$cytof_manual)
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)

clr_mat <- clr_transform(mat)

meta_lookup <- base$metadata |>
  dplyr::distinct(subject_id, cytof_id, cytof_plate, timepoint, group_feeding) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))
plate_vec <- meta_lookup$cytof_plate[match(rownames(clr_mat), meta_lookup$cytof_id)]
stopifnot(!any(is.na(plate_vec)))

clr_corrected <- apply(clr_mat, 2, function(col) {
  fit <- stats::lm(col ~ factor(plate_vec))
  stats::residuals(fit) + mean(col)
})
rownames(clr_corrected) <- rownames(clr_mat)

centroid <- colMeans(clr_corrected)
dist_from_centroid <- sqrt(rowSums(sweep(clr_corrected, 2, centroid, "-")^2))

med_d <- stats::median(dist_from_centroid)
mad_d <- stats::mad(dist_from_centroid)
mad_threshold <- med_d + 3 * mad_d
EXCLUSION_THRESHOLD <- 10 # fixed cutoff actually applied downstream (Petter's decision, 2026-09-04)

df <- tibble::tibble(
  cytof_id = names(dist_from_centroid),
  distance = as.numeric(dist_from_centroid)
) |>
  dplyr::left_join(meta_lookup, by = "cytof_id") |>
  dplyr::arrange(dplyr::desc(distance)) |>
  dplyr::mutate(
    flagged = distance > mad_threshold,
    excluded = distance > EXCLUSION_THRESHOLD,
    status = dplyr::case_when(
      excluded ~ "Excluded (dist>10)",
      flagged ~ "Flagged, kept (MAD rule only)",
      TRUE ~ "Not flagged"
    ),
    rank = dplyr::row_number()
  )

OLD_OUTLIER_IDS <- c("453612193", "453610960")
df <- df |> dplyr::mutate(was_old_outlier = cytof_id %in% OLD_OUTLIER_IDS)

readr::write_csv(df, file.path(root, "output", "tables", "manualgating_outlier_check.csv"))

status_colors <- c(
  "Excluded (dist>10)" = "red",
  "Flagged, kept (MAD rule only)" = "#E8A33D",
  "Not flagged" = "grey40"
)

p <- ggplot2::ggplot(df, ggplot2::aes(x = rank, y = distance)) +
  ggplot2::geom_point(ggplot2::aes(color = status), size = 2) +
  ggplot2::geom_hline(yintercept = mad_threshold, linetype = "dashed", color = "#E8A33D") +
  ggplot2::geom_hline(yintercept = EXCLUSION_THRESHOLD, linetype = "dashed", color = "red") +
  ggplot2::scale_color_manual(values = status_colors, name = "Status") +
  ggplot2::labs(
    title = "Distance from centroid, plate-corrected CLR space (v1.1 re-gated data)",
    subtitle = paste0(
      "red line = exclusion cutoff (dist>", EXCLUSION_THRESHOLD,
      "); orange line = median+3*MAD (", round(mad_threshold, 2), "), for reference only"
    ),
    x = "Rank (by distance, descending)", y = "Euclidean distance from centroid"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "qc", "Fig5_manualgating_outlier_check_v1.1.pdf"), width = 6, height = 4.5)

cat("n samples:", nrow(df), "\n")
cat("Exclusion cutoff (dist>10):", sum(df$excluded), "sample(s) excluded\n")
cat("MAD rule threshold (context only):", round(mad_threshold, 3), "\n")
cat("\nExcluded (dist>10):\n")
print(df |> dplyr::filter(excluded) |> dplyr::select(rank, cytof_id, subject_id, timepoint, distance, was_old_outlier))
cat("\nFlagged by MAD rule but kept (dist<=10):\n")
print(df |> dplyr::filter(flagged & !excluded) |> dplyr::select(rank, cytof_id, subject_id, timepoint, distance, was_old_outlier))
cat("\nTop 10 by distance (for context):\n")
print(df |> dplyr::slice_head(n = 10) |> dplyr::select(rank, cytof_id, subject_id, timepoint, distance, status, was_old_outlier))
cat("\nOld outlier IDs' new rank/distance:\n")
print(df |> dplyr::filter(was_old_outlier) |> dplyr::select(rank, cytof_id, subject_id, distance, status))
