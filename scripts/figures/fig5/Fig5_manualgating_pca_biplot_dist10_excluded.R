# PCA biplot on Kanth's re-gated v1.1 manual-gating base table, excluding
# samples whose Euclidean distance from centroid (plate-corrected CLR/
# Aitchison space) exceeds 10 -- per Petter's instruction, a fixed distance
# cutoff rather than the median+3*MAD statistical rule used in
# Fig5_manualgating_outlier_check.R. At this cutoff the 3 most extreme
# samples from that check are excluded (453612193 dist=22.6, 453611359
# dist=13.8, 453610960 dist=13.6); the two borderline ones flagged by the
# MAD rule but under distance 10 (453612170 dist=7.4, 453611485 dist=6.8)
# are kept in.
#
# Batch/plate correction is applied throughout (plate_4 is a documented
# outlier batch -- see Fig5_manualgating_mds_age_clean.R header), and the
# distance-from-centroid used for exclusion, the PCA, and the loadings are
# all computed on that same batch-corrected CLR matrix.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "ggrepel"))

root <- get_repo_root()
base <- load_base_tables(root)

DIST_THRESHOLD <- 10

mat <- as.data.frame(base$cytof_manual)
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)

clr_mat <- clr_transform(mat) # samples x 32 populations, log-ratio space

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
excluded_ids <- names(dist_from_centroid)[dist_from_centroid > DIST_THRESHOLD]

clr_kept <- clr_corrected[!(rownames(clr_corrected) %in% excluded_ids), , drop = FALSE]

pca <- stats::prcomp(clr_kept, center = TRUE, scale. = FALSE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

scores <- as.data.frame(pca$x[, 1:2]) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(meta_lookup, by = "cytof_id") |>
  dplyr::mutate(timepoint = relabel_timepoint(timepoint))

loadings <- as.data.frame(pca$rotation[, 1:2])
colnames(loadings) <- c("PC1", "PC2")
loadings$population <- rownames(loadings)

score_extent <- max(abs(c(scores$PC1, scores$PC2)))
loading_extent <- max(abs(c(loadings$PC1, loadings$PC2)))
arrow_scale <- 0.8 * score_extent / loading_extent
loadings$PC1_arrow <- loadings$PC1 * arrow_scale
loadings$PC2_arrow <- loadings$PC2 * arrow_scale

color_timepoint <- setNames(c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), TIMEPOINT_LABELS)

p_biplot <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = scores, ggplot2::aes(x = PC1, y = PC2, color = timepoint), size = 2, alpha = 0.6
  ) +
  ggplot2::stat_ellipse(
    data = scores, ggplot2::aes(x = PC1, y = PC2, color = timepoint), level = 0.95, alpha = 0.3
  ) +
  ggplot2::geom_segment(
    data = loadings,
    ggplot2::aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
    arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")), color = "grey20", linewidth = 0.3
  ) +
  ggrepel::geom_text_repel(
    data = loadings, ggplot2::aes(x = PC1_arrow, y = PC2_arrow, label = population),
    size = 2.3, color = "grey20", max.overlaps = 30, segment.size = 0.2
  ) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Manually-gated CyTOF (v1.1, plate-corrected, dist>10 excluded): PCA biplot",
    subtitle = paste0("Excluded ", length(excluded_ids), " sample(s); arrows = population loadings (scaled for display)"),
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5)
  )

save_pdf(p_biplot, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_pca_biplot_dist10_excluded.pdf"), width = 7.5, height = 6)

readr::write_csv(loadings |> dplyr::select(population, PC1, PC2), file.path(root, "output", "tables", "manualgating_pca_loadings_dist10_excluded.csv"))

cat("n samples kept:", nrow(scores), "(excluded", length(excluded_ids), "of", nrow(clr_corrected), ")\n")
cat("Excluded cytof_ids (distance > 10):", paste(excluded_ids, collapse = ", "), "\n")
cat("Variance explained PC1-PC5 (%):", paste(var_explained[1:5], collapse = ", "), "\n")
cat("\nTop 5 populations by |PC1 loading|:\n")
print(loadings |> dplyr::arrange(dplyr::desc(abs(PC1))) |> dplyr::slice_head(n = 5) |> dplyr::select(population, PC1, PC2))
cat("\nTop 5 populations by |PC2 loading|:\n")
print(loadings |> dplyr::arrange(dplyr::desc(abs(PC2))) |> dplyr::slice_head(n = 5) |> dplyr::select(population, PC1, PC2))
