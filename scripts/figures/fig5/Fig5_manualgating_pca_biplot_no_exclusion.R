# PCA biplot on Kanth's re-gated v1.1 manual-gating base table -- all 251
# samples, NO outlier exclusion (per Petter's instruction: outliers found
# so far are not dramatic enough to exclude outright; revisit later).
# Batch/plate correction IS applied -- plate_4 is a clear, well-documented
# outlier batch in this dataset (see Fig5_manualgating_mds_age_clean.R
# header), regressed out per-population in CLR (Aitchison) space, same as
# every other ordination script in this repo. All downstream PCA is on that
# batch-corrected CLR matrix. Adds loading arrows (all 32 populations) to
# show which populations drive each axis.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "ggrepel"))

root <- get_repo_root()
base <- load_base_tables(root)

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

# Batch (plate) correction, per population, in CLR space -- plate_4 is the
# documented outlier batch this removes.
clr_corrected <- apply(clr_mat, 2, function(col) {
  fit <- stats::lm(col ~ factor(plate_vec))
  stats::residuals(fit) + mean(col)
})
rownames(clr_corrected) <- rownames(clr_mat)

pca <- stats::prcomp(clr_corrected, center = TRUE, scale. = FALSE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

scores <- as.data.frame(pca$x[, 1:2]) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(meta_lookup, by = "cytof_id") |>
  dplyr::mutate(timepoint = relabel_timepoint(timepoint))

# Loadings (rotation matrix), scaled to occupy ~80% of the score cloud's
# range -- standard biplot convention so arrows and points share one plot
# legibly (loadings and scores otherwise live on very different scales).
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
    title = "Manually-gated CyTOF (v1.1, plate-corrected, NO outlier exclusion): PCA biplot",
    subtitle = "Arrows: population loadings (scaled for display)",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5)
  )

save_pdf(p_biplot, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_pca_biplot_no_exclusion.pdf"), width = 7.5, height = 6)

readr::write_csv(loadings |> dplyr::select(population, PC1, PC2), file.path(root, "output", "tables", "manualgating_pca_loadings_no_exclusion.csv"))

cat("n samples:", nrow(scores), "\n")
cat("Variance explained PC1-PC5 (%):", paste(var_explained[1:5], collapse = ", "), "\n")
cat("\nTop 5 populations by |PC1 loading|:\n")
print(loadings |> dplyr::arrange(dplyr::desc(abs(PC1))) |> dplyr::slice_head(n = 5) |> dplyr::select(population, PC1, PC2))
cat("\nTop 5 populations by |PC2 loading|:\n")
print(loadings |> dplyr::arrange(dplyr::desc(abs(PC2))) |> dplyr::slice_head(n = 5) |> dplyr::select(population, PC1, PC2))
