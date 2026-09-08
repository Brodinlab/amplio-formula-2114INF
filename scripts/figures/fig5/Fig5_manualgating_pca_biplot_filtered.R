# PCA + biplot on the manually-gated CyTOF table, restricted to a filtered
# subset of the 32 populations, on top of the standard preprocessing used
# throughout this repo's ordination scripts (MANUAL_GATING_OUTLIER_IDS
# excluded, plate regressed out per population in CLR/Aitchison space --
# see Fig5_manualgating_pca_biplot_dist10_excluded.R).
#
# Two independent filtering criteria (Petter's request: "overall variance or
# time-associated changes") -- both run, not just one, so they can be
# compared directly:
#   1. Overall variance: populations ranked by variance of their
#      plate-corrected CLR values, pooled across all samples/timepoints/
#      groups (group- and time-blind).
#   2. Time-associated: populations ranked by eta-squared (one-way ANOVA
#      vs. timepoint) on the same CLR values -- same method used to select
#      the pseudotime input populations in fig5_pseudotime_manualgating_by_feeding.R,
#      reused here for consistency.
# Both take the top 20 of 32 populations (Petter: 6/32 was too skewed).
# Change N_TOP below to adjust.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "ggrepel", "purrr"))

root <- get_repo_root()
base <- load_base_tables(root)

N_TOP <- 20

mat <- as.data.frame(base$cytof_manual) |>
  dplyr::filter(!as.character(cytof_id) %in% MANUAL_GATING_OUTLIER_IDS)
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

timepoint_vec <- meta_lookup$timepoint[match(rownames(clr_corrected), meta_lookup$cytof_id)]

pop_cols <- colnames(clr_corrected)
n_top <- min(N_TOP, length(pop_cols))

# ---- Criterion 1: overall variance (pooled, group/time-blind) ----
overall_variance <- apply(clr_corrected, 2, stats::var)
top_by_variance <- names(sort(overall_variance, decreasing = TRUE))[seq_len(n_top)]

# ---- Criterion 2: time-association (eta-squared vs. timepoint) ----
eta_squared <- purrr::map_dbl(pop_cols, function(pop) {
  d <- data.frame(value = clr_corrected[, pop], timepoint = timepoint_vec)
  fit <- aov(value ~ timepoint, data = d)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
})
names(eta_squared) <- pop_cols
top_by_time <- names(sort(eta_squared, decreasing = TRUE))[seq_len(n_top)]

cat("Top", n_top, "populations by overall variance:\n")
print(round(sort(overall_variance, decreasing = TRUE)[seq_len(n_top)], 3))
cat("\nTop", n_top, "populations by time-association (eta^2 vs timepoint):\n")
print(round(sort(eta_squared, decreasing = TRUE)[seq_len(n_top)], 3))

readr::write_csv(
  tibble::tibble(population = pop_cols, overall_variance = overall_variance, eta_squared_timepoint = eta_squared) |>
    dplyr::arrange(dplyr::desc(overall_variance)),
  file.path(root, "output", "tables", "manualgating_pca_population_filter_rankings.csv")
)

# ---- Shared plotting machinery ----
color_timepoint <- setNames(c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), TIMEPOINT_LABELS)

make_biplot <- function(selected_populations, title_suffix) {
  sub_mat <- clr_corrected[, selected_populations, drop = FALSE]
  pca <- stats::prcomp(sub_mat, center = TRUE, scale. = FALSE)
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

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = scores, ggplot2::aes(x = PC1, y = PC2, color = timepoint), size = 2, alpha = 0.6) +
    ggplot2::stat_ellipse(data = scores, ggplot2::aes(x = PC1, y = PC2, color = timepoint), level = 0.95, alpha = 0.3) +
    ggplot2::geom_segment(
      data = loadings, ggplot2::aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")), color = "grey20", linewidth = 0.3
    ) +
    ggrepel::geom_text_repel(
      data = loadings, ggplot2::aes(x = PC1_arrow, y = PC2_arrow, label = population),
      size = 2.6, color = "grey20", max.overlaps = 30, segment.size = 0.2
    ) +
    ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
    ggplot2::labs(
      title = paste0("Manually-gated CyTOF PCA biplot: ", title_suffix),
      subtitle = paste0(length(selected_populations), " of 32 populations (top ", N_TOP, ")"),
      x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5)
    )

  list(plot = p, loadings = loadings, var_explained = var_explained)
}

result_variance <- make_biplot(top_by_variance, "top populations by overall variance")
result_time <- make_biplot(top_by_time, "top populations by time-association")

save_pdf(result_variance$plot, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_pca_biplot_filtered_variance.pdf"), width = 7, height = 6)
save_pdf(result_time$plot, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_pca_biplot_filtered_timeassoc.pdf"), width = 7, height = 6)

cat("\n[Overall-variance filter] n samples:", nrow(base$cytof_manual) - length(MANUAL_GATING_OUTLIER_IDS),
    "| variance explained PC1-PC2:", paste(result_variance$var_explained[1:2], collapse = ", "), "%\n")
cat("[Time-association filter] n samples:", nrow(base$cytof_manual) - length(MANUAL_GATING_OUTLIER_IDS),
    "| variance explained PC1-PC2:", paste(result_time$var_explained[1:2], collapse = ", "), "%\n")
