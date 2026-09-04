# MDS and PCA on the full v1.1 re-gated dataset with NO outlier exclusion,
# for visual inspection before deciding what (if anything) to exclude --
# see Fig5_manualgating_outlier_check.R for the quantitative flagging (same
# preprocessing: plate-corrected CLR/Aitchison space, all 251 samples, all
# 32 populations).
#
# Points are colored by status against the two thresholds tracked in this
# repo: the fixed dist>10 cutoff Petter chose to actually exclude downstream
# (Fig5_manualgating_pca_biplot_dist10_excluded.R -- 3 samples: 453612193,
# 453611359, 453610960) vs. the 2 more that the median+3*MAD rule alone
# flags but that stay in (453612170, 453611485), shown for context.

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

# Flag outliers the same way as Fig5_manualgating_outlier_check.R, purely to
# annotate the plots below -- nothing is excluded here.
EXCLUSION_THRESHOLD <- 10 # fixed cutoff actually applied downstream
centroid <- colMeans(clr_corrected)
dist_from_centroid <- sqrt(rowSums(sweep(clr_corrected, 2, centroid, "-")^2))
mad_threshold <- stats::median(dist_from_centroid) + 3 * stats::mad(dist_from_centroid)
removed_ids <- names(dist_from_centroid)[dist_from_centroid > EXCLUSION_THRESHOLD]
flagged_kept_ids <- names(dist_from_centroid)[
  dist_from_centroid > mad_threshold & dist_from_centroid <= EXCLUSION_THRESHOLD
]

mds <- stats::cmdscale(stats::dist(clr_corrected), k = 2)
pca <- stats::prcomp(clr_corrected, center = TRUE, scale. = FALSE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

build_df <- function(coords, col1, col2) {
  as.data.frame(coords) |>
    tibble::rownames_to_column("cytof_id") |>
    dplyr::rename(AX1 = 1 + 1, AX2 = 2 + 1) |> # first two cols after rownames_to_column
    dplyr::left_join(meta_lookup, by = "cytof_id") |>
    dplyr::mutate(
      timepoint = relabel_timepoint(timepoint),
      status = dplyr::case_when(
        cytof_id %in% removed_ids ~ "Excluded (dist>10)",
        cytof_id %in% flagged_kept_ids ~ "Flagged, kept (MAD rule only)",
        TRUE ~ "Not flagged"
      )
    )
}

df_mds <- build_df(mds)
df_pca <- build_df(pca$x[, 1:2])

color_timepoint <- setNames(c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), TIMEPOINT_LABELS)
outline_colors <- c(
  "Excluded (dist>10)" = "red",
  "Flagged, kept (MAD rule only)" = "#E8A33D",
  "Not flagged" = NA
)

plot_ordination <- function(df, xlab, ylab, title) {
  highlighted <- dplyr::filter(df, status != "Not flagged")

  ggplot2::ggplot(df, ggplot2::aes(x = AX1, y = AX2)) +
    ggplot2::geom_point(ggplot2::aes(color = timepoint), size = 2, alpha = 0.7) +
    ggplot2::geom_point(
      data = highlighted,
      ggplot2::aes(shape = status), size = 4.5, stroke = 1.2, fill = NA,
      color = outline_colors[highlighted$status]
    ) +
    ggplot2::stat_ellipse(ggplot2::aes(color = timepoint), level = 0.95, alpha = 0.3) +
    ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
    ggplot2::scale_shape_manual(values = c("Excluded (dist>10)" = 21, "Flagged, kept (MAD rule only)" = 24), name = "Outlier status") +
    ggplot2::labs(
      title = title,
      subtitle = "No outlier exclusion; red circles = excluded downstream (dist>10), orange triangles = flagged but kept",
      x = xlab, y = ylab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5)
    )
}

p_mds <- plot_ordination(
  df_mds, "MDS1", "MDS2",
  "Manually-gated CyTOF (plate regressed out, NO outlier exclusion): MDS colored by age"
)
p_pca <- plot_ordination(
  df_pca,
  paste0("PC1 (", var_explained[1], "%)"), paste0("PC2 (", var_explained[2], "%)"),
  "Manually-gated CyTOF (plate regressed out, NO outlier exclusion): PCA colored by age"
)

save_pdf(p_mds, file.path(root, "output", "figures", "qc", "Fig5_manualgating_mds_age_no_outlier_exclusion.pdf"), width = 6, height = 4.8)
save_pdf(p_pca, file.path(root, "output", "figures", "qc", "Fig5_manualgating_pca_age_no_outlier_exclusion.pdf"), width = 6, height = 4.8)

cat("n samples:", nrow(df_mds), "\n")
cat("Excluded downstream (dist>10, red):", paste(removed_ids, collapse = ", "), "\n")
cat("Flagged, kept (MAD rule only, orange):", paste(flagged_kept_ids, collapse = ", "), "\n")
