# PCA on CLR-transformed manually-gated cell frequencies (all 32 populations,
# Kanth's re-gated v1.1 base table), plate effect regressed out first --
# identical preprocessing to Fig5_manualgating_mds_age_clean.R (same
# MANUAL_GATING_OUTLIER_IDS exclusions from common.R, same CLR +
# per-population plate regression; see that script's header for full
# rationale), so MDS and PCA are directly comparable views of the same
# corrected matrix. PCA additionally reports % variance explained per axis.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

mat <- base$cytof_manual |>
  dplyr::filter(!as.character(cytof_id) %in% MANUAL_GATING_OUTLIER_IDS) |>
  as.data.frame()
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)

clr_mat <- clr_transform(mat) # samples x populations, log-ratio space

meta_lookup <- base$metadata |>
  dplyr::distinct(cytof_id, cytof_plate, timepoint, group_feeding) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))
plate_vec <- meta_lookup$cytof_plate[match(rownames(clr_mat), meta_lookup$cytof_id)]
stopifnot(!any(is.na(plate_vec)))

# Regress out plate per population, keep residuals + global mean (preserves scale)
clr_corrected <- apply(clr_mat, 2, function(col) {
  fit <- stats::lm(col ~ factor(plate_vec))
  stats::residuals(fit) + mean(col)
})
rownames(clr_corrected) <- rownames(clr_mat)

pca <- stats::prcomp(clr_corrected, center = TRUE, scale. = FALSE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

df <- as.data.frame(pca$x[, 1:2]) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(meta_lookup, by = "cytof_id") |>
  dplyr::rename(PC1 = PC1, PC2 = PC2) |>
  dplyr::mutate(timepoint = relabel_timepoint(timepoint))

color_timepoint <- setNames(c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), TIMEPOINT_LABELS)

p_age <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::stat_ellipse(level = 0.95, alpha = 0.3) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Manually-gated CyTOF (outlier removed, plate regressed out): PCA colored by age",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

p_plate_check <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = plate_vec)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::labs(
    title = "Same, colored by plate (sanity check -- should now be mixed)",
    color = "Plate",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p_age, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_pca_age_clean.pdf"), width = 5.5, height = 4.5)
save_pdf(p_plate_check, file.path(root, "output", "figures", "qc", "Fig5_manualgating_pca_plate_clean_check.pdf"), width = 6, height = 4.5)

cat("n samples:", nrow(df), "\n")
cat("Variance explained PC1-PC5 (%):", paste(var_explained[1:5], collapse = ", "), "\n")
