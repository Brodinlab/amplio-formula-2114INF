# Olink NPX: PCA biplot with loading arrows for the top 20 age-associated
# proteins (by eta-squared vs timepoint, Fig5_olink_pca_age.R), in the same
# biplot style as Fig5_manualgating_pca_biplot_no_exclusion.R. Uses the full
# all-protein PCA (so PC1/PC2 reflect overall structure, not just the top
# proteins), but only labels/draws arrows for the top-ranked subset -- with
# 358 proteins, labeling all of them (as the 32-population CyTOF biplot does)
# would be unreadable.
#
# Output:
#   output/figures/manuscript/Fig5_olink_top_age_biplot.pdf
#   output/tables/olink_top_age_loadings.csv

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "ggrepel", "readr"))

root <- get_repo_root()
base <- load_base_tables(root)

N_TOP <- 20

meta <- base$olink_metadata |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("V1", "V3", "V5")))

mat <- as.data.frame(base$olink_npx_clean)
rownames(mat) <- mat$SampleID
mat$SampleID <- NULL

age_association <- readr::read_csv(file.path(root, "output", "tables", "olink_age_association.csv"), show_col_types = FALSE)
top_proteins <- age_association |> dplyr::arrange(dplyr::desc(eta_squared)) |> dplyr::slice_head(n = N_TOP) |> dplyr::pull(protein)

pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

scores <- as.data.frame(pca$x[, 1:2]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(meta |> dplyr::select(SampleID, timepoint), by = "SampleID")

loadings_all <- as.data.frame(pca$rotation[, 1:2])
colnames(loadings_all) <- c("PC1", "PC2")
loadings_all$protein <- rownames(loadings_all)
loadings <- loadings_all |> dplyr::filter(protein %in% top_proteins)

# Arrows scaled to occupy ~80% of the score cloud's range -- same convention
# as the CyTOF biplot -- computed from the top-15 subset only (scaling to the
# full 358-protein loading range would shrink these top hits' arrows to
# near-invisible, since a handful of low-ranked proteins can have much larger
# raw loadings on other axes).
score_extent <- max(abs(c(scores$PC1, scores$PC2)))
loading_extent <- max(abs(c(loadings$PC1, loadings$PC2)))
arrow_scale <- 0.8 * score_extent / loading_extent
loadings$PC1_arrow <- loadings$PC1 * arrow_scale
loadings$PC2_arrow <- loadings$PC2 * arrow_scale

color_timepoint <- c(V1 = "#F2AF4AFF", V3 = "#C36377FF", V5 = "#1D457FFF")

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
    data = loadings, ggplot2::aes(x = PC1_arrow, y = PC2_arrow, label = protein),
    size = 2.6, color = "grey20", max.overlaps = 30, segment.size = 0.2
  ) +
  ggplot2::scale_color_manual(values = color_timepoint, labels = TIMEPOINT_LABELS[c("V1", "V3", "V5")], name = "Timepoint") +
  ggplot2::labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  ggplot2::theme_bw(base_size = 10)

save_pdf(p_biplot, file.path(root, "output", "figures", "manuscript", "Fig5_olink_top_age_biplot.pdf"), width = 7, height = 6)

readr::write_csv(
  loadings_all |> dplyr::filter(protein %in% top_proteins) |>
    dplyr::left_join(age_association, by = "protein") |>
    dplyr::arrange(dplyr::desc(eta_squared)) |>
    dplyr::select(protein, eta_squared, p_value, p_fdr, PC1, PC2),
  file.path(root, "output", "tables", "olink_top_age_loadings.csv")
)

cat("Top", N_TOP, "age-associated proteins (loadings on PC1/PC2):\n")
print(loadings |> dplyr::select(protein, PC1, PC2) |> dplyr::arrange(dplyr::desc(abs(PC1))))
