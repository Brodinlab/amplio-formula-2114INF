# Visualizes the 40 features making up fine-grained WGCNA Module 3 ("brown":
# CD8 T central-memory, CD8 T effector-memory, CD8 T TEMRA + 37 Olink
# proteins) and highlights how tightly co-regulated they are relative to the
# rest of the network -- per Petter's request, 2026-09-08.
#
# Two panels:
#   A. Pairwise biweight-midcorrelation (bicor) heatmap of the 40 Module 3
#      features, hierarchically clustered -- the same correlation function
#      used to build the WGCNA network itself (adjacency(..., corFnc =
#      "bicor")), so this is a direct, un-thresholded view of the same
#      co-regulation structure the module detection algorithm acted on.
#   B. Density comparison of pairwise bicor values: all 780 within-Module-3
#      pairs vs. all pairs among the remaining (background) features. A
#      module built from real co-regulation should sit well to the right of
#      (more positive than) the background distribution, which should center
#      near zero.
#
# No formal significance test is run on panel B -- pairwise correlations
# among a fixed set of samples are not independent observations, so a
# Wilcoxon/t-test p-value across ~780 vs ~59,000 non-independent pairs would
# be misleadingly small. Medians are reported descriptively instead.
#
# Output:
#   output/figures/manuscript/wgcna_module3_feature_heatmap.pdf
#   output/figures/manuscript/wgcna_module3_vs_background_density.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
  library(WGCNA)
})
WGCNA::disableWGCNAThreads()

load_required_packages(c("dplyr", "tidyr", "tibble", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# ---- Rebuild the same combined CyTOF(CLR) + Olink(NPX) matrix used throughout ----
cytof_meta <- base$metadata |>
  dplyr::distinct(cytof_id, subject_id, timepoint) |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::filter(!is.na(cytof_id))

cytof_mat <- as.data.frame(base$cytof_manual_clean)
rownames(cytof_mat) <- cytof_mat$cytof_id
cytof_mat$cytof_id <- NULL
cytof_clr <- clr_transform(as.matrix(cytof_mat))
colnames(cytof_clr) <- paste0("CyTOF_", colnames(cytof_clr))

olink_mat <- as.data.frame(base$olink_npx_clean)
rownames(olink_mat) <- olink_mat$SampleID
olink_mat$SampleID <- NULL
colnames(olink_mat) <- paste0("Olink_", colnames(olink_mat))

link <- cytof_meta |>
  dplyr::inner_join(
    base$olink_metadata |> dplyr::distinct(SampleID, subject_id, timepoint),
    by = c("subject_id", "timepoint")
  ) |>
  dplyr::filter(cytof_id %in% rownames(cytof_clr), SampleID %in% rownames(olink_mat))

combined <- cbind(
  cytof_clr[link$cytof_id, , drop = FALSE],
  as.matrix(olink_mat[link$SampleID, , drop = FALSE])
)

# ---- Module 3 ("brown") feature set + numeric label ----
assignment <- readr::read_csv(
  file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"),
  show_col_types = FALSE
)
module_sizes <- assignment |> dplyr::filter(module != "grey") |> dplyr::count(module, sort = TRUE)
brown_number <- which(module_sizes$module == "brown")

module_features <- assignment |> dplyr::filter(module == "brown") |> dplyr::pull(feature)
background_features <- assignment |> dplyr::filter(module != "brown") |> dplyr::pull(feature)
feature_names <- setNames(assignment$feature_name, assignment$feature)

cat("Module", brown_number, "(\"brown\"):", length(module_features), "features.",
    "Background:", length(background_features), "features.\n")

# ---- Pairwise bicor, same correlation function used to build the network ----
bicor_mat <- WGCNA::bicor(combined)

module_cor <- bicor_mat[module_features, module_features]
dimnames(module_cor) <- list(feature_names[module_features], feature_names[module_features])

# ---- Panel A: Module 3 heatmap, hierarchically clustered ----
hc <- stats::hclust(stats::as.dist(1 - module_cor), method = "average")
ordered_features <- rownames(module_cor)[hc$order]

heatmap_df <- module_cor |>
  as.data.frame() |>
  tibble::rownames_to_column("feature_a") |>
  tidyr::pivot_longer(-feature_a, names_to = "feature_b", values_to = "bicor") |>
  dplyr::mutate(
    feature_a = factor(feature_a, levels = ordered_features),
    feature_b = factor(feature_b, levels = ordered_features)
  )

p_heatmap <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = feature_a, y = feature_b, fill = bicor)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "bicor") +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = NULL, y = NULL, title = paste0("Module ", brown_number, " (", length(module_features), " features): pairwise correlation")) +
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    axis.text.y = ggplot2::element_text(size = 5.5),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = 9, hjust = 0.5)
  )

save_pdf(p_heatmap, file.path(root, "output", "figures", "manuscript", "wgcna_module3_feature_heatmap.pdf"), width = 8, height = 8)

# ---- Panel B: within-module vs. background pairwise bicor ----
within_module <- module_cor[upper.tri(module_cor)]

background_cor <- bicor_mat[background_features, background_features]
background_vals <- background_cor[upper.tri(background_cor)]

cat("\nWithin-Module", brown_number, "pairwise bicor: n =", length(within_module),
    ", median =", round(median(within_module), 3), ", mean =", round(mean(within_module), 3), "\n")
cat("Background pairwise bicor: n =", length(background_vals),
    ", median =", round(median(background_vals), 3), ", mean =", round(mean(background_vals), 3), "\n")

comparison_df <- dplyr::bind_rows(
  tibble::tibble(bicor = within_module, set = paste0("Module ", brown_number)),
  tibble::tibble(bicor = background_vals, set = "Background")
) |>
  dplyr::mutate(set = factor(set, levels = c(paste0("Module ", brown_number), "Background")))

medians_df <- comparison_df |> dplyr::group_by(set) |> dplyr::summarise(median_bicor = median(bicor), .groups = "drop")

color_set <- setNames(c("#D2691E", "grey60"), c(paste0("Module ", brown_number), "Background"))

p_density <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = bicor, fill = set, color = set)) +
  ggplot2::geom_density(alpha = 0.4, linewidth = 0.6) +
  ggplot2::geom_vline(data = medians_df, ggplot2::aes(xintercept = median_bicor, color = set), linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = color_set, guide = "none") +
  ggplot2::scale_color_manual(values = color_set, guide = "none") +
  ggplot2::annotate("text",
    x = medians_df$median_bicor[medians_df$set == paste0("Module ", brown_number)], y = Inf,
    label = paste0("Module ", brown_number, "\n(median r=", round(medians_df$median_bicor[medians_df$set == paste0("Module ", brown_number)], 2), ")"),
    color = color_set[paste0("Module ", brown_number)], vjust = 1.3, size = 3, fontface = "bold"
  ) +
  ggplot2::annotate("text",
    x = medians_df$median_bicor[medians_df$set == "Background"], y = Inf,
    label = paste0("Background\n(median r=", round(medians_df$median_bicor[medians_df$set == "Background"], 2), ")"),
    color = "grey40", vjust = 3.2, size = 3, fontface = "bold"
  ) +
  ggplot2::labs(x = "Pairwise correlation (bicor)", y = "Density") +
  ggplot2::theme_classic(base_size = 9)

save_pdf(p_density, file.path(root, "output", "figures", "manuscript", "wgcna_module3_vs_background_density.pdf"), width = 6.5, height = 5)

cat("\nSaved:\n  output/figures/manuscript/wgcna_module3_feature_heatmap.pdf\n  output/figures/manuscript/wgcna_module3_vs_background_density.pdf\n")
