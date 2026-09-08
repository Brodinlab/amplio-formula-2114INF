suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "tidyr", "viridis"))

root <- get_repo_root()
base <- load_base_tables(root)

cytof <- base$cytof_cluster
stopifnot("cytof_id" %in% colnames(cytof))

meta <- base$metadata |>
  dplyr::filter(!is.na(cytof_id)) |>
  dplyr::mutate(
    cytof_id = as.character(cytof_id),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    group_FCT = factor(group_FCT, levels = c("T1", "T2", "T3", "T4"))
  ) |>
  dplyr::select(cytof_id, subject_id, timepoint, group_FCT, group_feeding)

# Align cytof rows to metadata
cytof <- cytof |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::semi_join(meta, by = "cytof_id")

mat <- as.data.frame(cytof)
rownames(mat) <- mat$cytof_id
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0

# Top time-varying clusters (same idea as amplio_merge extra section)
top_prop <- 0.2
cluster_long <- mat |>
  tibble::rownames_to_column("cytof_id") |>
  tidyr::pivot_longer(cols = -cytof_id, names_to = "Cluster", values_to = "Frequency") |>
  dplyr::left_join(meta, by = "cytof_id") |>
  dplyr::filter(!is.na(Frequency) & !is.na(timepoint))

time_var <- cluster_long |>
  dplyr::group_by(Cluster, timepoint) |>
  dplyr::summarise(mean_freq = mean(Frequency, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::summarise(time_var = stats::var(mean_freq, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(time_var))

top_clusters <- time_var |>
  dplyr::slice_head(prop = top_prop) |>
  dplyr::pull(Cluster)
top_clusters <- intersect(top_clusters, colnames(mat))
if (length(top_clusters) < 2) stop("Not enough top time-varying CyTOF clusters found.")

x <- mat[, top_clusters, drop = FALSE]
x_clr <- clr_transform(as.matrix(x))

pca <- stats::prcomp(x_clr, scale. = TRUE, center = TRUE)

pt_vals <- NULL
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("slingshot", quietly = TRUE)) {
  tryCatch({
    expr_matrix <- t(x_clr)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = expr_matrix, logcounts = expr_matrix),
      colData = meta[match(colnames(expr_matrix), meta$cytof_id), ]
    )
    n_pcs <- min(10, ncol(pca$x))
    SingleCellExperiment::reducedDim(sce, "PCA") <- pca$x[, seq_len(n_pcs), drop = FALSE]
    SingleCellExperiment::colData(sce)$timepoint <- factor(SingleCellExperiment::colData(sce)$timepoint, levels = c("V1", "V3", "V5"))
    sce <- slingshot::slingshot(sce, clusterLabels = "timepoint", reducedDim = "PCA", start.clus = "V1")
    pt <- slingshot::slingPseudotime(sce)
    pt_vals <- if (is.matrix(pt)) pt[, 1] else pt
  }, error = function(e) {
    pt_vals <<- NULL
  })
}
if (is.null(pt_vals)) {
  # Deterministic fallback
  pt_vals <- pca$x[, 1]
  names(pt_vals) <- rownames(pca$x)
}

df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(meta, by = "cytof_id") |>
  dplyr::mutate(pseudotime = pt_vals[match(cytof_id, names(pt_vals))]) |>
  dplyr::rename(PC1 = PC1, PC2 = PC2)

p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = pseudotime)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::scale_color_viridis_c(name = "Pseudotime") +
  ggplot2::labs(
    title = "CyTOF PCA colored by pseudotime",
    subtitle = paste0("Top ", round(top_prop * 100), "% clusters, n = ", length(top_clusters)),
    x = "PC1", y = "PC2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), aspect.ratio = 1) +
  ggplot2::coord_equal()

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig6D_CyTOF_topClusters_PCA_pseudotime.pdf"), plot = p, width = 4.5, height = 4.5, device = "pdf", dpi = 300)

