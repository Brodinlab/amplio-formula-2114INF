suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "tidyr", "viridis", "SingleCellExperiment", "slingshot"))

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
  dplyr::select(cytof_id, subject_id, timepoint, group_FCT)

cytof <- cytof |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::semi_join(meta, by = "cytof_id")

mat <- as.data.frame(cytof)
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0

# Top clusters
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

# Pseudotime: slingshot if available, else PC1
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
  pt_vals <- pca$x[, 1]
  names(pt_vals) <- rownames(pca$x)
}

pt_df <- meta |>
  dplyr::mutate(pseudotime = pt_vals[match(cytof_id, names(pt_vals))]) |>
  dplyr::filter(!is.na(pseudotime) & !is.na(subject_id) & !is.na(group_FCT) & !is.na(timepoint))

delta_df <- pt_df |>
  dplyr::group_by(subject_id, group_FCT, timepoint) |>
  dplyr::summarise(pt = mean(pseudotime, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = timepoint, values_from = pt)

df <- delta_df |>
  dplyr::filter(!is.na(V1)) |>
  dplyr::mutate(dV3 = V3 - V1, dV5 = V5 - V1) |>
  dplyr::select(subject_id, group_FCT, dV3, dV5) |>
  tidyr::pivot_longer(cols = c(dV3, dV5), names_to = "Delta_tp", values_to = "Delta_pt") |>
  dplyr::mutate(Delta_tp = dplyr::recode(Delta_tp, dV3 = "V3 - V1", dV5 = "V5 - V1"))

color_FCT <- c("T1" = "#33A190", "T2" = "#A1D5CF", "T3" = "#0F5096", "T4" = "#1C79E3")
df$group_FCT <- factor(df$group_FCT, levels = c("T1", "T2", "T3", "T4"))

p <- ggplot2::ggplot(df, ggplot2::aes(x = Delta_pt, fill = group_FCT)) +
  ggplot2::geom_density(alpha = 0.6, color = NA) +
  ggplot2::facet_grid(group_FCT ~ Delta_tp) +
  ggplot2::scale_fill_manual(values = color_FCT, name = "FCT trajectory") +
  ggplot2::labs(
    title = "CyTOF delta pseudotime distribution by FCT trajectory",
    x = "Delta pseudotime",
    y = "Density"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig6F_CyTOF_delta_pseudotime_density_faceted_FCT.pdf"), plot = p, width = 8, height = 4.5, device = "pdf", dpi = 300)

