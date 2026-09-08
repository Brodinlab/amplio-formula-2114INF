suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "tidyr", "viridis"))

root <- get_repo_root()
base <- load_base_tables(root)

olink <- base$olink_npx_wide
meta <- base$olink_metadata |>
  dplyr::mutate(
    SampleID = as.character(SampleID),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    group_FCT = factor(group_FCT, levels = c("T1", "T2", "T3", "T4"))
  )

stopifnot("SampleID" %in% colnames(olink))
x <- as.data.frame(olink)
rownames(x) <- as.character(x$SampleID)
x$SampleID <- NULL

# Top time-varying proteins by variance of timepoint means
top_prop <- 0.2
long <- x |>
  tibble::rownames_to_column("SampleID") |>
  tidyr::pivot_longer(cols = -SampleID, names_to = "Protein", values_to = "NPX") |>
  dplyr::left_join(meta |> dplyr::select(SampleID, timepoint), by = "SampleID") |>
  dplyr::filter(!is.na(NPX) & !is.na(timepoint))

prot_time_var <- long |>
  dplyr::group_by(Protein, timepoint) |>
  dplyr::summarise(mean_npx = mean(NPX, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(Protein) |>
  dplyr::summarise(time_var = stats::var(mean_npx, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(time_var))

top_proteins <- prot_time_var |>
  dplyr::slice_head(prop = top_prop) |>
  dplyr::pull(Protein)
top_proteins <- intersect(top_proteins, colnames(x))
if (length(top_proteins) < 2) stop("Not enough top time-varying Olink proteins found.")

x_top <- x[, top_proteins, drop = FALSE]
x_top[is.na(x_top)] <- 0
pca <- stats::prcomp(x_top, scale. = TRUE, center = TRUE)

# Pseudotime: slingshot if available, else PC1
pt_vals <- NULL
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("slingshot", quietly = TRUE)) {
  tryCatch({
    expr_matrix <- t(as.matrix(x_top))
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = expr_matrix, logcounts = expr_matrix),
      colData = meta[match(colnames(expr_matrix), meta$SampleID), ]
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

df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(meta, by = "SampleID") |>
  dplyr::mutate(pseudotime = pt_vals[match(SampleID, names(pt_vals))]) |>
  dplyr::rename(PC1 = PC1, PC2 = PC2)

p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = pseudotime)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::scale_color_viridis_c(name = "Pseudotime") +
  ggplot2::labs(
    title = "Olink PCA colored by pseudotime",
    subtitle = paste0("Top ", round(top_prop * 100), "% proteins, n = ", length(top_proteins)),
    x = "PC1", y = "PC2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), aspect.ratio = 1) +
  ggplot2::coord_equal()

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig6E_Olink_topProteins_PCA_pseudotime.pdf"), plot = p, width = 4.5, height = 4.5, device = "pdf", dpi = 300)

