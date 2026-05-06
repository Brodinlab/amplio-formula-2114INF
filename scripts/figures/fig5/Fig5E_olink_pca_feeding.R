suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

olink <- base$olink_npx_wide
stopifnot("SampleID" %in% colnames(olink))
x <- as.data.frame(olink)
rownames(x) <- x$SampleID
x$SampleID <- NULL
x[is.na(x)] <- 0

pca <- stats::prcomp(x, scale. = TRUE, center = TRUE)
df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(base$olink_metadata, by = "SampleID") |>
  dplyr::rename(PC1 = PC1, PC2 = PC2)

color_feeding <- c("CF" = "#39AE71", "EF" = "#33AEFA")

p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = group_feeding)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::scale_color_manual(values = color_feeding, name = "Feeding group") +
  ggplot2::labs(title = "Olink PCA (colored by feeding group)", x = "PC1", y = "PC2") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
  ggplot2::coord_equal()

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig5E_olink_pca_feeding.pdf"), plot = p, width = 6, height = 5, device = "pdf", dpi = 300)

