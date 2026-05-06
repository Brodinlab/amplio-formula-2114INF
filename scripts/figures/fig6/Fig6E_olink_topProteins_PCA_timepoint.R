suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "tidyr"))

root <- get_repo_root()
base <- load_base_tables(root)

olink <- base$olink_npx_wide
meta <- base$olink_metadata |>
  dplyr::mutate(
    SampleID = as.character(SampleID),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  ) |>
  dplyr::select(SampleID, timepoint)

stopifnot("SampleID" %in% colnames(olink))
x <- as.data.frame(olink)
rownames(x) <- as.character(x$SampleID)
x$SampleID <- NULL

top_prop <- 0.2
long <- x |>
  tibble::rownames_to_column("SampleID") |>
  tidyr::pivot_longer(cols = -SampleID, names_to = "Protein", values_to = "NPX") |>
  dplyr::left_join(meta, by = "SampleID") |>
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

df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(meta, by = "SampleID") |>
  dplyr::rename(PC1 = PC1, PC2 = PC2)

color_timepoint <- c("V1" = "#F2AF4AFF", "V3" = "#C36377FF", "V5" = "#1D457FFF")

p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Olink PCA colored by timepoint",
    subtitle = paste0("Top ", round(top_prop * 100), "% proteins, n = ", length(top_proteins)),
    x = "PC1", y = "PC2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), aspect.ratio = 1) +
  ggplot2::coord_equal()

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig6E_Olink_topProteins_PCA_timepoint.pdf"), plot = p, width = 4.5, height = 4.5, device = "pdf", dpi = 300)

