suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "robCompositions"))

root <- get_repo_root()
base <- load_base_tables(root)

cytof_lineage <- base$cytof_lineage
stopifnot("cytof_id" %in% colnames(cytof_lineage))
mat <- cytof_lineage
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- mat + 0.001

aitch <- robCompositions::aDist(as.matrix(mat))
aitch[aitch < 0] <- 0
aitch <- (aitch + t(aitch)) / 2
aitch[is.na(aitch) | is.infinite(aitch) | is.nan(aitch)] <- max(aitch[is.finite(aitch)], na.rm = TRUE)

mds <- stats::cmdscale(aitch, k = 2)
df <- as.data.frame(mds) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(
    base$metadata |>
      dplyr::select(cytof_id, group_feeding),
    by = "cytof_id"
  ) |>
  dplyr::rename(MDS1 = V1, MDS2 = V2)

color_feeding <- c("CF" = "#39AE71", "EF" = "#33AEFA")

p <- ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = group_feeding)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::stat_ellipse(level = 0.95, alpha = 0.3) +
  ggplot2::scale_color_manual(values = color_feeding, name = "Feeding Group") +
  ggplot2::labs(title = "Panel A: Colored by Feeding Group", x = "MDS1", y = "MDS2") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"))

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig5C_mds_panelA_feeding.pdf"), plot = p, width = 5, height = 4, device = "pdf", dpi = 300)

