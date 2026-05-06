suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr"))

root <- get_repo_root()
base <- load_base_tables(root)
df <- base$metadata |>
  dplyr::distinct(subject_id, group_FCT, group_feeding) |>
  dplyr::filter(!is.na(group_FCT) & !is.na(group_feeding)) |>
  dplyr::group_by(group_FCT, group_feeding) |>
  dplyr::summarise(N = dplyr::n(), .groups = "drop")

color_feeding <- c("CF" = "#39AE71", "EF" = "#33AEFA")

p <- ggplot2::ggplot(df, ggplot2::aes(x = group_FCT, y = N, fill = group_feeding)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  ggplot2::scale_fill_manual(values = color_feeding, name = "Feeding Group") +
  ggplot2::labs(x = "FCT Trajectory", y = "Number of Subjects", title = "FCT Trajectory Distribution") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig6B_fct_distribution.pdf"), plot = p, width = 6, height = 4.5, device = "pdf", dpi = 300)

