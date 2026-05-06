suppressPackageStartupMessages({
  source("code/lib/common.R")
})

load_required_packages(c("ggplot2", "dplyr", "readr"))

root <- get_repo_root()
in_path <- file.path(root, "data/figure_data/fig2/fig2E_igsq13_gi_tolerance.csv")
df <- readr::read_csv(in_path, show_col_types = FALSE)

age_months <- c(0.75, 1.5, 2, 3, 4, 6)
visit_labels <- c("V1", "V2", "V3", "V4", "V5", "V6")

p <- ggplot2::ggplot(df, ggplot2::aes(x = age_months, y = mean, color = group, group = group)) +
  ggplot2::geom_hline(yintercept = 30, linetype = "dashed", color = "red", alpha = 0.7) +
  ggplot2::geom_hline(yintercept = 23, linetype = "dashed", color = "orange", alpha = 0.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.8) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c("EF" = "#4A90D9", "CF" = "#5CB85C")) +
  ggplot2::scale_x_continuous(breaks = age_months, labels = visit_labels, limits = c(0, 6.5)) +
  ggplot2::scale_y_continuous(limits = c(13, 32), breaks = seq(13, 31, 2)) +
  ggplot2::labs(
    title = "Mean IGSQ-13 Scores of GI Tolerance",
    x = "Visit",
    y = "IGSQ-13 Score (95% CI)",
    color = ""
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    aspect.ratio = 1,
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.background = ggplot2::element_rect(fill = "white", color = "gray80"),
    axis.text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 12)
  ) +
  ggplot2::annotate(
    "text", x = 0.5, y = 30.5,
    label = "Clinically meaningful GI discomfort",
    size = 3, hjust = 0, color = "red"
  ) +
  ggplot2::annotate(
    "text", x = 0.5, y = 23.5,
    label = "Potential GI discomfort",
    size = 3, hjust = 0, color = "orange"
  )

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig2E_igsq13_gi_tolerance.pdf"), plot = p, width = 8, height = 6)

