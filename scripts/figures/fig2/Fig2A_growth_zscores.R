suppressPackageStartupMessages({
  source("code/lib/common.R")
})

load_required_packages(c("ggplot2", "dplyr", "readr", "gridExtra"))

root <- get_repo_root()
in_path <- file.path(root, "data/figure_data/fig2/fig2A_growth_zscores.csv")
df <- readr::read_csv(in_path, show_col_types = FALSE)

age_months <- c(0.75, 1.5, 2, 3, 4, 6)
visit_labels <- c("V1", "V2", "V3", "V4", "V5", "V6")

common_theme <- ggplot2::theme_classic() +
  ggplot2::theme(
    aspect.ratio = 1,
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = c(0.15, 0.95),
    legend.background = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.4, "cm"),
    legend.text = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 11),
    axis.title.y = ggplot2::element_text(size = 10)
  )

create_plot <- function(d, title, ylab, ylim_range, ref_line = 0) {
  ggplot2::ggplot(d, ggplot2::aes(x = age_months, y = mean, color = group, group = group)) +
    ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed", color = "gray50", alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15, linewidth = 0.6
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c("EF" = "#4A90D9", "CF" = "#5CB85C")) +
    ggplot2::scale_x_continuous(
      breaks = age_months,
      labels = visit_labels,
      limits = c(0, 6.5)
    ) +
    ggplot2::scale_y_continuous(limits = ylim_range, breaks = seq(ylim_range[1], ylim_range[2], 0.5)) +
    ggplot2::labs(title = title, x = "Visit", y = ylab, color = "") +
    common_theme
}

plots <- list(
  wfa = list(title = "Weight-for-age", ylab = "Weight-for-age z-score\n(95% CI)", ylim = c(-1.0, 1.0)),
  lfa = list(title = "Length-for-age", ylab = "Length-for-age z-score\n(95% CI)", ylim = c(-1.5, 1.0)),
  wfl = list(title = "Weight-for-length", ylab = "Weight-for-length z-score\n(95% CI)", ylim = c(-1.5, 1.0)),
  hcfa = list(title = "Head circumference-for-age", ylab = "Head circumference-for-age\nz-score (95% CI)", ylim = c(-1.0, 1.0))
)

p1 <- create_plot(df |> dplyr::filter(metric == "Weight-for-age"), plots$wfa$title, plots$wfa$ylab, plots$wfa$ylim)
p2 <- create_plot(df |> dplyr::filter(metric == "Length-for-age"), plots$lfa$title, plots$lfa$ylab, plots$lfa$ylim)
p3 <- create_plot(df |> dplyr::filter(metric == "Weight-for-length"), plots$wfl$title, plots$wfl$ylab, plots$wfl$ylim)
p4 <- create_plot(df |> dplyr::filter(metric == "Head circumference-for-age"), plots$hcfa$title, plots$hcfa$ylab, plots$hcfa$ylim)

combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig2A_growth_zscores.pdf"), plot = combined, width = 10, height = 10)

