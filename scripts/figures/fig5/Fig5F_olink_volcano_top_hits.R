suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "patchwork", "ggrepel", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

olink <- base$olink_npx_wide
meta <- base$olink_metadata
stopifnot("SampleID" %in% colnames(olink))

x <- as.data.frame(olink)
rownames(x) <- x$SampleID
x$SampleID <- NULL

timepoints <- c("V1", "V3", "V5")
protein_names <- colnames(x)

rows <- list()
for (tp in timepoints) {
  tp_ids <- meta$SampleID[meta$timepoint == tp]
  tp_x <- x[tp_ids, , drop = FALSE]
  tp_meta <- meta[meta$timepoint == tp, , drop = FALSE]
  for (protein in protein_names) {
    vals <- tp_x[[protein]]
    cf <- vals[tp_meta$group_feeding == "CF"]; cf <- cf[!is.na(cf)]
    ef <- vals[tp_meta$group_feeding == "EF"]; ef <- ef[!is.na(ef)]
    if (length(cf) >= 3 && length(ef) >= 3) {
      rows[[length(rows) + 1]] <- data.frame(
        Timepoint = tp,
        Protein = protein,
        log2FC = mean(ef, na.rm = TRUE) - mean(cf, na.rm = TRUE),
        p_value = stats::wilcox.test(cf, ef)$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}
df <- do.call(rbind, rows)
df$neg_log10_p <- -log10(df$p_value)
df$Category <- ifelse(df$p_value < 0.05, "p < 0.05", "Not significant")

timepoints <- c("V1", "V3", "V5")
df$Timepoint <- factor(df$Timepoint, levels = timepoints)

x_max <- max(abs(df$log2FC), na.rm = TRUE)
x_lim <- c(-x_max, x_max)
y_max <- max(df$neg_log10_p, na.rm = TRUE)
y_lim <- c(0, max(y_max * 1.05, 1))

plots <- lapply(seq_along(timepoints), function(i) {
  tp <- timepoints[i]
  panel_label <- LETTERS[i]
  d <- df |> dplyr::filter(Timepoint == tp)
  to_label <- d |> dplyr::filter(Category == "p < 0.05")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = log2FC, y = neg_log10_p)) +
    ggplot2::geom_point(ggplot2::aes(color = Category), size = 1.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c("Not significant" = "gray70", "p < 0.05" = "orange"),
      name = "Category"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    ggplot2::coord_cartesian(xlim = x_lim, ylim = y_lim) +
    ggplot2::labs(
      x = "log2 Fold Change (EF/CF)",
      y = "-log10(p-value)",
      title = paste0("Panel ", panel_label, ": ", tp)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (nrow(to_label) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = to_label,
      ggplot2::aes(label = Protein),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.5
    )
  }
  p
})

combined <- plots[[1]] + plots[[2]] + plots[[3]] +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(
    title = "Figure 4.13: Volcano Plots - Proteins by Feeding Group"
  )

out_dir <- file.path(root, "output/figures/manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(out_dir, "Fig5F_olink_volcano_top_hits.pdf"), plot = combined, width = 15, height = 5, units = "in", dpi = 300, device = "pdf")

