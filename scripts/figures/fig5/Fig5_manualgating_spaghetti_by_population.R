# Spaghetti plot, one panel per manually-gated population: individual subject
# trajectories (light background lines, colored by feeding group per Petter's
# request 2026-09-03) with SynF/CtrlF group medians overlaid as bold colored
# lines in the same palette. Uses cytof_manual_clean (outlier excluded, plate
# regressed out per population) -- same input as the other Fig5/Fig6
# population-level analyses.
#
# Output: output/figures/manuscript/Fig5_manualgating_spaghetti_by_population.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

pop_cols <- setdiff(colnames(base$cytof_manual_clean), "cytof_id")

meta <- base$metadata |>
  dplyr::mutate(
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    age_days = TIMEPOINT_DAYS[as.character(timepoint)]
  )

long <- base$cytof_manual_clean |>
  tidyr::pivot_longer(cols = dplyr::all_of(pop_cols), names_to = "population", values_to = "value") |>
  dplyr::inner_join(meta |> dplyr::select(cytof_id, subject_id, group_feeding, timepoint, age_days), by = "cytof_id") |>
  tidyr::drop_na(group_feeding, timepoint, value) |>
  dplyr::mutate(population = factor(population, levels = pop_cols))

medians <- long |>
  dplyr::group_by(population, group_feeding, age_days, timepoint) |>
  dplyr::summarise(value = median(value, na.rm = TRUE), .groups = "drop")

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

p <- ggplot2::ggplot(long, ggplot2::aes(x = age_days, y = value, color = group_feeding)) +
  ggplot2::geom_line(ggplot2::aes(group = subject_id), linewidth = 0.25, alpha = 0.25) +
  ggplot2::geom_line(
    data = medians,
    ggplot2::aes(x = age_days, y = value, color = group_feeding, group = group_feeding),
    linewidth = 1.3
  ) +
  ggplot2::geom_point(
    data = medians,
    ggplot2::aes(x = age_days, y = value, color = group_feeding),
    size = 1.8
  ) +
  ggplot2::scale_color_manual(values = color_group, name = "Feeding Group") +
  ggplot2::scale_x_continuous(breaks = TIMEPOINT_DAYS, labels = TIMEPOINT_LABELS) +
  ggplot2::facet_wrap(~population, scales = "free_y", ncol = 5) +
  ggplot2::labs(
    title = "Manually-gated CyTOF populations over time: individual trajectories + group medians",
    x = "Timepoint", y = "% of CD45+ (plate-corrected)"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
    strip.text = ggplot2::element_text(size = 6.5)
  )

save_pdf(p, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_spaghetti_by_population.pdf"), width = 15, height = 18)

cat("n populations:", length(pop_cols), " | n samples:", dplyr::n_distinct(long$cytof_id), "\n")
