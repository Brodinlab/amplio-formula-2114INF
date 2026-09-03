# Neutrophil-to-lymphocyte ratio (NLR), computed from the manually-gated CyTOF
# populations (cytof_manual_clean: outlier excluded, plate regressed out).
#
# NLR = Neutrophils / (B.cells + T.cells + NK)
# All three lymphocyte lineages (B cells, T cells, NK) and Neutrophils are
# each expressed as % of CD45+ in this table, so the CD45+ denominator cancels
# in the ratio -- this is equivalent to the ratio of absolute cell counts.
#
# Output:
#   output/tables/manualgating_nlr_values.csv
#   output/figures/manuscript/Fig5_manualgating_nlr_by_time.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)
timepoint_days <- c(V1 = 0, V3 = 60, V5 = 120)

meta <- base$metadata |>
  dplyr::mutate(
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    age_days = timepoint_days[as.character(timepoint)]
  )

nlr_cols <- c("Neutrophils", "B.cells", "T.cells", "NK")
stopifnot(all(nlr_cols %in% colnames(base$cytof_manual_clean)))

nlr_df <- base$cytof_manual_clean |>
  dplyr::transmute(
    cytof_id,
    lymphocytes = B.cells + T.cells + NK,
    neutrophils = Neutrophils,
    nlr = neutrophils / lymphocytes
  ) |>
  dplyr::inner_join(meta |> dplyr::select(cytof_id, subject_id, group_feeding, timepoint, age_days), by = "cytof_id") |>
  tidyr::drop_na(group_feeding, timepoint, nlr)

readr::write_csv(nlr_df, file.path(root, "output", "tables", "manualgating_nlr_values.csv"))

medians <- nlr_df |>
  dplyr::group_by(group_feeding, age_days, timepoint) |>
  dplyr::summarise(nlr = median(nlr, na.rm = TRUE), .groups = "drop")

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

p <- ggplot2::ggplot(nlr_df, ggplot2::aes(x = age_days, y = nlr)) +
  ggplot2::geom_line(ggplot2::aes(group = subject_id), color = "grey75", linewidth = 0.3, alpha = 0.5) +
  ggplot2::geom_line(
    data = medians,
    ggplot2::aes(x = age_days, y = nlr, color = group_feeding, group = group_feeding),
    linewidth = 1.3
  ) +
  ggplot2::geom_point(data = medians, ggplot2::aes(x = age_days, y = nlr, color = group_feeding), size = 2.2) +
  ggplot2::scale_color_manual(values = color_group, name = "Feeding Group") +
  ggplot2::scale_x_continuous(breaks = timepoint_days, labels = names(timepoint_days)) +
  ggplot2::labs(
    title = "Neutrophil-to-lymphocyte ratio over time: individual trajectories + group medians",
    x = "Timepoint", y = "NLR (Neutrophils / [B + T + NK])"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_nlr_by_time.pdf"), width = 6, height = 5)

cat("NLR summary by group x timepoint (median [IQR]):\n")
print(
  nlr_df |>
    dplyr::group_by(group_feeding, timepoint) |>
    dplyr::summarise(
      n = dplyr::n(), median_nlr = median(nlr), q25 = quantile(nlr, 0.25), q75 = quantile(nlr, 0.75),
      .groups = "drop"
    )
)
