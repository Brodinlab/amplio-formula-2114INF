# Two ranked summaries across all 32 manually-gated populations:
#   1. Which cell types change most with age (V1/V3/V5), regardless of group.
#      eta^2 from a one-way ANOVA vs timepoint, pooled across both feeding
#      groups.
#   2. Which cell types differ most by feeding group, BY EFFECT SIZE regardless
#      of significance (none of these survive FDR at this sample size -- this
#      is for hypothesis generation / prioritizing what a larger study should
#      look at, not a claim of a real effect). Uses the change-from-baseline
#      Cohen's d (V1->V3, V1->V5) as the primary metric; the raw cross-sectional
#      d is included alongside for reference only.
#
# Uses cytof_manual_clean (scripts/export/export_cytof_manual_clean.R):
# outlier sample excluded and cytof_plate regressed out per population, per
# Petter's instruction 2026-09-03 -- both rankings now use fully batch-
# corrected, outlier-free data (previously only the MDS scripts were
# batch-corrected).
#
# Outputs:
#   output/tables/manualgating_age_effect_ranking.csv
#   output/tables/manualgating_feeding_group_effect_ranking.csv
#   output/figures/manuscript/Fig5_manualgating_age_effect_ranking.pdf
#   output/figures/manuscript/Fig5_manualgating_feeding_group_effect_ranking.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)
pop_cols <- setdiff(colnames(base$cytof_manual_clean), "cytof_id")

meta <- base$metadata |>
  dplyr::mutate(
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  )

df <- base$cytof_manual_clean |>
  dplyr::inner_join(meta |> dplyr::select(cytof_id, subject_id, group_feeding, timepoint), by = "cytof_id")

# ---- 1. Age effect ranking (eta^2 vs timepoint, pooled across groups) ----
age_ranking <- purrr::map_dfr(pop_cols, function(pop) {
  d <- df |> dplyr::select(value = dplyr::all_of(pop), timepoint) |> tidyr::drop_na()
  fit <- aov(value ~ timepoint, data = d)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  tibble::tibble(population = pop, eta_squared = ss[1] / sum(ss))
}) |>
  dplyr::arrange(dplyr::desc(eta_squared)) |>
  dplyr::mutate(rank = dplyr::row_number())

readr::write_csv(age_ranking, file.path(root, "output", "tables", "manualgating_age_effect_ranking.csv"))

p_age <- ggplot2::ggplot(age_ranking, ggplot2::aes(x = eta_squared, y = reorder(population, eta_squared))) +
  ggplot2::geom_col(fill = "#1D457FFF") +
  ggplot2::labs(
    title = "Which cell types change most with age (baseline/2 months/4 months)?",
    x = expression(eta^2 ~ "(one-way ANOVA vs. timepoint, both groups pooled)"), y = NULL
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"))

save_pdf(p_age, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_age_effect_ranking.pdf"), width = 6, height = 7)

# ---- 2. Feeding-group effect ranking (change-from-baseline Cohen's d, max |d| across V3/V5) ----
wide_by_subject <- df |>
  dplyr::select(subject_id, group_feeding, timepoint, dplyr::all_of(pop_cols)) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = dplyr::all_of(pop_cols), names_sep = "__")

delta_d <- purrr::map_dfr(pop_cols, function(pop) {
  v1_col <- paste0(pop, "__V1")
  purrr::map_dfr(c("V3", "V5"), function(tp) {
    tp_col <- paste0(pop, "__", tp)
    d <- wide_by_subject |> dplyr::transmute(group_feeding, delta = .data[[tp_col]] - .data[[v1_col]]) |> tidyr::drop_na()
    x <- d$delta[d$group_feeding == "SynF"]; y <- d$delta[d$group_feeding == "CtrlF"]
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(population = pop, follow_up_timepoint = tp, cohens_d = cd$d)
  })
})

feeding_group_ranking <- delta_d |>
  dplyr::group_by(population) |>
  dplyr::slice_max(abs(cohens_d), n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::rename(max_abs_d_timepoint = follow_up_timepoint, cohens_d_change_from_baseline = cohens_d) |>
  dplyr::arrange(dplyr::desc(abs(cohens_d_change_from_baseline))) |>
  dplyr::mutate(rank = dplyr::row_number())

# reference only: raw cross-sectional d at the same follow-up timepoint, not used for ranking
cross_sectional_for_ref <- purrr::map_dfr(pop_cols, function(pop) {
  purrr::map_dfr(c("V1", "V3", "V5"), function(tp) {
    d <- df |> dplyr::filter(timepoint == tp) |> dplyr::select(group_feeding, value = dplyr::all_of(pop)) |> tidyr::drop_na()
    x <- d$value[d$group_feeding == "SynF"]; y <- d$value[d$group_feeding == "CtrlF"]
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(population = pop, timepoint = tp, cohens_d_cross_sectional = cd$d)
  })
})

feeding_group_ranking <- feeding_group_ranking |>
  dplyr::left_join(
    cross_sectional_for_ref |> dplyr::rename(max_abs_d_timepoint = timepoint),
    by = c("population", "max_abs_d_timepoint")
  )

readr::write_csv(feeding_group_ranking, file.path(root, "output", "tables", "manualgating_feeding_group_effect_ranking.csv"))

plot_df <- feeding_group_ranking |>
  dplyr::mutate(population = factor(population, levels = rev(population)))

p_group <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cohens_d_change_from_baseline, y = population, fill = max_abs_d_timepoint)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = 0, color = "grey40") +
  ggplot2::scale_fill_manual(
    values = c(V3 = "#C36377FF", V5 = "#1D457FFF"),
    labels = TIMEPOINT_LABELS[c("V3", "V5")],
    name = "Follow-up"
  ) +
  ggplot2::labs(
    title = "Which cell types differ most by feeding group?\n(by effect size, change from baseline -- not gated on significance)",
    x = "Cohen's d (SynF vs CtrlF, change from baseline)", y = NULL
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p_group, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_feeding_group_effect_ranking.pdf"), width = 6.5, height = 7)

cat("Top 5 age-associated populations:\n")
print(age_ranking |> dplyr::slice_head(n = 5))
cat("\nTop 5 feeding-group-associated populations (by |Cohen's d|, change from baseline):\n")
print(feeding_group_ranking |> dplyr::slice_head(n = 5) |> dplyr::select(population, max_abs_d_timepoint, cohens_d_change_from_baseline, cohens_d_cross_sectional))
