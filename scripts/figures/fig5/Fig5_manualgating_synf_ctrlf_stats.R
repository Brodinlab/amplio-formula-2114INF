# Manually-gated CyTOF (analysis_20260831_v1.1_relative_frequencies.csv) vs. SynF/CtrlF feeding group.
#
# Input gating: Kanth's QC'd v1.1 pass (2026-09-03), superseding the preliminary
# analysis_20260831.csv run (2026-09-02). Results were consistent in direction
# and magnitude between the two passes -- see Fig6_manualgating_pseudotime.R
# header for the same comparison on the pseudotime side.
#
# Uses cytof_manual_clean (scripts/export/export_cytof_manual_clean.R):
# outlier sample excluded (subject PHL001-0073, V1, CtrlF) and cytof_plate
# regressed out per population (percentage-scale residuals + global mean) --
# per Petter's instruction, 2026-09-03. group_delivery is kept as an
# additional covariate in the LME (a separate, already-checked-non-confounded
# factor); plate is NOT also added as a model term since it is now already
# regressed out of the input values -- doing both would double-correct.
#
# Statistical design mirrors the manuscript's own "Blood immune marker analyses" Methods:
#   - Cross-sectional: Wilcoxon rank-sum test per population, per timepoint (V1/V3/V5)
#   - Longitudinal: linear mixed-effects model, freq ~ group * timepoint + delivery + (1|subject),
#     with an omnibus likelihood-ratio test on the group:timepoint interaction (categorical
#     timepoint, not a single linear slope) -- tests whether the group-difference pattern
#     differs anywhere across the full V1/V3/V5 series, not just a linear trend.
#   - Complementary subject-level summary: trapezoidal AUC across each subject's available
#     timepoints, per population, compared SynF vs CtrlF by Wilcoxon -- easier to visualize
#     (spaghetti + AUC boxplot) than the LME itself.
#   - Benjamini-Hochberg FDR within each test family (cross-sectional per timepoint; omnibus
#     interaction; AUC) -- FDR<0.05 significant, nominal p<0.05 exploratory, per manuscript
#     convention.
#
# Outputs:
#   output/tables/manualgating_synf_ctrlf_stats.csv   (all per-population results)
#   output/figures/manuscript/Fig5_manualgating_effect_sizes.pdf

suppressPackageStartupMessages({
  # NOTE: other scripts in this repo source("code/lib/...") but the actual
  # directory is scripts/lib/ (pre-existing path mismatch upstream, not fixed
  # here since that's out of scope for this script) -- using the real path.
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "lme4", "lmerTest", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

pop_cols <- setdiff(colnames(base$cytof_manual_clean), "cytof_id")

meta <- base$metadata |>
  dplyr::mutate(
    group_feeding = factor(
      dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"),
      levels = c("CtrlF", "SynF")
    ),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  )

df <- base$cytof_manual_clean |>
  dplyr::inner_join(
    meta |> dplyr::select(cytof_id, subject_id, group_feeding, group_delivery, timepoint),
    by = "cytof_id"
  )

stopifnot(nrow(df) == nrow(base$cytof_manual_clean)) # every manually-gated sample should have metadata

# ---- 1. Cross-sectional: Wilcoxon + Cohen's d per population per timepoint ----
cross_sectional <- purrr::map_dfr(pop_cols, function(pop) {
  purrr::map_dfr(levels(df$timepoint), function(tp) {
    d <- df |> dplyr::filter(timepoint == tp) |> dplyr::select(group_feeding, value = dplyr::all_of(pop))
    x <- d$value[d$group_feeding == "SynF"]
    y <- d$value[d$group_feeding == "CtrlF"]
    if (length(x) < 3 || length(y) < 3) {
      return(tibble::tibble(
        population = pop, timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
        p_value = NA_real_, cohens_d = NA_real_, d_ci_lower = NA_real_, d_ci_upper = NA_real_
      ))
    }
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(
      population = pop, timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
      p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
    )
  })
}) |>
  dplyr::group_by(timepoint) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

# ---- 2. Omnibus group:timepoint interaction (LRT, categorical timepoint) ----
omnibus <- purrr::map_dfr(pop_cols, function(pop) {
  d <- df |> dplyr::select(value = dplyr::all_of(pop), group_feeding, timepoint, group_delivery, subject_id) |>
    tidyr::drop_na()
  if (dplyr::n_distinct(d$subject_id) < 10 || dplyr::n_distinct(d$timepoint) < 2) {
    return(tibble::tibble(population = pop, lrt_chisq = NA_real_, lrt_df = NA_integer_, p_value = NA_real_))
  }
  full_model <- tryCatch(
    lme4::lmer(value ~ group_feeding * timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE),
    error = function(e) NULL
  )
  reduced_model <- tryCatch(
    lme4::lmer(value ~ group_feeding + timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(full_model) || is.null(reduced_model)) {
    return(tibble::tibble(population = pop, lrt_chisq = NA_real_, lrt_df = NA_integer_, p_value = NA_real_))
  }
  a <- anova(reduced_model, full_model)
  tibble::tibble(
    population = pop,
    lrt_chisq = a$Chisq[2],
    lrt_df = a$Df[2],
    p_value = a$`Pr(>Chisq)`[2]
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

# ---- 3. Subject-level AUC (trapezoidal, across each subject's available timepoints) ----
timepoint_days <- c(V1 = 0, V3 = 60, V5 = 120) # nominal days-from-enrollment spacing

trapezoidal_auc <- function(t, v) {
  if (length(t) < 2) return(NA_real_)
  o <- order(t)
  t <- t[o]; v <- v[o]
  sum(diff(t) * (head(v, -1) + tail(v, -1)) / 2) / (max(t) - min(t)) # time-normalized AUC
}

auc_by_subject <- purrr::map_dfr(pop_cols, function(pop) {
  df |>
    dplyr::mutate(t_days = timepoint_days[as.character(timepoint)]) |>
    dplyr::select(subject_id, group_feeding, t_days, value = dplyr::all_of(pop)) |>
    tidyr::drop_na() |>
    dplyr::group_by(subject_id, group_feeding) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::summarise(auc = trapezoidal_auc(t_days, value), .groups = "drop") |>
    dplyr::mutate(population = pop)
})

auc_test <- auc_by_subject |>
  dplyr::group_by(population) |>
  dplyr::summarise(
    n_synf = sum(group_feeding == "SynF"),
    n_ctrlf = sum(group_feeding == "CtrlF"),
    p_value = tryCatch(
      suppressWarnings(wilcox.test(auc ~ group_feeding)$p.value),
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

# ---- 4. Change from baseline (V1): does SynF vs CtrlF diverge relative to each
# subject's OWN starting point, rather than in raw cross-sectional level? This
# is the right test for whether a difference reflects the intervention, since a
# raw V3/V5 group difference can't distinguish "the formula caused this" from
# "the groups already differed at V1" -- V1 is the baseline/enrollment visit,
# before the synbiotic can have acted (see Fig5 discussion, 2026-09-03).
wide_by_subject <- df |>
  dplyr::select(subject_id, group_feeding, group_delivery, timepoint, dplyr::all_of(pop_cols)) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = dplyr::all_of(pop_cols), names_sep = "__")

delta_from_baseline <- purrr::map_dfr(pop_cols, function(pop) {
  v1_col <- paste0(pop, "__V1")
  purrr::map_dfr(c("V3", "V5"), function(tp) {
    tp_col <- paste0(pop, "__", tp)
    if (!all(c(v1_col, tp_col) %in% colnames(wide_by_subject))) {
      return(tibble::tibble(population = pop, follow_up_timepoint = tp, n_synf = 0L, n_ctrlf = 0L,
                             p_value = NA_real_, cohens_d = NA_real_, d_ci_lower = NA_real_, d_ci_upper = NA_real_))
    }
    d <- wide_by_subject |>
      dplyr::transmute(group_feeding, delta = .data[[tp_col]] - .data[[v1_col]]) |>
      tidyr::drop_na()
    x <- d$delta[d$group_feeding == "SynF"]
    y <- d$delta[d$group_feeding == "CtrlF"]
    if (length(x) < 3 || length(y) < 3) {
      return(tibble::tibble(population = pop, follow_up_timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
                             p_value = NA_real_, cohens_d = NA_real_, d_ci_lower = NA_real_, d_ci_upper = NA_real_))
    }
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(
      population = pop, follow_up_timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
      p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
    )
  })
}) |>
  dplyr::group_by(follow_up_timepoint) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

# ---- Save results table ----
out_table_dir <- file.path(root, "output", "tables")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

results <- list(
  cross_sectional = cross_sectional,
  omnibus_interaction = omnibus,
  subject_level_auc = auc_test
)
readr::write_csv(cross_sectional, file.path(out_table_dir, "manualgating_synf_ctrlf_cross_sectional.csv"))
readr::write_csv(omnibus, file.path(out_table_dir, "manualgating_synf_ctrlf_omnibus_interaction.csv"))
readr::write_csv(auc_test, file.path(out_table_dir, "manualgating_synf_ctrlf_auc.csv"))
readr::write_csv(delta_from_baseline, file.path(out_table_dir, "manualgating_synf_ctrlf_change_from_baseline.csv"))

# ---- Summary figure: effect size (Cohen's d) per population, faceted by timepoint ----
plot_df <- cross_sectional |>
  dplyr::filter(!is.na(cohens_d)) |>
  dplyr::mutate(population = factor(population, levels = rev(pop_cols)), timepoint = relabel_timepoint(timepoint))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cohens_d, y = population, color = p_fdr < 0.05)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = d_ci_lower, xmax = d_ci_upper), height = 0.2) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~timepoint, nrow = 1) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#C36377FF", `FALSE` = "grey40"), name = "FDR<0.05") +
  ggplot2::labs(
    x = "Cohen's d (SynF vs CtrlF)", y = NULL,
    title = "Manually-gated CyTOF: SynF vs CtrlF effect size by timepoint"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_effect_sizes.pdf"), width = 10, height = 7)

# ---- Summary figure: change-from-baseline effect size, faceted by follow-up timepoint ----
delta_plot_df <- delta_from_baseline |>
  dplyr::filter(!is.na(cohens_d)) |>
  dplyr::mutate(
    population = factor(population, levels = rev(pop_cols)),
    follow_up_timepoint = factor(
      follow_up_timepoint,
      levels = c("V3", "V5"),
      labels = paste0("Baseline -> ", TIMEPOINT_LABELS[c("V3", "V5")])
    )
  )

p_delta <- ggplot2::ggplot(delta_plot_df, ggplot2::aes(x = cohens_d, y = population, color = p_fdr < 0.05)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = d_ci_lower, xmax = d_ci_upper), height = 0.2) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~follow_up_timepoint, nrow = 1) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#C36377FF", `FALSE` = "grey40"), name = "FDR<0.05") +
  ggplot2::labs(
    x = "Cohen's d, change from baseline (SynF vs CtrlF)", y = NULL,
    title = "Manually-gated CyTOF: SynF vs CtrlF, change from baseline"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"))

save_pdf(p_delta, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_change_from_baseline.pdf"), width = 7, height = 7)

cat("Cross-sectional: ", sum(cross_sectional$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(cross_sectional),
    " (population x timepoint) comparisons FDR-significant; ",
    sum(cross_sectional$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Omnibus interaction: ", sum(omnibus$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(omnibus),
    " populations FDR-significant; ", sum(omnibus$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Subject-level AUC: ", sum(auc_test$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(auc_test),
    " populations FDR-significant; ", sum(auc_test$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Change from baseline: ", sum(delta_from_baseline$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(delta_from_baseline),
    " (population x follow-up timepoint) comparisons FDR-significant; ",
    sum(delta_from_baseline$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
