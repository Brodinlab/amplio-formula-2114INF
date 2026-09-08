# Olink NPX vs. SynF/CtrlF feeding group. Replaces the broken
# Fig5F_olink_volcano_top_hits.R (wrong source path; raw NPX mean-difference
# instead of Cohen's d; no multiple-testing correction; no omnibus/AUC/
# change-from-baseline tests) -- upgraded to the same statistical framework
# used for the manually-gated CyTOF populations (Fig5_manualgating_synf_ctrlf_stats.R),
# applied here to all 358 Olink proteins (olink_npx_clean).
#
# Statistical design (identical to the CyTOF script, see its header for full
# rationale):
#   - Cross-sectional: Wilcoxon rank-sum test per protein, per time point (V1/V3/V5)
#   - Longitudinal: linear mixed-effects model, NPX ~ group * timepoint + delivery + (1|subject),
#     with an omnibus likelihood-ratio test on the group:timepoint interaction
#   - Subject-level summary: trapezoidal AUC across each subject's available time points,
#     compared SynF vs CtrlF by Wilcoxon
#   - Change from baseline (V1): each subject's own delta at V3/V5, compared SynF vs CtrlF
#   - Benjamini-Hochberg FDR within each test family; FDR<0.05 significant, nominal p<0.05
#     exploratory, per manuscript convention
#
# Primary figure: volcano plot per time point (Cohen's d vs -log10 nominal p),
# the natural visualization for ~358 features (a per-protein forest plot, as
# used for the 27-population CyTOF figure, would be unreadable at this scale).
#
# Outputs:
#   output/tables/olink_synf_ctrlf_cross_sectional.csv
#   output/tables/olink_synf_ctrlf_omnibus_interaction.csv
#   output/tables/olink_synf_ctrlf_auc.csv
#   output/tables/olink_synf_ctrlf_change_from_baseline.csv
#   output/figures/manuscript/Fig5_olink_volcano.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "lme4", "lmerTest", "ggplot2", "ggrepel"))

root <- get_repo_root()
base <- load_base_tables(root)

protein_cols <- setdiff(colnames(base$olink_npx_clean), "SampleID")

meta <- base$olink_metadata |>
  dplyr::mutate(
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  )

df <- base$olink_npx_clean |>
  dplyr::inner_join(
    meta |> dplyr::select(SampleID, subject_id, group_feeding, group_delivery, timepoint),
    by = "SampleID"
  )

stopifnot(nrow(df) == nrow(base$olink_npx_clean)) # every Olink sample should have metadata

# ---- 1. Cross-sectional: Wilcoxon + Cohen's d per protein per timepoint ----
cross_sectional <- purrr::map_dfr(protein_cols, function(protein) {
  purrr::map_dfr(levels(df$timepoint), function(tp) {
    d <- df |> dplyr::filter(timepoint == tp) |> dplyr::select(group_feeding, value = dplyr::all_of(protein))
    x <- d$value[d$group_feeding == "SynF"]
    y <- d$value[d$group_feeding == "CtrlF"]
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(
      protein = protein, timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
      p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
    )
  })
}) |>
  dplyr::group_by(timepoint) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

# ---- 2. Omnibus group:timepoint interaction (LRT, categorical timepoint) ----
omnibus <- purrr::map_dfr(protein_cols, function(protein) {
  d <- df |> dplyr::select(value = dplyr::all_of(protein), group_feeding, timepoint, group_delivery, subject_id) |>
    tidyr::drop_na()
  full_model <- lme4::lmer(value ~ group_feeding * timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE)
  reduced_model <- lme4::lmer(value ~ group_feeding + timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE)
  a <- anova(reduced_model, full_model)
  tibble::tibble(protein = protein, lrt_chisq = a$Chisq[2], lrt_df = a$Df[2], p_value = a$`Pr(>Chisq)`[2])
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

# ---- 3. Subject-level AUC (trapezoidal, across each subject's available timepoints) ----
auc_by_subject <- purrr::map_dfr(protein_cols, function(protein) {
  df |>
    dplyr::mutate(t_days = TIMEPOINT_DAYS[as.character(timepoint)]) |>
    dplyr::select(subject_id, group_feeding, t_days, value = dplyr::all_of(protein)) |>
    tidyr::drop_na() |>
    dplyr::group_by(subject_id, group_feeding) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::summarise(
      auc = { o <- order(t_days); sum(diff(t_days[o]) * (head(value[o], -1) + tail(value[o], -1)) / 2) / (max(t_days) - min(t_days)) },
      .groups = "drop"
    ) |>
    dplyr::mutate(protein = protein)
})

auc_test <- auc_by_subject |>
  dplyr::group_by(protein) |>
  dplyr::summarise(
    n_synf = sum(group_feeding == "SynF"), n_ctrlf = sum(group_feeding == "CtrlF"),
    p_value = tryCatch(suppressWarnings(wilcox.test(auc ~ group_feeding)$p.value), error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

# ---- 4. Change from baseline (V1) ----
wide_by_subject <- df |>
  dplyr::select(subject_id, group_feeding, timepoint, dplyr::all_of(protein_cols)) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = dplyr::all_of(protein_cols), names_sep = "__")

delta_from_baseline <- purrr::map_dfr(protein_cols, function(protein) {
  v1_col <- paste0(protein, "__V1")
  purrr::map_dfr(c("V3", "V5"), function(tp) {
    tp_col <- paste0(protein, "__", tp)
    d <- wide_by_subject |> dplyr::transmute(group_feeding, delta = .data[[tp_col]] - .data[[v1_col]]) |> tidyr::drop_na()
    x <- d$delta[d$group_feeding == "SynF"]
    y <- d$delta[d$group_feeding == "CtrlF"]
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(
      protein = protein, follow_up_timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
      p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
    )
  })
}) |>
  dplyr::group_by(follow_up_timepoint) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

# ---- Save results tables ----
out_table_dir <- file.path(root, "output", "tables")
readr::write_csv(cross_sectional, file.path(out_table_dir, "olink_synf_ctrlf_cross_sectional.csv"))
readr::write_csv(omnibus, file.path(out_table_dir, "olink_synf_ctrlf_omnibus_interaction.csv"))
readr::write_csv(auc_test, file.path(out_table_dir, "olink_synf_ctrlf_auc.csv"))
readr::write_csv(delta_from_baseline, file.path(out_table_dir, "olink_synf_ctrlf_change_from_baseline.csv"))

# ---- Primary figure: volcano plot per timepoint (Cohen's d vs nominal p) ----
volcano_df <- cross_sectional |>
  dplyr::mutate(
    neg_log10_p = -log10(p_value),
    timepoint = relabel_timepoint(timepoint),
    significant = p_fdr < 0.05
  )

# Label the top 10 most differentially regulated proteins per panel (lowest
# nominal p-value), rather than a fixed p-value cutoff -- a fixed threshold
# left some panels (e.g. 2 months, where nothing crosses p<0.01) with no
# labels at all, per Petter's request, 2026-09-06.
N_LABEL <- 10
to_label <- volcano_df |>
  dplyr::group_by(timepoint) |>
  dplyr::slice_min(p_value, n = N_LABEL) |>
  dplyr::ungroup()

p_volcano <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = cohens_d, y = neg_log10_p)) +
  ggplot2::geom_point(ggplot2::aes(color = significant), size = 1.3, alpha = 0.6) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = to_label, ggplot2::aes(label = protein), size = 2.2, max.overlaps = 20, segment.size = 0.2) +
  ggplot2::facet_wrap(~timepoint, nrow = 1) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#C36377FF", `FALSE` = "grey40"), name = "FDR<0.05") +
  ggplot2::labs(x = "Cohen's d (SynF vs CtrlF)", y = expression(-log[10]~"(nominal p-value)")) +
  ggplot2::theme_bw(base_size = 9)

save_pdf(p_volcano, file.path(root, "output", "figures", "manuscript", "Fig5_olink_volcano.pdf"), width = 11, height = 4.5)

# ---- Baseline-normalized figure: volcano on CHANGE FROM BASELINE (Petter's request, 2026-09-06) ----
# The cross-sectional volcano above compares raw NPX levels at each timepoint,
# which conflates a true treatment effect with any pre-existing baseline (V1)
# imbalance between groups -- several proteins are already nominally different
# at Baseline in the plot above (e.g. PLAUR), before the synbiotic could have
# acted. Using each subject's own delta from their V1 value as the comparator
# (delta_from_baseline, computed above) isolates the treatment-period response
# and is immune to a baseline offset by construction -- same rationale as the
# CyTOF change-from-baseline analysis (Fig5_manualgating_synf_ctrlf_stats.R).
delta_volcano_df <- delta_from_baseline |>
  dplyr::mutate(
    neg_log10_p = -log10(p_value),
    follow_up_timepoint = factor(
      follow_up_timepoint, levels = c("V3", "V5"),
      labels = paste0("Baseline -> ", TIMEPOINT_LABELS[c("V3", "V5")])
    ),
    # Direction, not FDR significance, is the informative encoding here (0
    # proteins are FDR-significant, so that color mapping was flat/uninformative
    # -- per the recorded figure-style convention, a legend should only encode
    # something that actually varies). Colored to match the same CtrlF/SynF
    # palette used everywhere else in this repo (color_group), so "which group
    # is higher" reads consistently across figures rather than introducing a
    # new ad hoc color scheme -- per Petter's request, 2026-09-06.
    higher_in = ifelse(cohens_d < 0, "CtrlF", "SynF")
  )

N_LABEL_DELTA <- 15 # more labels than the raw cross-sectional volcano's N_LABEL, per Petter's request, 2026-09-06
to_label_delta <- delta_volcano_df |>
  dplyr::group_by(follow_up_timepoint) |>
  dplyr::slice_min(p_value, n = N_LABEL_DELTA) |>
  dplyr::ungroup()

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

# Tufte-style cleanup (per the figure-style conventions established today,
# observations/figure-style-conventions.md): panel border and gridlines are
# non-data ink and are minimized/removed rather than left at ggplot defaults;
# facet strip labels are kept since they are essential wayfinding (which
# panel is which timepoint), not decorative chrome.
p_volcano_delta <- ggplot2::ggplot(delta_volcano_df, ggplot2::aes(x = cohens_d, y = neg_log10_p)) +
  ggplot2::geom_point(ggplot2::aes(color = higher_in), size = 1.3, alpha = 0.7) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = to_label_delta, ggplot2::aes(label = protein), size = 2.2, color = "grey20", max.overlaps = 20, segment.size = 0.2) +
  ggplot2::facet_wrap(~follow_up_timepoint, nrow = 1) +
  ggplot2::scale_color_manual(values = color_group, name = "Higher in") +
  ggplot2::labs(x = "Cohen's d, change from baseline (SynF vs CtrlF)", y = expression(-log[10]~"(nominal p-value)")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    # Further Tufte cleanup, per Petter's request, 2026-09-06: the panel/plot
    # background fill is itself non-data ink -- removed (transparent) rather
    # than left white. Gridlines are dropped entirely rather than just
    # muted: the two dashed reference lines (x=0, p=0.05) already give the
    # reader the reference points gridlines would otherwise provide, so a
    # full grid is redundant ink on top of that. The facet strip keeps its
    # text (essential wayfinding -- which panel is which timepoint) but
    # loses its background fill.
    panel.border = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, color = NA),
    plot.background = ggplot2::element_rect(fill = NA, color = NA),
    axis.line = ggplot2::element_line(color = "grey40", linewidth = 0.3),
    axis.ticks = ggplot2::element_line(color = "grey40", linewidth = 0.3),
    strip.background = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(fill = NA, color = NA),
    legend.key = ggplot2::element_rect(fill = NA, color = NA)
  )

save_pdf(p_volcano_delta, file.path(root, "output", "figures", "manuscript", "Fig5_olink_volcano_change_from_baseline.pdf"), width = 8, height = 4.5)

cat("Cross-sectional: ", sum(cross_sectional$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(cross_sectional),
    " (protein x timepoint) comparisons FDR-significant; ",
    sum(cross_sectional$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Omnibus interaction: ", sum(omnibus$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(omnibus),
    " proteins FDR-significant; ", sum(omnibus$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Subject-level AUC: ", sum(auc_test$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(auc_test),
    " proteins FDR-significant; ", sum(auc_test$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
cat("Change from baseline: ", sum(delta_from_baseline$p_fdr < 0.05, na.rm = TRUE), " of ", nrow(delta_from_baseline),
    " (protein x follow-up timepoint) comparisons FDR-significant; ",
    sum(delta_from_baseline$p_value < 0.05, na.rm = TRUE), " nominal p<0.05\n", sep = "")
