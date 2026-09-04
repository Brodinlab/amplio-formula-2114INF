# Immune maturation pseudotime from the manually-gated CyTOF table, mirroring the
# manuscript's existing FlowSOM-cluster-based pseudotime analysis (Fig 6d/f):
#
# Input gating: Kanth's QC'd v1.1 pass (2026-09-03), superseding the preliminary
# analysis_20260831.csv run (2026-09-02). Age-correlation was essentially
# unchanged after QC (0.473 -> 0.467), and the same 6 populations top the
# time-varying ranking in both passes (just reordered) -- so the weaker/noisier
# signal relative to the published FlowSOM-based pseudotime is NOT explained by
# the preliminary-gating hypothesis; it more likely reflects manual gating
# capturing less age-related variance than the 94 unsupervised FlowSOM clusters
# (fewer, coarser population definitions), independent of QC completeness.
#
# Uses the same outlier-exclusion + plate-correction as cytof_manual_clean
# (scripts/export/export_cytof_manual_clean.R) but reconstructed locally
# over ALL 32 raw populations plus NLR as a 33rd candidate (Petter's
# instruction, 2026-09-04) -- NOT the repo-wide 28-population filtered set
# (MANUAL_GATING_EXCLUDED_POPULATIONS, common.R), which is reversed here
# specifically for the pseudotime's time-association population selection.
#
#   "Top 20% of time varying immune cell clusters used to construct a pseudotime
#    metric ... Embedding using PCA ... coloring by pseudotime or actual age."
#   "We compared pseudotime distributions ... EF-associated trajectories were on
#    average lower in pseudotime at both 2 and 4 months."
#
# This rebuilds the same pseudotime construction from the 32 manually-gated
# populations instead of the 94 unsupervised FlowSOM clusters, as a complementary/
# validating result for whether SynF vs CtrlF differ across the full V1->V3->V5
# trajectory (not just at single timepoints -- see Fig5_manualgating_synf_ctrlf_stats.R
# for the population-by-population cross-sectional/omnibus-interaction results).
#
# Method:
#   1. Rank the 32 populations by how strongly they vary with (nominal) age --
#      one-way ANOVA eta-squared across V1/V3/V5, computed on ALL samples pooled
#      (group-blind), matching the manuscript's "time varying" (not "group
#      varying") selection criterion.
#   2. Take the top 20% (~6-7 populations), scale, PCA.
#   3. Fit a principal curve (princurve) through the top-2-PC embedding -- this is
#      the same core single-lineage algorithm Slingshot uses when there is no
#      branching, which is the case here (one lineage: age).
#   4. Orient pseudotime to increase with age (flip sign if anti-correlated).
#   5. Test SynF vs CtrlF: cross-sectional Wilcoxon per timepoint (mirrors the
#      manuscript's own comparison), plus the same omnibus group:timepoint LRT
#      used in the population-level analysis, for consistency.
#
# Outputs:
#   output/tables/manualgating_pseudotime_values.csv               (per-sample pseudotime)
#   output/tables/manualgating_pseudotime_group_tests.csv          (cross-sectional + omnibus)
#   output/tables/manualgating_pseudotime_change_from_baseline.csv (per-subject delta vs V1)
#   output/figures/manuscript/Fig6_manualgating_pseudotime_pca.pdf
#   output/figures/manuscript/Fig6_manualgating_pseudotime_by_group.pdf
#   output/figures/manuscript/Fig6_manualgating_pseudotime_change_from_baseline.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "lme4", "princurve", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# Reverse the global 28-population filter (common.R
# MANUAL_GATING_EXCLUDED_POPULATIONS) for pseudotime specifically -- Petter's
# instruction, 2026-09-04: the time-association ranking below should draw
# from all 32 raw populations, not the reduced set used elsewhere in this
# repo. Reconstructed locally from the raw table (never modified) via the
# same outlier-exclusion + plate-correction used to build the shared "clean"
# table (scripts/export/export_cytof_manual_clean.R, plate_correct_log_scale
# in common.R), so this is on identical footing otherwise. Also adds NLR
# (Neutrophils / [B.cells+T.cells+NK], same definition as
# Fig5_manualgating_nlr.R) as an additional candidate population, per
# Petter's instruction -- so it can compete for inclusion in the top-20%
# time-varying set like any gated population.
raw_all_pops <- readr::read_csv(file.path(root, "data", "tables", "cytof_manual_gating_frequency.csv"), show_col_types = FALSE) |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::filter(!cytof_id %in% MANUAL_GATING_OUTLIER_IDS) |>
  dplyr::left_join(base$metadata |> dplyr::select(cytof_id, cytof_plate), by = "cytof_id")
stopifnot(!any(is.na(raw_all_pops$cytof_plate)))

all_pop_cols <- setdiff(colnames(raw_all_pops), c("cytof_id", "cytof_plate"))
cytof_manual_clean_full <- plate_correct_log_scale(raw_all_pops, all_pop_cols, raw_all_pops$cytof_plate) |>
  dplyr::select(cytof_id, dplyr::all_of(all_pop_cols))

cytof_manual_clean_full$NLR <- cytof_manual_clean_full$Neutrophils /
  (cytof_manual_clean_full$B.cells + cytof_manual_clean_full$T.cells + cytof_manual_clean_full$NK)
cytof_manual_clean_full$NLR[!is.finite(cytof_manual_clean_full$NLR)] <- NA_real_

pop_cols <- setdiff(colnames(cytof_manual_clean_full), "cytof_id")
timepoint_days <- c(V1 = 0, V3 = 60, V5 = 120)

meta <- base$metadata |>
  dplyr::mutate(
    group_feeding = factor(
      dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"),
      levels = c("CtrlF", "SynF")
    ),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    age_days = timepoint_days[as.character(timepoint)]
  )

df <- cytof_manual_clean_full |>
  dplyr::inner_join(
    meta |> dplyr::select(cytof_id, subject_id, group_feeding, group_delivery, timepoint, age_days),
    by = "cytof_id"
  ) |>
  tidyr::drop_na(group_feeding, timepoint)

# ---- 1. Rank populations by age-association (eta-squared, one-way ANOVA vs timepoint) ----
eta_squared <- purrr::map_dbl(pop_cols, function(pop) {
  d <- df |> dplyr::select(value = dplyr::all_of(pop), timepoint) |> tidyr::drop_na()
  fit <- tryCatch(aov(value ~ timepoint, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  ss[1] / sum(ss) # eta^2 = SS_between / SS_total
})
names(eta_squared) <- pop_cols

n_top <- max(2, round(0.20 * length(pop_cols)))

# Redundancy guard: walking the eta^2-sorted list top-down, skip a candidate
# whose correlation with an already-selected population exceeds 0.9 --
# without this, near-perfect complement pairs from the gating hierarchy
# (e.g. Neutrophils/NonNeutrophils, r=-0.9995: they sum to ~100% by
# construction) can both get selected, feeding the PCA two duplicate
# (negated) dimensions. That degeneracy destabilized the principal-curve
# fit and produced a bimodal pseudotime distribution (Petter, 2026-09-04) --
# caught by checking this correlation directly once he flagged the bimodal
# shapes. Ranking still considers all populations in pop_cols (unchanged
# per Petter's instruction); only the greedy pick is now redundancy-aware.
REDUNDANCY_THRESHOLD <- 0.9
ranked_populations <- names(sort(eta_squared, decreasing = TRUE))
top_populations <- character(0)
for (pop in ranked_populations) {
  if (length(top_populations) >= n_top) break
  is_redundant <- length(top_populations) > 0 && any(purrr::map_lgl(top_populations, function(p) {
    abs(cor(df[[pop]], df[[p]], use = "complete.obs")) >= REDUNDANCY_THRESHOLD
  }))
  if (!is_redundant) top_populations <- c(top_populations, pop)
}

cat("Top", n_top, "time-varying populations (eta^2 vs timepoint, redundancy-filtered):\n")
print(round(eta_squared[top_populations], 3))

# ---- 2. PCA on the top time-varying populations, all samples pooled ----
mat <- df |>
  dplyr::select(cytof_id, dplyr::all_of(top_populations)) |>
  tidyr::drop_na()

pca_input <- scale(as.matrix(mat |> dplyr::select(-cytof_id)))
pca <- prcomp(pca_input, center = FALSE, scale. = FALSE) # already scaled above

# ---- 3. Principal curve through PC1-PC2 (single-lineage pseudotime, as in Slingshot) ----
curve <- princurve::principal_curve(pca$x[, 1:2])
pseudotime_raw <- curve$lambda

pt_df <- tibble::tibble(cytof_id = mat$cytof_id, pseudotime_raw = pseudotime_raw) |>
  dplyr::left_join(df |> dplyr::select(cytof_id, subject_id, group_feeding, group_delivery, timepoint, age_days),
                    by = "cytof_id")

# ---- 4. Orient pseudotime to increase with age ----
orientation <- sign(cor(pt_df$pseudotime_raw, pt_df$age_days, use = "complete.obs", method = "spearman"))
if (orientation == 0) orientation <- 1
pt_df$pseudotime <- if (orientation < 0) max(pt_df$pseudotime_raw) - pt_df$pseudotime_raw else pt_df$pseudotime_raw
pt_df$pseudotime <- as.numeric(scale(pt_df$pseudotime)) # z-score for interpretability

age_pseudotime_cor <- cor(pt_df$pseudotime, pt_df$age_days, method = "spearman")
cat("\nSpearman correlation, pseudotime vs. nominal age (sanity check, should be strongly positive): ",
    round(age_pseudotime_cor, 3), "\n", sep = "")

readr::write_csv(pt_df, file.path(root, "output", "tables", "manualgating_pseudotime_values.csv"))

# ---- 5a. Cross-sectional: pseudotime SynF vs CtrlF, per timepoint (mirrors Fig 6f) ----
# p_value (Wilcoxon) tests for a general stochastic/rank shift; p_value_median
# (Mood's median test, common.R) is the more literal test of whether the two
# groups' medians differ -- both reported since they answer different questions.
cross_sectional <- purrr::map_dfr(levels(pt_df$timepoint), function(tp) {
  d <- pt_df |> dplyr::filter(timepoint == tp)
  x <- d$pseudotime[d$group_feeding == "SynF"]
  y <- d$pseudotime[d$group_feeding == "CtrlF"]
  wt <- suppressWarnings(wilcox.test(x, y))
  mt <- moods_median_test(x, y)
  cd <- cohens_d_with_ci(x, y)
  tibble::tibble(
    timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
    median_synf = median(x), median_ctrlf = median(y),
    p_value = wt$p.value, p_value_median = mt$p_value, median_test_method = mt$method,
    cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"), p_fdr_median = p.adjust(p_value_median, method = "BH"))

# ---- 5b. Omnibus group:timepoint interaction (LRT) on pseudotime ----
d_full <- pt_df |> tidyr::drop_na(pseudotime, group_feeding, timepoint, group_delivery, subject_id)
full_model <- lme4::lmer(pseudotime ~ group_feeding * timepoint + group_delivery + (1 | subject_id), data = d_full, REML = FALSE)
reduced_model <- lme4::lmer(pseudotime ~ group_feeding + timepoint + group_delivery + (1 | subject_id), data = d_full, REML = FALSE)
a <- anova(reduced_model, full_model)
omnibus <- tibble::tibble(lrt_chisq = a$Chisq[2], lrt_df = a$Df[2], p_value = a$`Pr(>Chisq)`[2])

# ---- 5c. Change from baseline (V1): per-subject delta pseudotime, SynF vs CtrlF ----
# Same rationale as Fig5_manualgating_synf_ctrlf_stats.R -- V1 is the
# enrollment/baseline visit, before the synbiotic can have acted, so a raw
# cross-sectional group difference at V1 can't be a treatment effect. Testing
# each subject's own change from their V1 pseudotime isolates whatever
# happens AFTER baseline, and is immune to a baseline offset (real or
# batch-driven) by construction.
pt_wide <- pt_df |>
  dplyr::select(subject_id, group_feeding, timepoint, pseudotime) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = pseudotime, names_prefix = "pt_")

pseudotime_change <- purrr::map_dfr(c("V3", "V5"), function(tp) {
  tp_col <- paste0("pt_", tp)
  if (!tp_col %in% colnames(pt_wide)) return(NULL)
  d <- pt_wide |> dplyr::transmute(group_feeding, delta = .data[[tp_col]] - pt_V1) |> tidyr::drop_na()
  x <- d$delta[d$group_feeding == "SynF"]
  y <- d$delta[d$group_feeding == "CtrlF"]
  wt <- suppressWarnings(wilcox.test(x, y))
  cd <- cohens_d_with_ci(x, y)
  tibble::tibble(
    follow_up_timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
    median_delta_synf = median(x), median_delta_ctrlf = median(y),
    p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

readr::write_csv(pseudotime_change, file.path(root, "output", "tables", "manualgating_pseudotime_change_from_baseline.csv"))

readr::write_csv(
  dplyr::bind_rows(
    cross_sectional |> dplyr::mutate(test = "cross_sectional_wilcoxon"),
    omnibus |> dplyr::mutate(test = "omnibus_interaction_lrt", timepoint = "V1-V3-V5")
  ),
  file.path(root, "output", "tables", "manualgating_pseudotime_group_tests.csv")
)

# ---- 6. Figures ----
pca_plot_df <- pt_df |>
  dplyr::mutate(PC1 = pca$x[match(cytof_id, mat$cytof_id), 1], PC2 = pca$x[match(cytof_id, mat$cytof_id), 2])

p1 <- ggplot2::ggplot(pca_plot_df, ggplot2::aes(x = PC1, y = PC2, color = pseudotime)) +
  ggplot2::geom_point(size = 1.8, alpha = 0.8) +
  ggplot2::geom_path(data = as.data.frame(curve$s[order(curve$lambda), ]), ggplot2::aes(x = PC1, y = PC2), color = "black", inherit.aes = FALSE) +
  ggplot2::scale_color_viridis_c(name = "Pseudotime") +
  ggplot2::labs(title = "Manually-gated CyTOF pseudotime: colored by pseudotime") +
  ggplot2::theme_bw(base_size = 9)

p2 <- ggplot2::ggplot(pca_plot_df, ggplot2::aes(x = PC1, y = PC2, color = factor(age_days))) +
  ggplot2::geom_point(size = 1.8, alpha = 0.8) +
  ggplot2::scale_color_manual(
    values = c("#F2AF4AFF", "#C36377FF", "#1D457FFF"),
    labels = TIMEPOINT_LABELS[c("V1", "V3", "V5")],
    name = "Timepoint"
  ) +
  ggplot2::labs(title = "colored by actual age (baseline/2 months/4 months)") +
  ggplot2::theme_bw(base_size = 9)

save_pdf(p1, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_pca.pdf"), width = 5.5, height = 4.5)
save_pdf(p2, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_pca_age.pdf"), width = 5.5, height = 4.5)

# ---- Horizontal raincloud: pseudotime ("immune maturation age") on x-axis,
# one row per timepoint x feeding-group combination on y. Built manually
# (half-violin polygon + boxplot + jittered points) since no raincloud
# package (ggrain/PupillometryR/ggdist/gghalves) is installed. Cohen's d +
# p-value from the cross-sectional test (computed above) annotated per
# timepoint.
row_positions <- tidyr::expand_grid(
  timepoint = factor(c("V1", "V3", "V5"), levels = c("V1", "V3", "V5")),
  group_feeding = factor(c("CtrlF", "SynF"), levels = c("CtrlF", "SynF"))
) |>
  dplyr::mutate(
    # Baseline (V1) at top, 4 months (V5) at bottom; within each timepoint,
    # CtrlF above SynF -- reading order top-to-bottom: Baseline (CtrlF, SynF),
    # 2 months (CtrlF, SynF), 4 months (CtrlF, SynF). Gaps sized so a row's
    # rain (extends down to y0-0.85) can't reach the next row's cloud
    # (extends up to y0+0.9): needs >1.75 clearance, so 1.8 within a
    # timepoint block and 4.2 between blocks (> 1.8+1.75) -- avoids the
    # overplotting seen at the previous 1.2/2.6 spacing.
    timepoint_block = 3 - as.integer(timepoint), # V1->2 (top), V3->1, V5->0 (bottom)
    group_offset = 2 - as.integer(group_feeding), # CtrlF->1 (top), SynF->0 (bottom)
    y0 = timepoint_block * 4.2 + group_offset * 1.8
  )

raincloud_df <- pt_df |>
  dplyr::inner_join(row_positions, by = c("timepoint", "group_feeding"))

set.seed(42) # reproducible jitter for the "rain" points
raincloud_df$y_rain <- raincloud_df$y0 - stats::runif(nrow(raincloud_df), 0.20, 0.85)

cloud_polys <- raincloud_df |>
  dplyr::group_by(timepoint, group_feeding, y0) |>
  dplyr::group_modify(function(d, key) {
    dens <- stats::density(d$pseudotime, n = 200)
    height <- 0.9 * dens$y / max(dens$y) # normalize each row's cloud to the same max height
    tibble::tibble(
      x = c(dens$x, rev(dens$x)),
      y = c(key$y0 + height, rep(key$y0, length(dens$x)))
    )
  }) |>
  dplyr::ungroup()

row_medians <- raincloud_df |>
  dplyr::group_by(timepoint, group_feeding, y0) |>
  dplyr::summarise(median_pt = stats::median(pseudotime), .groups = "drop")

stat_annotations <- cross_sectional |>
  dplyr::left_join(
    row_positions |> dplyr::distinct(timepoint) |>
      dplyr::mutate(y_mid = purrr::map_dbl(timepoint, function(tp) mean(row_positions$y0[row_positions$timepoint == tp]))),
    by = "timepoint"
  ) |>
  dplyr::mutate(label = sprintf("d=%.2f, p(Wilcoxon)=%.3f\np(median)=%.3f", cohens_d, p_value, p_value_median))

x_range <- range(raincloud_df$pseudotime)
x_annotation <- x_range[2] + 0.15 * diff(x_range)

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

p3 <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = cloud_polys, ggplot2::aes(x = x, y = y, group = interaction(timepoint, group_feeding), fill = group_feeding),
    alpha = 0.6, color = NA
  ) +
  ggplot2::geom_boxplot(
    data = raincloud_df,
    ggplot2::aes(x = pseudotime, y = y0 - 0.15, group = interaction(timepoint, group_feeding), fill = group_feeding),
    orientation = "y", width = 0.25, outlier.shape = NA, linewidth = 0.3
  ) +
  ggplot2::geom_point(
    data = raincloud_df,
    ggplot2::aes(x = pseudotime, y = y_rain, color = group_feeding),
    size = 0.9, alpha = 0.6
  ) +
  ggplot2::geom_segment(
    data = row_medians,
    ggplot2::aes(x = median_pt, xend = median_pt, y = y0, yend = y0 + 0.9),
    color = "black", linewidth = 0.5
  ) +
  ggplot2::geom_text(
    data = stat_annotations, ggplot2::aes(x = x_annotation, y = y_mid, label = label),
    hjust = 0, size = 2.6, color = "grey20"
  ) +
  ggplot2::scale_fill_manual(values = color_group, name = "Feeding Group") +
  ggplot2::scale_color_manual(values = color_group, guide = "none") +
  ggplot2::scale_y_continuous(
    breaks = row_positions$y0,
    labels = paste0(TIMEPOINT_LABELS[as.character(row_positions$timepoint)], " - ", row_positions$group_feeding)
  ) +
  ggplot2::coord_cartesian(xlim = c(x_range[1] - 0.1 * diff(x_range), x_annotation + 0.35 * diff(x_range)), clip = "off") +
  ggplot2::labs(
    title = "Manually-gated CyTOF: immune maturation age (pseudotime) by feeding group",
    x = "Immune maturation age (pseudotime, z-scored)", y = NULL
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
    plot.margin = ggplot2::margin(5.5, 45, 5.5, 5.5)
  )

save_pdf(p3, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_by_group.pdf"), width = 7, height = 7)

# ---- Change-from-baseline figure ----
pt_change_long <- pt_wide |>
  tidyr::pivot_longer(cols = dplyr::starts_with("pt_V") & !dplyr::any_of("pt_V1"), names_to = "timepoint", values_to = "pt_followup") |>
  dplyr::mutate(timepoint = sub("pt_", "", timepoint), delta = pt_followup - pt_V1) |>
  tidyr::drop_na(delta) |>
  dplyr::mutate(timepoint = factor(
    timepoint,
    levels = c("V3", "V5"),
    labels = paste0("Baseline -> ", TIMEPOINT_LABELS[c("V3", "V5")])
  ))

p4 <- ggplot2::ggplot(pt_change_long, ggplot2::aes(x = timepoint, y = delta, fill = group_feeding)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  ggplot2::geom_boxplot(outlier.size = 0.5, position = ggplot2::position_dodge(width = 0.75), width = 0.6) +
  ggplot2::scale_fill_manual(values = c(CtrlF = "#39AE71", SynF = "#33AEFA"), name = "Feeding Group") +
  ggplot2::labs(
    title = "Manually-gated CyTOF pseudotime: change from baseline",
    x = NULL, y = "Change in pseudotime (z-scored, follow-up - baseline)"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))

save_pdf(p4, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_change_from_baseline.pdf"), width = 5, height = 4)

cat("\nCross-sectional pseudotime SynF vs CtrlF:\n")
print(cross_sectional)
cat("\nOmnibus interaction:\n")
print(omnibus)
cat("\nChange from baseline (V1):\n")
print(pseudotime_change)
