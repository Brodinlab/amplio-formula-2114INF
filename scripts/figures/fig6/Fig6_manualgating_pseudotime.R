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
#   output/tables/manualgating_pseudotime_values.csv        (per-sample pseudotime)
#   output/tables/manualgating_pseudotime_group_tests.csv   (cross-sectional + omnibus)
#   output/figures/manuscript/Fig6_manualgating_pseudotime_pca.pdf
#   output/figures/manuscript/Fig6_manualgating_pseudotime_by_group.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "lme4", "princurve", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

pop_cols <- setdiff(colnames(base$cytof_manual), "cytof_id")
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

df <- base$cytof_manual |>
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
top_populations <- names(sort(eta_squared, decreasing = TRUE))[seq_len(n_top)]

cat("Top", n_top, "time-varying populations (eta^2 vs timepoint):\n")
print(round(sort(eta_squared, decreasing = TRUE)[seq_len(n_top)], 3))

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
cross_sectional <- purrr::map_dfr(levels(pt_df$timepoint), function(tp) {
  d <- pt_df |> dplyr::filter(timepoint == tp)
  x <- d$pseudotime[d$group_feeding == "SynF"]
  y <- d$pseudotime[d$group_feeding == "CtrlF"]
  wt <- suppressWarnings(wilcox.test(x, y))
  cd <- cohens_d_with_ci(x, y)
  tibble::tibble(
    timepoint = tp, n_synf = length(x), n_ctrlf = length(y),
    median_synf = median(x), median_ctrlf = median(y),
    p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

# ---- 5b. Omnibus group:timepoint interaction (LRT) on pseudotime ----
d_full <- pt_df |> tidyr::drop_na(pseudotime, group_feeding, timepoint, group_delivery, subject_id)
full_model <- lme4::lmer(pseudotime ~ group_feeding * timepoint + group_delivery + (1 | subject_id), data = d_full, REML = FALSE)
reduced_model <- lme4::lmer(pseudotime ~ group_feeding + timepoint + group_delivery + (1 | subject_id), data = d_full, REML = FALSE)
a <- anova(reduced_model, full_model)
omnibus <- tibble::tibble(lrt_chisq = a$Chisq[2], lrt_df = a$Df[2], p_value = a$`Pr(>Chisq)`[2])

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
  ggplot2::scale_color_manual(values = c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), name = "Timepoint (days)") +
  ggplot2::labs(title = "colored by actual age (V1/V3/V5)") +
  ggplot2::theme_bw(base_size = 9)

save_pdf(p1, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_pca.pdf"), width = 5.5, height = 4.5)
save_pdf(p2, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_pca_age.pdf"), width = 5.5, height = 4.5)

p3 <- ggplot2::ggplot(pt_df, ggplot2::aes(x = timepoint, y = pseudotime, fill = group_feeding)) +
  ggplot2::geom_boxplot(outlier.size = 0.5, position = ggplot2::position_dodge(width = 0.75), width = 0.6) +
  ggplot2::scale_fill_manual(values = c(CtrlF = "#39AE71", SynF = "#33AEFA"), name = "Feeding Group") +
  ggplot2::labs(title = "Manually-gated CyTOF pseudotime by feeding group", x = "Timepoint", y = "Pseudotime (z-scored)") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))

save_pdf(p3, file.path(root, "output", "figures", "manuscript", "Fig6_manualgating_pseudotime_by_group.pdf"), width = 5, height = 4)

cat("\nCross-sectional pseudotime SynF vs CtrlF:\n")
print(cross_sectional)
cat("\nOmnibus interaction:\n")
print(omnibus)
