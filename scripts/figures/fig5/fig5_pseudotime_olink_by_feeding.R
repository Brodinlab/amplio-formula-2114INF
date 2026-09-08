# Immune maturation pseudotime from Olink NPX data, mirroring
# Fig6_manualgating_pseudotime.R (CyTOF) and replacing the broken/ad hoc
# Fig6E_olink_topProteins_PCA_*.R and Fig6G_olink_delta_pseudotime_*.R
# (wrong source path; slingshot/SingleCellExperiment dependency; "top 20% by
# variance of timepoint means" ranking instead of a formal per-protein test).
#
# Protein selection: top age-associated proteins by eta-squared vs timepoint
# (already computed in Fig5_olink_pca_age.R -> output/tables/olink_age_association.csv),
# NOT restricted to cytokines -- per Petter's instruction, 2026-09-06 -- using
# the top 20% (~72 of 358 proteins), matching the proportion already cited in
# the manuscript's Figure 6e caption ("top 20%, n=71") and the CyTOF
# pseudotime script's own selection fraction, for consistency across both.
#
# Method (identical to the CyTOF script -- see its header for full rationale):
#   1. Rank all proteins by eta-squared (one-way ANOVA vs timepoint, all
#      samples pooled, group-blind) -- reuses the already-computed ranking
#      rather than re-running 358 ANOVAs.
#   2. Take the top 20%, with a redundancy filter dropping any candidate
#      correlated with an already-selected protein at |r|>=0.9 (guards
#      against near-duplicate/co-regulated proteins destabilizing the
#      principal-curve fit, as happened for the CyTOF Neutrophils/
#      NonNeutrophils pair).
#   3. Scale, PCA, fit a principal curve through PC1-PC2 (single-lineage,
#      same algorithm Slingshot uses with no branching).
#   4. Orient pseudotime to increase with age (flip sign if anti-correlated); z-score.
#   5. Test SynF vs CtrlF: cross-sectional Wilcoxon + Mood's median test per
#      timepoint, omnibus group:timepoint interaction LRT, and change from
#      baseline (V1) by Wilcoxon.
#
# Output:
#   output/tables/olink_pseudotime_values.csv
#   output/tables/olink_pseudotime_group_tests.csv
#   output/tables/olink_pseudotime_change_from_baseline.csv
#   output/figures/manuscript/fig5_pseudotime_olink_pca.pdf
#   output/figures/manuscript/fig5_pseudotime_olink_pca_age.pdf
#   output/figures/manuscript/fig5_pseudotime_olink_by_group.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "lme4", "princurve", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

meta <- base$olink_metadata |>
  dplyr::mutate(
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    age_days = TIMEPOINT_DAYS[as.character(timepoint)]
  )

mat <- as.data.frame(base$olink_npx_clean)
rownames(mat) <- mat$SampleID
mat$SampleID <- NULL
protein_cols <- colnames(mat)

df <- mat |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::inner_join(meta |> dplyr::select(SampleID, subject_id, group_feeding, group_delivery, timepoint, age_days), by = "SampleID") |>
  tidyr::drop_na(group_feeding, timepoint)

# ---- 1. Top 20% age-associated proteins (eta^2, reused from Fig5_olink_pca_age.R) ----
age_association <- readr::read_csv(file.path(root, "output", "tables", "olink_age_association.csv"), show_col_types = FALSE)
n_top <- max(2, round(0.20 * length(protein_cols)))

REDUNDANCY_THRESHOLD <- 0.9
ranked_proteins <- age_association |> dplyr::arrange(dplyr::desc(eta_squared)) |> dplyr::pull(protein)
top_proteins <- character(0)
for (protein in ranked_proteins) {
  if (length(top_proteins) >= n_top) break
  is_redundant <- length(top_proteins) > 0 && any(purrr::map_lgl(top_proteins, function(p) {
    abs(cor(df[[protein]], df[[p]], use = "complete.obs")) >= REDUNDANCY_THRESHOLD
  }))
  if (!is_redundant) top_proteins <- c(top_proteins, protein)
}

cat("Top", n_top, "age-varying proteins (eta^2 vs timepoint, redundancy-filtered):\n")
print(age_association |> dplyr::filter(protein %in% top_proteins) |> dplyr::arrange(dplyr::desc(eta_squared)))

# ---- 2. PCA on the top time-varying proteins, all samples pooled ----
pca_input <- scale(as.matrix(df |> dplyr::select(dplyr::all_of(top_proteins))))
pca <- prcomp(pca_input, center = FALSE, scale. = FALSE) # already scaled above

# ---- 3. Principal curve through PC1-PC2 (single-lineage pseudotime, as in Slingshot) ----
curve <- princurve::principal_curve(pca$x[, 1:2])
pseudotime_raw <- curve$lambda

pt_df <- tibble::tibble(SampleID = df$SampleID, pseudotime_raw = pseudotime_raw) |>
  dplyr::left_join(df |> dplyr::select(SampleID, subject_id, group_feeding, group_delivery, timepoint, age_days), by = "SampleID")

# ---- 4. Orient pseudotime to increase with age ----
orientation <- sign(cor(pt_df$pseudotime_raw, pt_df$age_days, use = "complete.obs", method = "spearman"))
if (orientation == 0) orientation <- 1
pt_df$pseudotime <- if (orientation < 0) max(pt_df$pseudotime_raw) - pt_df$pseudotime_raw else pt_df$pseudotime_raw
pt_df$pseudotime <- as.numeric(scale(pt_df$pseudotime)) # z-score for interpretability

age_pseudotime_cor <- cor(pt_df$pseudotime, pt_df$age_days, method = "spearman")
cat("\nSpearman correlation, pseudotime vs. nominal age (sanity check, should be strongly positive): ",
    round(age_pseudotime_cor, 3), "\n", sep = "")

readr::write_csv(pt_df, file.path(root, "output", "tables", "olink_pseudotime_values.csv"))

# ---- 5a. Cross-sectional: pseudotime SynF vs CtrlF, per timepoint ----
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

readr::write_csv(pseudotime_change, file.path(root, "output", "tables", "olink_pseudotime_change_from_baseline.csv"))
readr::write_csv(
  dplyr::bind_rows(
    cross_sectional |> dplyr::mutate(test = "cross_sectional_wilcoxon"),
    omnibus |> dplyr::mutate(test = "omnibus_interaction_lrt", timepoint = "V1-V3-V5")
  ),
  file.path(root, "output", "tables", "olink_pseudotime_group_tests.csv")
)

# ---- 6. Figures ----
pca_plot_df <- pt_df |>
  dplyr::mutate(PC1 = pca$x[match(SampleID, df$SampleID), 1], PC2 = pca$x[match(SampleID, df$SampleID), 2])

p1 <- ggplot2::ggplot(pca_plot_df, ggplot2::aes(x = PC1, y = PC2, color = pseudotime)) +
  ggplot2::geom_point(size = 1.8, alpha = 0.8) +
  ggplot2::geom_path(data = as.data.frame(curve$s[order(curve$lambda), ]), ggplot2::aes(x = PC1, y = PC2), color = "black", inherit.aes = FALSE) +
  ggplot2::scale_color_viridis_c(name = "Pseudotime") +
  ggplot2::labs(title = "Olink pseudotime: colored by pseudotime") +
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

save_pdf(p1, file.path(root, "output", "figures", "manuscript", "fig5_pseudotime_olink_pca.pdf"), width = 5.5, height = 4.5)
save_pdf(p2, file.path(root, "output", "figures", "manuscript", "fig5_pseudotime_olink_pca_age.pdf"), width = 5.5, height = 4.5)

# ---- Horizontal raincloud: pseudotime on x-axis, one row per timepoint x
# feeding-group combination on y (same construction as Fig6_manualgating_pseudotime.R:
# half-violin polygon + boxplot + jittered points, since no raincloud
# package is installed). Cohen's d + p-values from the cross-sectional test
# annotated per timepoint.
row_positions <- tidyr::expand_grid(
  timepoint = factor(c("V1", "V3", "V5"), levels = c("V1", "V3", "V5")),
  group_feeding = factor(c("CtrlF", "SynF"), levels = c("CtrlF", "SynF"))
) |>
  dplyr::mutate(
    timepoint_block = 3 - as.integer(timepoint),
    group_offset = 2 - as.integer(group_feeding),
    y0 = timepoint_block * 4.2 + group_offset * 1.8
  )

raincloud_df <- pt_df |>
  dplyr::inner_join(row_positions, by = c("timepoint", "group_feeding"))

set.seed(42)
raincloud_df$y_rain <- raincloud_df$y0 - stats::runif(nrow(raincloud_df), 0.20, 0.85)

cloud_polys <- raincloud_df |>
  dplyr::group_by(timepoint, group_feeding, y0) |>
  dplyr::group_modify(function(d, key) {
    dens <- stats::density(d$pseudotime, n = 200)
    height <- 0.9 * dens$y / max(dens$y)
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
  ggplot2::labs(x = "Immune maturation age (pseudotime, z-scored)", y = NULL) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(plot.margin = ggplot2::margin(5.5, 45, 5.5, 5.5))

save_pdf(p3, file.path(root, "output", "figures", "manuscript", "fig5_pseudotime_olink_by_group.pdf"), width = 7, height = 7)

cat("\nCross-sectional pseudotime SynF vs CtrlF:\n")
print(cross_sectional)
cat("\nOmnibus interaction:\n")
print(omnibus)
cat("\nChange from baseline (V1):\n")
print(pseudotime_change)
