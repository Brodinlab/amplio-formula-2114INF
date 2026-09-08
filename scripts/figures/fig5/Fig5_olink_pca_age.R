# Olink NPX: PCA (all proteins) colored by age, plus a formal per-protein
# age-association test. Replaces the broken Fig5D_olink_pca_timepoint.R
# (wrong source path: "code/lib/..." instead of "scripts/lib/...") and adds
# the rigor missing from it (no age-association test was previously done at
# the individual-protein level, only an ad hoc "top 20% by variance of
# timepoint means" ranking used downstream for pseudotime -- see
# Fig5_olink_top_age_biplot.R for that ranking's replacement).
#
# Method: age-association per protein is quantified as eta-squared from a
# one-way ANOVA of each protein's NPX value against time point (V1/V3/V5),
# computed across all samples irrespective of feeding group -- the same
# method used for the CyTOF manual-gating pseudotime population ranking
# (Fig6_manualgating_pseudotime.R), for consistency across modalities.
# Benjamini-Hochberg FDR is applied across all 358 proteins.
#
# Output:
#   output/tables/olink_age_association.csv   (per-protein eta^2, ANOVA p, FDR)
#   output/figures/manuscript/Fig5_olink_pca_age.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "purrr", "tidyr", "readr", "tibble", "ggplot2"))

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
  dplyr::inner_join(meta |> dplyr::select(SampleID, subject_id, group_feeding, timepoint, age_days), by = "SampleID")

# ---- 1. PCA on all proteins, colored by age/timepoint ----
pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(meta |> dplyr::select(SampleID, group_feeding, timepoint), by = "SampleID")

var_explained <- summary(pca)$importance["Proportion of Variance", 1:2] * 100

color_timepoint <- c(V1 = "#F2AF4AFF", V3 = "#C36377FF", V5 = "#1D457FFF")

p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::scale_color_manual(values = color_timepoint, labels = TIMEPOINT_LABELS[c("V1", "V3", "V5")], name = "Timepoint") +
  ggplot2::labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::coord_equal()

save_pdf(p_pca, file.path(root, "output", "figures", "manuscript", "Fig5_olink_pca_age.pdf"), width = 5.5, height = 4.5)

# ---- 2. Per-protein age-association: eta-squared, one-way ANOVA vs timepoint ----
age_association <- purrr::map_dfr(protein_cols, function(protein) {
  d <- df |> dplyr::select(value = dplyr::all_of(protein), timepoint) |> tidyr::drop_na()
  fit <- stats::aov(value ~ timepoint, data = d)
  ss <- summary(fit)[[1]][["Sum Sq"]]
  p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
  tibble::tibble(protein = protein, eta_squared = ss[1] / sum(ss), p_value = p_val)
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(dplyr::desc(eta_squared))

readr::write_csv(age_association, file.path(root, "output", "tables", "olink_age_association.csv"))

cat("PCA: PC1 explains", round(var_explained[1], 1), "%, PC2 explains", round(var_explained[2], 1), "% of variance\n")
cat("\nAge-association (eta^2 vs timepoint), top 15 proteins:\n")
print(age_association |> dplyr::slice_head(n = 15))
cat("\n", sum(age_association$p_fdr < 0.05), " of ", nrow(age_association),
    " proteins FDR-significant for age association; ",
    sum(age_association$p_value < 0.05), " nominal p<0.05\n", sep = "")
