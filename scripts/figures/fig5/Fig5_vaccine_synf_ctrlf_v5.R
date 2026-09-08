# Vaccine-specific IgG (Luminex, 36-antigen panel) vs. SynF/CtrlF feeding
# group, restricted to the latest timepoint (V5, 4 months), per Petter's
# request, 2026-09-06: "assess vaccine response data at latest timepoint to
# see if feeding groups are associated with any difference in systemic
# immune function." This matches the manuscript's own stated Methods
# exactly: "For vaccine-specific antibodies, cross-sectional comparisons
# between feeding groups were performed at 4 months (V5) using Wilcoxon
# rank-sum tests" -- unlike CyTOF/Olink, no longitudinal LME/omnibus-
# interaction framework is used here, matching the narrower method already
# established for this data type in the manuscript.
#
# All 36 non-control antigens are included (export_vaccine.R), matching the
# manuscript's literal "36 vaccine-related antigens" count -- see that
# script's header for the Hib-antigen-missing / CMV-EBV-RSV-HPV-included
# discrepancy versus the manuscript's prose antigen-category list.
#
# Statistics: Wilcoxon rank-sum test + Cohen's d per antigen at V5,
# Benjamini-Hochberg FDR across the 36 antigens (FDR<0.05 significant,
# nominal p<0.05 exploratory, per manuscript convention).
#
# Output:
#   output/tables/vaccine_synf_ctrlf_v5.csv
#   output/figures/manuscript/Fig5_vaccine_volcano_v5.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("dplyr", "tidyr", "purrr", "readr", "ggplot2", "ggrepel"))

root <- get_repo_root()
base <- load_base_tables(root)

df <- base$vaccine_igg |>
  dplyr::filter(timepoint == "V5") |>
  dplyr::mutate(group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")))

results <- df |>
  dplyr::group_by(antigen) |>
  dplyr::group_modify(function(d, key) {
    x <- d$z_score[d$group_feeding == "SynF"]
    y <- d$z_score[d$group_feeding == "CtrlF"]
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(
      n_synf = length(x), n_ctrlf = length(y),
      p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    category = classify_vaccine_antigen(antigen)
  ) |>
  dplyr::arrange(p_value)
stopifnot(!any(is.na(results$category)))

readr::write_csv(results, file.path(root, "output", "tables", "vaccine_synf_ctrlf_v5.csv"))

cat("\nNominal hits (p<0.05) by category:\n")
print(results |> dplyr::filter(p_value < 0.05) |> dplyr::count(category))

# ---- Volcano plot: Cohen's d vs nominal p, colored by antigen category
# (infant-schedule vaccine vs maternal/MMR-family vs non-vaccine pathogen) --
# these answer different biological questions and were easily conflated in a
# flat antigen list; top 10 labeled regardless of significance threshold ----
volcano_df <- results |> dplyr::mutate(neg_log10_p = -log10(p_value))

N_LABEL <- 10
to_label <- volcano_df |> dplyr::slice_min(p_value, n = N_LABEL)

category_colors <- c(
  "Infant-schedule vaccine" = "#1D457FFF",
  "Maternal/MMR-family" = "#C36377FF",
  "Non-vaccine pathogen" = "grey60"
)

p_volcano <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = cohens_d, y = neg_log10_p)) +
  ggplot2::geom_point(ggplot2::aes(color = category), size = 1.8, alpha = 0.8) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = to_label, ggplot2::aes(label = antigen), size = 2.4, max.overlaps = 20, segment.size = 0.2) +
  ggplot2::scale_color_manual(values = category_colors, name = "Antigen category") +
  ggplot2::labs(x = "Cohen's d (SynF vs CtrlF), IgG z-score at 4 months", y = expression(-log[10]~"(nominal p-value)")) +
  ggplot2::theme_bw(base_size = 9)

save_pdf(p_volcano, file.path(root, "output", "figures", "manuscript", "Fig5_vaccine_volcano_v5.pdf"), width = 7, height = 6)

cat("Vaccine IgG, SynF vs CtrlF at V5 (4 months):", sum(results$p_fdr < 0.05), "of", nrow(results),
    "antigens FDR-significant;", sum(results$p_value < 0.05), "nominal p<0.05\n")
cat("\nTop 10 by nominal p-value:\n")
print(to_label |> dplyr::select(antigen, n_synf, n_ctrlf, cohens_d, p_value, p_fdr))
