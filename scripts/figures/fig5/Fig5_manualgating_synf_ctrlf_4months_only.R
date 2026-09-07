# Manually-gated CyTOF (cytof_manual_gating_frequency_clean.csv, Kanth's v1.1
# re-gated pass) vs. SynF/CtrlF feeding group, AT 4 MONTHS (V5) ONLY.
#
# Requested by Petter, 2026-09-07: a dedicated cross-sectional comparison
# restricted to the 4-month timepoint, rather than reading it off the
# per-timepoint-family results in Fig5_manualgating_synf_ctrlf_stats.R (where
# BH-FDR is computed within each of V1/V3/V5 separately). Same input data,
# same QC (3-outlier exclusion baked into the "_clean" table, 5 non-informative
# intermediate/complement populations dropped by load_base_tables()), same
# per-population test (Wilcoxon rank-sum + Cohen's d with CI) as that script's
# cross-sectional arm -- isolated here as its own single-timepoint family so
# BH-FDR is computed across the 27 populations at V5 alone, not pooled with
# V1/V3.
#
# Output: output/tables/manualgating_synf_ctrlf_4months_only.csv
#         output/figures/manuscript/Fig5_manualgating_4months_only_effect_sizes.pdf

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
    group_feeding = factor(
      dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"),
      levels = c("CtrlF", "SynF")
    )
  )

df <- base$cytof_manual_clean |>
  dplyr::inner_join(
    meta |> dplyr::select(cytof_id, subject_id, group_feeding, timepoint),
    by = "cytof_id"
  ) |>
  dplyr::filter(timepoint == "V5")

cat("4-month (V5) samples after QC: n_synf =", sum(df$group_feeding == "SynF"),
    ", n_ctrlf =", sum(df$group_feeding == "CtrlF"), "\n")

results <- purrr::map_dfr(pop_cols, function(pop) {
  d <- df |> dplyr::select(group_feeding, value = dplyr::all_of(pop)) |> tidyr::drop_na()
  x <- d$value[d$group_feeding == "SynF"]
  y <- d$value[d$group_feeding == "CtrlF"]
  if (length(x) < 3 || length(y) < 3) {
    return(tibble::tibble(
      population = pop, n_synf = length(x), n_ctrlf = length(y),
      median_synf = NA_real_, median_ctrlf = NA_real_,
      p_value = NA_real_, cohens_d = NA_real_, d_ci_lower = NA_real_, d_ci_upper = NA_real_
    ))
  }
  wt <- suppressWarnings(wilcox.test(x, y))
  cd <- cohens_d_with_ci(x, y)
  tibble::tibble(
    population = pop, n_synf = length(x), n_ctrlf = length(y),
    median_synf = median(x), median_ctrlf = median(y),
    p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

out_table_dir <- file.path(root, "output", "tables")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(results, file.path(out_table_dir, "manualgating_synf_ctrlf_4months_only.csv"))

cat("\n4-month-only cross-sectional results (", nrow(results), " populations):\n", sep = "")
cat("  Nominal p<0.05: ", sum(results$p_value < 0.05, na.rm = TRUE), "\n", sep = "")
cat("  BH-FDR<0.05:    ", sum(results$p_fdr < 0.05, na.rm = TRUE), "\n\n", sep = "")
print(results |> dplyr::filter(p_value < 0.05) |>
        dplyr::select(population, n_synf, n_ctrlf, median_synf, median_ctrlf, cohens_d, p_value, p_fdr))

# ---- Forest plot, ranked by effect size ----
plot_df <- results |>
  dplyr::filter(!is.na(cohens_d)) |>
  dplyr::arrange(cohens_d) |>
  dplyr::mutate(population = factor(population, levels = population))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cohens_d, y = population, color = p_value < 0.05)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = d_ci_lower, xmax = d_ci_upper), height = 0.2) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c(`TRUE` = "#C36377FF", `FALSE` = "grey40"), name = "nominal p<0.05") +
  ggplot2::labs(
    x = "Cohen's d (SynF vs CtrlF)", y = NULL,
    title = "Manually-gated CyTOF at 4 months only: SynF vs CtrlF"
  ) +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_4months_only_effect_sizes.pdf"),
         width = 7, height = 7)
cat("\nSaved: output/tables/manualgating_synf_ctrlf_4months_only.csv\n")
cat("Saved: output/figures/manuscript/Fig5_manualgating_4months_only_effect_sizes.pdf\n")
