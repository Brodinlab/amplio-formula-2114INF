# Tests whether fine-grained WGCNA module eigengenes
# (wgcna_cytof_olink_modules_finegrained.R) show a feeding-group-by-timepoint
# INTERACTION, per Petter's request, 2026-09-06: the module-trait
# correlations already run test feeding group as a flat main effect
# (ignoring time), which is the wrong test given feeding exposure is ~0 at
# baseline (V1, both groups pre-intervention) and only diverges by 2/4
# months -- a flat correlation across all three timepoints could dilute or
# miss an effect that is genuinely absent at V1 and present later.
#
# NOTE (2026-09-08, per Petter): the primary (coarse, 3-module) WGCNA run
# and its dependent scripts/outputs have been removed -- only the
# fine-grained module set is retained going forward.
#
# Reuses the already-computed module ASSIGNMENTS
# (wgcna_finegrained_module_assignments.csv) rather than re-running the
# expensive adjacency/TOM computation -- module eigengenes (first PC per
# module) only need the assignment + the same combined CyTOF-CLR/Olink-NPX
# matrix, via WGCNA::moduleEigengenes(), which is cheap.
#
# For each fine-grained module (excluding grey/unassigned):
#   - Omnibus test: LRT comparing eigengene ~ group_feeding * timepoint +
#     group_delivery + (1|subject_id) against the same model without the
#     interaction term -- mirrors the omnibus-interaction approach already
#     used for individual proteins/populations (Fig5_olink_synf_ctrlf_stats.R,
#     Fig5_manualgating_synf_ctrlf_stats.R).
#   - Cross-sectional Wilcoxon + Cohen's d per timepoint (V1/V3/V5), to show
#     directly whether V1 is null (as expected) while V3/V5 diverge.
#
# Output:
#   output/tables/wgcna_module_feeding_interaction.csv (omnibus LRT per module)
#   output/tables/wgcna_module_feeding_cross_sectional.csv (per-timepoint Wilcoxon/d)
#   output/figures/manuscript/wgcna_module_trajectories_by_group.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
  library(WGCNA)
})
WGCNA::disableWGCNAThreads()

load_required_packages(c("dplyr", "tidyr", "purrr", "tibble", "readr", "lme4", "lmerTest", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# ---- Rebuild the same combined CyTOF(CLR) + Olink(NPX) matrix used by both WGCNA scripts ----
cytof_meta <- base$metadata |>
  dplyr::distinct(cytof_id, subject_id, timepoint, group_feeding, group_delivery) |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::filter(!is.na(cytof_id))

cytof_mat <- as.data.frame(base$cytof_manual_clean)
rownames(cytof_mat) <- cytof_mat$cytof_id
cytof_mat$cytof_id <- NULL
cytof_clr <- clr_transform(as.matrix(cytof_mat))
colnames(cytof_clr) <- paste0("CyTOF_", colnames(cytof_clr))

olink_mat <- as.data.frame(base$olink_npx_clean)
rownames(olink_mat) <- olink_mat$SampleID
olink_mat$SampleID <- NULL
colnames(olink_mat) <- paste0("Olink_", colnames(olink_mat))

link <- cytof_meta |>
  dplyr::inner_join(
    base$olink_metadata |> dplyr::distinct(SampleID, subject_id, timepoint),
    by = c("subject_id", "timepoint")
  ) |>
  dplyr::filter(cytof_id %in% rownames(cytof_clr), SampleID %in% rownames(olink_mat))

combined <- cbind(
  cytof_clr[link$cytof_id, , drop = FALSE],
  as.matrix(olink_mat[link$SampleID, , drop = FALSE])
)
row_id <- paste(link$subject_id, link$timepoint, sep = "__")
rownames(combined) <- row_id

meta_df <- link |>
  dplyr::transmute(
    row_id = row_id,
    subject_id, group_delivery,
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  ) |>
  tibble::column_to_rownames("row_id")
meta_df <- meta_df[rownames(combined), ]

# ---- Compute fine-grained module eigengenes and run both tests ----
assignment <- readr::read_csv(file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"), show_col_types = FALSE)
colors <- assignment$module[match(colnames(combined), assignment$feature)]
stopifnot(!any(is.na(colors)))

MEs <- WGCNA::moduleEigengenes(combined, colors)$eigengenes
modules <- setdiff(colnames(MEs), "MEgrey")

eigen_df <- MEs |>
  tibble::rownames_to_column("row_id") |>
  dplyr::mutate(row_id = rownames(combined)) |> # moduleEigengenes drops rownames; combined's order is preserved
  dplyr::left_join(meta_df |> tibble::rownames_to_column("row_id"), by = "row_id")

omnibus_all <- purrr::map_dfr(modules, function(mod) {
  d <- eigen_df |> dplyr::select(value = dplyr::all_of(mod), group_feeding, timepoint, group_delivery, subject_id) |> tidyr::drop_na()
  full_model <- lme4::lmer(value ~ group_feeding * timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE)
  reduced_model <- lme4::lmer(value ~ group_feeding + timepoint + group_delivery + (1 | subject_id), data = d, REML = FALSE)
  a <- anova(reduced_model, full_model)
  tibble::tibble(module = sub("^ME", "", mod), n_features = sum(colors == sub("^ME", "", mod)),
                 lrt_chisq = a$Chisq[2], lrt_df = a$Df[2], p_value = a$`Pr(>Chisq)`[2])
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

cross_sectional_all <- purrr::map_dfr(modules, function(mod) {
  purrr::map_dfr(levels(eigen_df$timepoint), function(tp) {
    d <- eigen_df |> dplyr::filter(timepoint == tp)
    x <- d[[mod]][d$group_feeding == "SynF"]
    y <- d[[mod]][d$group_feeding == "CtrlF"]
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(module = sub("^ME", "", mod), timepoint = tp,
                   n_synf = length(x), n_ctrlf = length(y),
                   p_value = wt$p.value, cohens_d = cd$d, d_ci_lower = cd$ci_lower, d_ci_upper = cd$ci_upper)
  })
}) |>
  dplyr::group_by(timepoint) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup()

readr::write_csv(omnibus_all, file.path(root, "output", "tables", "wgcna_module_feeding_interaction.csv"))
readr::write_csv(cross_sectional_all, file.path(root, "output", "tables", "wgcna_module_feeding_cross_sectional.csv"))

cat("Omnibus feeding x timepoint interaction, fine-grained modules, ranked by p-value:\n")
print(omnibus_all)

cat("\nCross-sectional SynF vs CtrlF per module per timepoint (expect V1 ~ null):\n")
print(cross_sectional_all |> dplyr::arrange(module, timepoint))

# ---- Figure: module eigengene trajectories by feeding group, all fine-grained modules ----
traj_df <- eigen_df |>
  tidyr::pivot_longer(dplyr::all_of(modules), names_to = "module", values_to = "eigengene") |>
  dplyr::mutate(module = sub("^ME", "", module), age_days = TIMEPOINT_DAYS[as.character(timepoint)])

medians <- traj_df |>
  dplyr::group_by(module, group_feeding, age_days) |>
  dplyr::summarise(eigengene = median(eigengene, na.rm = TRUE), .groups = "drop")

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

p <- ggplot2::ggplot(traj_df, ggplot2::aes(x = age_days, y = eigengene, color = group_feeding)) +
  ggplot2::geom_line(ggplot2::aes(group = interaction(subject_id, group_feeding)), linewidth = 0.2, alpha = 0.15) +
  ggplot2::geom_line(data = medians, ggplot2::aes(group = group_feeding), linewidth = 1.1) +
  ggplot2::geom_point(data = medians, size = 1.8) +
  ggplot2::facet_wrap(~module, scales = "free_y") +
  ggplot2::scale_color_manual(values = color_group, name = "Feeding Group") +
  ggplot2::scale_x_continuous(breaks = TIMEPOINT_DAYS, labels = TIMEPOINT_LABELS) +
  ggplot2::labs(x = "Timepoint", y = "Module eigengene") +
  ggplot2::theme_bw(base_size = 8)

save_pdf(p, file.path(root, "output", "figures", "manuscript", "wgcna_module_trajectories_by_group.pdf"), width = 11, height = 7)

cat("\nDone.\n")
