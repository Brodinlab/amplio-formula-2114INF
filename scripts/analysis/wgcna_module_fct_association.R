# Tests whether the fine-grained WGCNA co-regulation modules
# (wgcna_cytof_olink_modules_finegrained.R) relate to fecal community type
# (FCT) trajectory group (T1-T4, metadata$group_FCT), per Petter's request,
# 2026-09-06 -- the recommended next step from the SynF/CtrlF interpretation:
# connecting the CD8 effector/memory module finding to the manuscript's own
# central microbiome-trajectory narrative (Figure 6a/b), rather than leaving
# it as an isolated immune-only result.
#
# NOTE (2026-09-08, per Petter): originally run against the primary (coarse,
# 3-module) set, which has since been removed -- rerun here against the
# fine-grained module set, which is now the only one retained. The
# fine-grained "brown" module is a different (smaller, 40- vs 68-feature)
# reconstruction than the primary "brown" module, so its FCT association
# p-value is expected to differ from the earlier primary-set result.
#
# Context already established: FCT trajectory is itself strongly associated
# with feeding group at the subject level (group_feeding x group_FCT,
# metadata) -- SynF subjects are concentrated in T2 (24/38, 63%), CtrlF in
# T1+T4 (32/42, 76%). This script tests the next link in that chain: does
# the immune module (already shown to diverge by feeding group at 4 months)
# also differ across these same trajectory groups?
#
# Reuses the already-computed module ASSIGNMENTS
# (wgcna_finegrained_module_assignments.csv) rather than re-running the
# network construction -- eigengenes only need WGCNA::moduleEigengenes() on
# the same combined CyTOF-CLR/Olink-NPX matrix used throughout.
#
# Method: at V5 (4 months, the timepoint where the feeding-group effect was
# found), each module's eigengene is compared across the 4 FCT trajectories
# by Kruskal-Wallis omnibus test; where omnibus p<0.10, pairwise Wilcoxon
# comparisons (BH-FDR across pairs) follow, mirroring the manuscript's own
# stated approach to FCT-trajectory comparisons (Figure 6c: "reference
# trajectories" vs "alternative trajectories", effect sizes, adjusted
# p-values).
#
# Output:
#   output/tables/wgcna_module_fct_association.csv (omnibus KW per module)
#   output/tables/wgcna_module_fct_pairwise.csv (pairwise, where run)
#   output/figures/manuscript/wgcna_module_by_fct_boxplot.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
  library(WGCNA)
})
WGCNA::disableWGCNAThreads()

load_required_packages(c("dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# ---- Rebuild the same combined CyTOF(CLR) + Olink(NPX) matrix used throughout ----
cytof_meta <- base$metadata |>
  dplyr::distinct(cytof_id, subject_id, timepoint, group_feeding, group_FCT) |>
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
  dplyr::filter(cytof_id %in% rownames(cytof_clr), SampleID %in% rownames(olink_mat), !is.na(group_FCT))

combined <- cbind(
  cytof_clr[link$cytof_id, , drop = FALSE],
  as.matrix(olink_mat[link$SampleID, , drop = FALSE])
)
row_id <- paste(link$subject_id, link$timepoint, sep = "__")
rownames(combined) <- row_id

meta_df <- link |>
  dplyr::transmute(
    row_id = row_id, subject_id,
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5")),
    group_FCT = factor(group_FCT, levels = c("T1", "T2", "T3", "T4"))
  ) |>
  tibble::column_to_rownames("row_id")
meta_df <- meta_df[rownames(combined), ]

# ---- Module eigengenes (fine-grained set) ----
assignment <- readr::read_csv(file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"), show_col_types = FALSE)
colors <- assignment$module[match(colnames(combined), assignment$feature)]
stopifnot(!any(is.na(colors)))

MEs <- WGCNA::moduleEigengenes(combined, colors)$eigengenes
modules <- setdiff(colnames(MEs), "MEgrey")

eigen_df <- MEs |>
  dplyr::mutate(row_id = rownames(combined)) |>
  dplyr::left_join(meta_df |> tibble::rownames_to_column("row_id"), by = "row_id")

# ---- At V5: module eigengene across FCT trajectory (Kruskal-Wallis + pairwise) ----
v5_df <- eigen_df |> dplyr::filter(timepoint == "V5")

omnibus <- purrr::map_dfr(modules, function(mod) {
  d <- v5_df |> dplyr::select(value = dplyr::all_of(mod), group_FCT) |> tidyr::drop_na()
  kt <- kruskal.test(value ~ group_FCT, data = d)
  tibble::tibble(
    module = sub("^ME", "", mod), n = nrow(d),
    kruskal_chisq = unname(kt$statistic), kruskal_df = unname(kt$parameter), p_value = kt$p.value
  )
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(p_value)

readr::write_csv(omnibus, file.path(root, "output", "tables", "wgcna_module_fct_association.csv"))

cat("Module eigengene (V5) across FCT trajectory, Kruskal-Wallis omnibus:\n")
print(omnibus)

pairwise_all <- purrr::map_dfr(modules, function(mod) {
  mod_name <- sub("^ME", "", mod)
  omni_p <- omnibus$p_value[omnibus$module == mod_name]
  if (omni_p >= 0.10) return(NULL) # only follow up modules with at least a suggestive omnibus signal
  d <- v5_df |> dplyr::select(value = dplyr::all_of(mod), group_FCT) |> tidyr::drop_na()
  combos <- utils::combn(levels(d$group_FCT), 2, simplify = FALSE)
  purrr::map_dfr(combos, function(pair) {
    x <- d$value[d$group_FCT == pair[1]]
    y <- d$value[d$group_FCT == pair[2]]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    wt <- suppressWarnings(wilcox.test(x, y))
    cd <- cohens_d_with_ci(x, y)
    tibble::tibble(module = mod_name, trajectory_a = pair[1], trajectory_b = pair[2],
                   n_a = length(x), n_b = length(y), cohens_d = cd$d, p_value = wt$p.value)
  })
})

if (nrow(pairwise_all) > 0) {
  pairwise_all <- pairwise_all |>
    dplyr::group_by(module) |>
    dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()
  readr::write_csv(pairwise_all, file.path(root, "output", "tables", "wgcna_module_fct_pairwise.csv"))
  cat("\nPairwise trajectory comparisons (modules with omnibus p<0.10):\n")
  print(pairwise_all |> dplyr::arrange(module, p_value))
} else {
  cat("\nNo module reached the omnibus p<0.10 threshold for pairwise follow-up.\n")
}

# ---- Figure: module eigengene by FCT trajectory, V5 only, all fine-grained modules ----
plot_df <- v5_df |>
  tidyr::pivot_longer(dplyr::all_of(modules), names_to = "module", values_to = "eigengene") |>
  dplyr::mutate(module = sub("^ME", "", module)) |>
  dplyr::filter(!is.na(group_FCT))

color_fct <- c(T1 = "#33A190", T2 = "#A1D5CF", T3 = "#0F5096", T4 = "#1C79E3")

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group_FCT, y = eigengene, fill = group_FCT)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  ggplot2::facet_wrap(~module, scales = "free_y") +
  ggplot2::scale_fill_manual(values = color_fct, name = "FCT trajectory") +
  ggplot2::labs(x = "FCT trajectory", y = "Module eigengene (V5, 4 months)") +
  ggplot2::theme_bw(base_size = 9)

save_pdf(p, file.path(root, "output", "figures", "manuscript", "wgcna_module_by_fct_boxplot.pdf"), width = 9, height = 4)

cat("\nDone.\n")
