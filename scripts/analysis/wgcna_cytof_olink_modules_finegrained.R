# FINE-GRAINED VARIANT of wgcna_cytof_olink_modules.R -- per Petter's
# follow-up, 2026-09-06: the primary run (minClusterSize=15, merge
# cutHeight=0.25) found only 3 large modules, none differing by feeding
# group -- but that coarse tree cut could be absorbing smaller, more
# specific sub-modules (particularly within the large 208-feature
# "turquoise" catch-all) into the dominant age-driven structure. This variant
# uses a smaller minClusterSize (6) and deeper split (deepSplit=4), plus much
# lighter module merging (cutHeight=0.10, vs 0.25 in the primary run), to see
# whether any smaller module shows even a nominal trend for a feeding-group
# difference that the coarser cut washed out. Everything else (data,
# network construction, soft power) is identical to the primary run.
#
# Co-regulation module analysis (WGCNA) across manually-gated CyTOF
# populations AND Olink NPX proteins together, per Petter's request,
# 2026-09-06 -- the follow-up check flagged in the SynF/CtrlF interpretation
# summary: do the nominally differentially-abundant cytokines (KLRB1, CD160,
# TNFSF10, IL12RB1, FASLG, GZMA, CD70, CTSO, LY9, IL2RB, LAMP3, LGALS9) form
# one co-regulated biological axis (ideally alongside a matching CyTOF cell
# population), or are they independent coincidental nominal hits?
#
# WGCNA (Langfelder & Horvath 2008, BMC Bioinformatics, PMID 19114008) was
# chosen as it is already a vetted method in the lab (see
# Claude/verified-references.md) and is purpose-built for exactly this --
# unsupervised co-expression/co-regulation MODULE detection from a
# correlation network, rather than a simple correlation heatmap.
#
# IMPORTANT CAVEAT, stated up front: samples are pooled across all three
# timepoints (V1/V3/V5) and across subjects (repeated measures, NOT
# independent observations) to get enough samples (n<=248) for a stable
# correlation/adjacency structure. Age is already known to be the dominant
# axis of variation in both modalities (Fig5_olink_pca_age.R: 214/358
# proteins FDR-significant vs timepoint; Fig6 pseudotime). Consequently, the
# modules found here are expected to be substantially driven by shared
# age-related co-drift, not necessarily direct regulatory relationships --
# module-trait correlations against age_days are reported explicitly so this
# is visible rather than hidden, and should be read as descriptive/
# hypothesis-generating, not causal.
#
# Method:
#   1. Combine cytof_manual_clean (27 populations, CLR-transformed -- same
#      compositional-data rationale as every other CyTOF ordination script in
#      this repo, since raw percentages are not free-standing for
#      correlation) with olink_npx_clean (358 proteins, already log2-scale
#      NPX), joined per subject_id x timepoint (n=248 samples have both).
#   2. WGCNA signed network: biweight midcorrelation (robust to outliers),
#      soft-thresholding power chosen via pickSoftThreshold (targeting
#      scale-free topology R^2 ~ 0.8-0.9), topological overlap matrix (TOM),
#      average-linkage hierarchical clustering, dynamic tree cut
#      (minModuleSize=15, appropriate for a ~385-feature network vs WGCNA's
#      usual genome-wide scale), close-module merging (dissimilarity<0.25).
#   3. Module eigengenes (first PC per module) correlated against age_days,
#      feeding group (SynF=1/CtrlF=0), and timepoint (numeric V1/V3/V5).
#
# Output:
#   output/tables/wgcna_module_assignments.csv   (feature, module, module membership kME)
#   output/tables/wgcna_module_trait_correlations.csv
#   output/figures/manuscript/wgcna_dendrogram_modules.pdf
#   output/figures/manuscript/wgcna_module_trait_heatmap.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
  library(WGCNA)
})
WGCNA::disableWGCNAThreads() # avoid parallel-cluster overhead for this sample size

load_required_packages(c("dplyr", "tidyr", "tibble", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# ---- 1. Combine CyTOF (CLR) + Olink (log2 NPX) into one per-sample feature matrix ----
cytof_meta <- base$metadata |>
  dplyr::distinct(cytof_id, subject_id, timepoint, group_feeding) |>
  dplyr::mutate(cytof_id = as.character(cytof_id)) |>
  dplyr::filter(!is.na(cytof_id))

cytof_mat <- as.data.frame(base$cytof_manual_clean)
rownames(cytof_mat) <- cytof_mat$cytof_id
cytof_mat$cytof_id <- NULL
cytof_clr <- clr_transform(as.matrix(cytof_mat)) # common.R -- same Aitchison-space rationale used throughout this repo
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
rownames(combined) <- paste(link$subject_id, link$timepoint, sep = "__")

sample_traits <- link |>
  dplyr::mutate(
    row_id = paste(subject_id, timepoint, sep = "__"),
    age_days = TIMEPOINT_DAYS[as.character(timepoint)],
    feeding_numeric = ifelse(group_feeding == "EF", 1, 0), # EF=SynF=1, CF=CtrlF=0
    timepoint_numeric = as.integer(factor(timepoint, levels = c("V1", "V3", "V5")))
  ) |>
  dplyr::select(row_id, age_days, feeding_numeric, timepoint_numeric)
stopifnot(!any(duplicated(sample_traits$row_id)))
traits <- sample_traits |> tibble::column_to_rownames("row_id")
traits <- traits[rownames(combined), ]

cat("Combined matrix:", nrow(combined), "samples x", ncol(combined), "features (",
    sum(grepl("^CyTOF_", colnames(combined))), "CyTOF +", sum(grepl("^Olink_", colnames(combined))), "Olink)\n")

# ---- 2. WGCNA: soft-thresholding power ----
powers <- c(1:10, seq(12, 20, 2))
sft <- WGCNA::pickSoftThreshold(combined, powerVector = powers, corFnc = "bicor", networkType = "signed", verbose = 0)
soft_power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq > 0.8)[1]]
if (is.na(soft_power)) soft_power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
cat("\nSoft-thresholding power chosen:", soft_power,
    "(scale-free R^2 =", round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == soft_power], 3), ")\n")

adjacency_mat <- WGCNA::adjacency(combined, power = soft_power, corFnc = "bicor", type = "signed")
tom <- WGCNA::TOMsimilarity(adjacency_mat, TOMType = "signed", verbose = 0)
dissTOM <- 1 - tom
rownames(dissTOM) <- colnames(dissTOM) <- colnames(combined)

gene_tree <- fastcluster::hclust(as.dist(dissTOM), method = "average")
dynamic_mods <- dynamicTreeCut::cutreeDynamic(
  dendro = gene_tree, distM = dissTOM, method = "hybrid",
  deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 6
)
module_colors <- WGCNA::labels2colors(dynamic_mods)

# Merge modules whose eigengenes are highly correlated (dissimilarity < 0.10,
# much lighter than the primary run's 0.25 -- deliberately preserves more
# distinct small modules instead of collapsing them into the big age-driven ones)
merge_result <- WGCNA::mergeCloseModules(combined, module_colors, cutHeight = 0.10, corFnc = "bicor", verbose = 0)
module_colors_merged <- merge_result$colors
module_eigengenes <- merge_result$newMEs

n_modules <- length(unique(module_colors_merged)) - as.integer("grey" %in% module_colors_merged)
cat("\n", n_modules, "modules found (after merging close modules); 'grey' = unassigned.\n", sep = "")
cat("Module sizes:\n")
print(table(module_colors_merged))

# Numeric module labels (largest module = 1, etc.), used in place of WGCNA's
# arbitrary color names for the dendrogram figure only -- per Petter's
# request, 2026-09-08. "grey" (unassigned features) gets no number. The
# underlying module identity elsewhere in this script and its downstream
# consumers (assignment/correlation CSVs, other WGCNA scripts) is unchanged
# and still keyed by color name.
module_size_rank <- sort(table(module_colors_merged[module_colors_merged != "grey"]), decreasing = TRUE)
module_number_map <- setNames(seq_along(module_size_rank), names(module_size_rank))
cat("\nNumeric module labels (dendrogram figure only), by size:\n")
print(module_number_map)

# ---- 3. Module-trait correlations ----
module_trait_cor <- WGCNA::cor(module_eigengenes, traits, use = "pairwise.complete.obs")
module_trait_p <- WGCNA::corPvalueStudent(module_trait_cor, nrow(combined))

module_trait_df <- as.data.frame(module_trait_cor) |>
  tibble::rownames_to_column("module") |>
  tidyr::pivot_longer(-module, names_to = "trait", values_to = "correlation") |>
  dplyr::left_join(
    as.data.frame(module_trait_p) |> tibble::rownames_to_column("module") |>
      tidyr::pivot_longer(-module, names_to = "trait", values_to = "p_value"),
    by = c("module", "trait")
  ) |>
  dplyr::arrange(trait, dplyr::desc(abs(correlation)))

readr::write_csv(module_trait_df, file.path(root, "output", "tables", "wgcna_finegrained_module_trait_correlations.csv"))

module_sizes <- table(module_colors_merged)
cat("\nModules ranked by |correlation| with feeding group (looking for any nominal trend):\n")
print(
  module_trait_df |>
    dplyr::filter(trait == "feeding_numeric") |>
    dplyr::mutate(module = sub("^ME", "", module), module_size = as.integer(module_sizes[module])) |>
    dplyr::arrange(p_value) |>
    dplyr::select(module, module_size, correlation, p_value)
)

# ---- Module membership (kME) + assignment table, with cytokine/CyTOF-of-interest flagged ----
kme <- WGCNA::signedKME(combined, module_eigengenes, corFnc = "bicor")
assignment <- tibble::tibble(
  feature = colnames(combined),
  modality = ifelse(grepl("^CyTOF_", colnames(combined)), "CyTOF", "Olink"),
  feature_name = sub("^(CyTOF_|Olink_)", "", colnames(combined)),
  module = module_colors_merged
) |>
  dplyr::mutate(kME = kme[cbind(seq_along(module), match(paste0("kME", module), colnames(kme)))])

readr::write_csv(assignment, file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"))

# Proteins previously flagged in the SynF/CtrlF interpretation (cross-sectional
# V5 + change-from-baseline top hits) -- report which module(s) they land in.
FLAGGED_PROTEINS <- c("KLRB1", "CD160", "TNFSF10", "IL12RB1", "FASLG", "GZMA", "CD70", "CTSO", "LY9", "IL2RB", "LAMP3", "LGALS9")
cat("\nModule assignment of previously-flagged SynF/CtrlF proteins:\n")
print(assignment |> dplyr::filter(feature_name %in% FLAGGED_PROTEINS) |> dplyr::arrange(module))

flagged_modules <- assignment |> dplyr::filter(feature_name %in% FLAGGED_PROTEINS) |> dplyr::pull(module) |> unique()
cat("\nAll features (CyTOF + Olink) in the same module(s) as the flagged proteins:\n")
print(assignment |> dplyr::filter(module %in% flagged_modules) |> dplyr::arrange(module, modality, dplyr::desc(kME)))

# ---- Figures ----
# Dendrogram + module band, with modules identified by NUMBER (largest = 1)
# instead of WGCNA's arbitrary color names -- per Petter's request,
# 2026-09-08. Contiguous same-module runs along the leaf order (as cut by
# the dendrogram) are shaded in alternating neutral greys purely to separate
# blocks visually; the shade itself carries no meaning -- the module number,
# printed centered in each run, is the label. Unassigned ("grey"/no module)
# runs are left unlabeled.
leaf_order <- gene_tree$order
n_leaves <- length(leaf_order)
leaf_modules <- module_colors_merged[leaf_order]
run_id <- cumsum(c(TRUE, leaf_modules[-1] != leaf_modules[-n_leaves]))
runs <- stats::aggregate(seq_len(n_leaves), by = list(run_id = run_id, module = leaf_modules), FUN = function(x) c(min(x), max(x)))
run_bounds <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
  data.frame(module = runs$module[i], start = runs$x[i, 1], end = runs$x[i, 2])
}))
run_bounds <- run_bounds[order(run_bounds$start), ]
run_bounds$number <- ifelse(run_bounds$module == "grey", NA, module_number_map[run_bounds$module])
run_bounds$shade <- rep(c("grey88", "grey96"), length.out = nrow(run_bounds))
run_bounds$shade[is.na(run_bounds$number)] <- "white"
# PAM reassignment (dynamicTreeCut, method="hybrid") does not guarantee
# contiguous module blocks along the leaf order -- a handful of leaves can
# be interleaved outliers. Suppress the number label (background shading is
# kept) on any run narrower than 1% of all leaves so those don't overlap
# into illegible clutter; only the dominant contiguous block(s) per module
# are labeled.
min_run_width <- max(1, round(0.01 * n_leaves))
run_bounds$label <- ifelse(is.na(run_bounds$number) | (run_bounds$end - run_bounds$start + 1) < min_run_width,
                            "", run_bounds$number)

pdf(file.path(root, "output", "figures", "manuscript", "wgcna_finegrained_dendrogram_modules.pdf"), width = 10, height = 5.5)
graphics::layout(matrix(1:2, nrow = 2), heights = c(4, 1))
graphics::par(mar = c(0, 4, 2, 1))
plot(gene_tree, labels = FALSE, hang = 0.03, xlab = "", sub = "", ylab = "Height",
     main = "CyTOF + Olink co-regulation modules (fine-grained)")
graphics::par(mar = c(1, 4, 0, 1))
plot(1, type = "n", xlim = c(0.5, n_leaves + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "Module")
graphics::rect(run_bounds$start - 0.5, 0, run_bounds$end + 0.5, 1, col = run_bounds$shade, border = "white")
graphics::text(
  x = (run_bounds$start + run_bounds$end) / 2, y = 0.5,
  labels = run_bounds$label,
  cex = 0.75, font = 2
)
# Highlight Module 3's block (the module examined in detail in Fig5 panels
# f/g) with a black border -- per Petter's manual edit, 2026-09-08. Uses
# the single widest run for module 3, matching the largest contiguous block
# that carries its printed number label.
brown_run <- run_bounds[run_bounds$module == "brown", ]
brown_run <- brown_run[which.max(brown_run$end - brown_run$start), ]
graphics::rect(brown_run$start - 0.5, 0, brown_run$end + 0.5, 1, col = NA, border = "black", lwd = 2)
dev.off()

trait_labels <- c(age_days = "Age (days)", feeding_numeric = "Feeding (SynF=1)", timepoint_numeric = "Timepoint (1-3)")
text_matrix <- paste0(round(module_trait_cor, 2), "\n(", signif(module_trait_p, 2), ")")
dim(text_matrix) <- dim(module_trait_cor)

pdf(file.path(root, "output", "figures", "manuscript", "wgcna_finegrained_module_trait_heatmap.pdf"), width = 6, height = max(4, 0.35 * nrow(module_trait_cor)))
WGCNA::labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = trait_labels[colnames(module_trait_cor)],
  yLabels = rownames(module_trait_cor),
  ySymbols = rownames(module_trait_cor),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = TRUE,
  cex.text = 0.6,
  main = "Module-trait relationships, fine-grained (correlation, p-value)"
)
dev.off()

cat("\nDone. See output/tables/wgcna_finegrained_module_assignments.csv and wgcna_finegrained_module_trait_correlations.csv.\n")
