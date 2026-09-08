# Split-violin plot of the fine-grained WGCNA "brown" module eigengene
# (CD8 effector/memory axis: CD8 TCM/TEM/TEMRA + 37 Olink proteins including
# the originally-flagged KLRB1/CD160/TNFSF10/IL12RB1/FASLG/GZMA/LGALS9/LAMP3)
# by feeding group at each timepoint -- per Petter's request, 2026-09-08.
#
# Supersedes wgcna_finegrained_brown_raincloud.R (removed): that raincloud
# construction was visually identical to the existing pseudotime raincloud
# (Fig6_manualgating_pseudotime.R), so this uses a different chart type --
# one violin per timepoint, split in half (left = SynF, right = CtrlF),
# rather than a stacked cloud/box/jitter per group/timepoint row.
#
# No split-violin geom package is installed (gghalves/see/ggdist) -- built
# manually as mirrored density polygons (standard technique: geom_polygon
# with the density-estimate curve on one side of the category's x position
# and a flat edge on the other, mirrored for the second group), sharing one
# global width scale across all 6 half-violins (3 timepoints x 2 groups) so
# widths are comparable panel-to-panel. Nominal p-value/effect-size
# annotated per timepoint per the wiki's figure-style convention
# (observations/figure-style-conventions.md).
#
# Reuses the fine-grained module eigengene computation from
# wgcna_module_feeding_interaction.R / wgcna_finegrained_module_assignments.csv
# rather than the expensive adjacency/TOM step.
#
# Output:
#   output/figures/manuscript/wgcna_finegrained_brown_splitviolin.pdf

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
  library(WGCNA)
})
WGCNA::disableWGCNAThreads()

load_required_packages(c("dplyr", "tidyr", "purrr", "tibble", "readr", "ggplot2"))

root <- get_repo_root()
base <- load_base_tables(root)

# ---- Rebuild the same combined CyTOF(CLR) + Olink(NPX) matrix used by the WGCNA scripts ----
cytof_meta <- base$metadata |>
  dplyr::distinct(cytof_id, subject_id, timepoint, group_feeding) |>
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
    subject_id,
    group_feeding = factor(dplyr::recode(group_feeding, EF = "SynF", CF = "CtrlF"), levels = c("CtrlF", "SynF")),
    timepoint = factor(timepoint, levels = c("V1", "V3", "V5"))
  ) |>
  tibble::column_to_rownames("row_id")
meta_df <- meta_df[rownames(combined), ]

# ---- Fine-grained module eigengene, "brown" only ----
assignment <- readr::read_csv(
  file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"),
  show_col_types = FALSE
)
colors <- assignment$module[match(colnames(combined), assignment$feature)]
stopifnot(!any(is.na(colors)))

MEs <- WGCNA::moduleEigengenes(combined, colors)$eigengenes
stopifnot("MEbrown" %in% colnames(MEs))

eigen_df <- MEs |>
  dplyr::select(eigengene = MEbrown) |>
  dplyr::mutate(row_id = rownames(combined)) |>
  dplyr::left_join(meta_df |> tibble::rownames_to_column("row_id"), by = "row_id") |>
  dplyr::mutate(timepoint_x = as.integer(timepoint)) # V1=1, V3=2, V5=3

# ---- Cross-sectional SynF vs CtrlF per timepoint (for inline annotation) ----
cross_sectional <- purrr::map_dfr(levels(eigen_df$timepoint), function(tp) {
  d <- eigen_df |> dplyr::filter(timepoint == tp)
  x <- d$eigengene[d$group_feeding == "SynF"]
  y <- d$eigengene[d$group_feeding == "CtrlF"]
  wt <- suppressWarnings(wilcox.test(x, y))
  cd <- cohens_d_with_ci(x, y)
  tibble::tibble(timepoint = tp, p_value = wt$p.value, cohens_d = cd$d)
}) |>
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

cat("Fine-grained 'brown' module, SynF vs CtrlF per timepoint:\n")
print(cross_sectional)

# ---- Build mirrored half-violin polygons: left = SynF, right = CtrlF ----
HALF_WIDTH_MAX <- 0.42 # max horizontal extent of a half-violin from its timepoint's x position

densities <- eigen_df |>
  dplyr::group_by(timepoint, timepoint_x, group_feeding) |>
  dplyr::group_modify(function(d, key) {
    # Trimmed to the observed data range (from/to) -- otherwise density()'s
    # default tail padding beyond min/max creates a near-zero-width sliver
    # that, combined with coord_cartesian(clip = "off") for the annotation
    # text above the panel, renders as a stray vertical line artifact.
    dens <- stats::density(d$eigengene, n = 200, from = min(d$eigengene), to = max(d$eigengene))
    tibble::tibble(value = dens$x, density = dens$y)
  }) |>
  dplyr::ungroup()

# One shared width scale across all 6 half-violins, so widths are comparable
global_max_density <- max(densities$density)
densities$half_width <- HALF_WIDTH_MAX * densities$density / global_max_density

violin_polys <- densities |>
  dplyr::group_by(timepoint, timepoint_x, group_feeding) |>
  dplyr::group_modify(function(d, key) {
    d <- d[order(d$value), ]
    x0 <- key$timepoint_x
    if (key$group_feeding == "SynF") {
      # left half: density curve bulges left, flat edge on the center line
      poly_x <- c(x0 - d$half_width, rev(rep(x0, nrow(d))))
    } else {
      # right half: density curve bulges right, flat edge on the center line
      poly_x <- c(rep(x0, nrow(d)), rev(x0 + d$half_width))
    }
    poly_y <- c(d$value, rev(d$value))
    tibble::tibble(x = poly_x, y = poly_y)
  }) |>
  dplyr::ungroup()

medians <- eigen_df |>
  dplyr::group_by(timepoint, timepoint_x, group_feeding) |>
  dplyr::summarise(median_eigengene = stats::median(eigengene), .groups = "drop") |>
  dplyr::mutate(
    x_start = timepoint_x,
    x_end = ifelse(group_feeding == "SynF", timepoint_x - HALF_WIDTH_MAX * 0.7, timepoint_x + HALF_WIDTH_MAX * 0.7)
  )

stat_annotations <- cross_sectional |>
  dplyr::left_join(eigen_df |> dplyr::distinct(timepoint, timepoint_x), by = "timepoint") |>
  dplyr::mutate(
    # Simplified to d/p only (no p_fdr) to match the final assembled Fig5 --
    # per Petter's manual edit, 2026-09-08. Full p_fdr values remain in
    # wgcna_module_feeding_cross_sectional.csv for reference.
    label = sprintf("d=%.2f\np=%.3f", cohens_d, p_value),
    y = max(eigen_df$eigengene) + 0.15 * diff(range(eigen_df$eigengene))
  )

# Group labels shown as direct colored text above the Baseline pair only,
# replacing a color-key legend entirely -- per Petter's manual edit,
# 2026-09-08 (observations/figure-style-conventions.md: prefer direct
# labeling over a separate legend box where the mapping is unambiguous).
group_labels <- tibble::tibble(
  group_feeding = factor(c("SynF", "CtrlF"), levels = c("CtrlF", "SynF")),
  x = 1 + c(-1, 1) * HALF_WIDTH_MAX * 0.5,
  y = max(eigen_df$eigengene) + 0.05 * diff(range(eigen_df$eigengene))
)

# Numeric module label (largest fine-grained module = 1, etc.), matching
# the numbering convention used in the dendrogram figure
# (wgcna_cytof_olink_modules_finegrained.R) -- computed from the same
# assignment table so it stays correct if module sizes change.
module_sizes <- assignment |> dplyr::filter(module != "grey") |> dplyr::count(module, sort = TRUE)
brown_number <- which(module_sizes$module == "brown")

# Individual points, jittered within each half only (never crossing the
# center line into the other group's half) -- offset kept close to center
# (0.03-0.32 of HALF_WIDTH_MAX) so points sit over the violin body.
set.seed(42)
points_df <- eigen_df |>
  dplyr::mutate(
    jitter = stats::runif(dplyr::n(), 0.03, 0.32),
    x_point = ifelse(group_feeding == "SynF", timepoint_x - jitter, timepoint_x + jitter)
  )

color_group <- c(CtrlF = "#39AE71", SynF = "#33AEFA")

p <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = violin_polys, ggplot2::aes(x = x, y = y, group = interaction(timepoint, group_feeding), fill = group_feeding),
    alpha = 0.75, color = "grey30", linewidth = 0.2
  ) +
  ggplot2::geom_point(
    data = points_df, ggplot2::aes(x = x_point, y = eigengene, color = group_feeding),
    size = 0.9, alpha = 0.6
  ) +
  ggplot2::geom_segment(
    data = medians,
    ggplot2::aes(x = x_start, xend = x_end, y = median_eigengene, yend = median_eigengene),
    color = "black", linewidth = 0.8
  ) +
  ggplot2::geom_text(
    data = stat_annotations, ggplot2::aes(x = timepoint_x, y = y, label = label),
    size = 2.6, color = "grey20", lineheight = 0.9
  ) +
  ggplot2::geom_text(
    data = group_labels, ggplot2::aes(x = x, y = y, label = group_feeding, color = group_feeding),
    size = 3.6, fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = color_group, guide = "none") +
  ggplot2::scale_color_manual(values = color_group, guide = "none") +
  ggplot2::scale_x_continuous(breaks = 1:3, labels = TIMEPOINT_LABELS[c("V1", "V3", "V5")], limits = c(0.5, 3.5)) +
  ggplot2::coord_cartesian(ylim = c(min(eigen_df$eigengene), max(stat_annotations$y) + 0.05 * diff(range(eigen_df$eigengene))), clip = "off") +
  ggplot2::labs(x = NULL, y = paste0("Module ", brown_number, " eigengene")) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(plot.margin = ggplot2::margin(20, 5.5, 5.5, 5.5))

save_pdf(p, file.path(root, "output", "figures", "manuscript", "wgcna_finegrained_brown_splitviolin.pdf"), width = 6.5, height = 5.5)

cat("\nSaved: output/figures/manuscript/wgcna_finegrained_brown_splitviolin.pdf\n")
