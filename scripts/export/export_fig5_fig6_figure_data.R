suppressPackageStartupMessages({
  source("code/lib/common.R")
  source("code/lib/load_data.R")
})

load_required_packages(c(
  "dplyr", "tibble", "readr", "tidyr", "ggrepel",
  "robCompositions", "OlinkAnalyze"
))

root <- get_repo_root()
dat <- load_amplio_data(root)

out_root <- file.path(root, "data/figure_data")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## -----------------------
## Fig5B/Fig5C: MDS panels
## Source: scripts/amplio_merge.Rmd chunk `mds-plots-dimred`
## Fig5C = panel A (feeding), Fig5B = panel B (timepoint)
## -----------------------
{
  fig_dir <- file.path(out_root, "fig5")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  cytof_lineage_data_mds <- dat$cytof_lineage_data
  cytof_lineage_data_mds[is.na(cytof_lineage_data_mds)] <- 0
  cytof_lineage_data_mds <- cytof_lineage_data_mds + 0.001

  aitch <- robCompositions::aDist(cytof_lineage_data_mds)
  aitch[aitch < 0] <- 0
  aitch <- (aitch + t(aitch)) / 2
  aitch[is.na(aitch) | is.infinite(aitch) | is.nan(aitch)] <- max(aitch[is.finite(aitch)], na.rm = TRUE)

  mds <- stats::cmdscale(aitch, k = 2)
  mds_df <- as.data.frame(mds) |>
    tibble::rownames_to_column("cytof_id") |>
    dplyr::left_join(dat$cytof_metadata, by = "cytof_id") |>
    dplyr::rename(MDS1 = V1, MDS2 = V2)

  readr::write_csv(
    mds_df |> dplyr::select(cytof_id, MDS1, MDS2, subject_id, group_feeding, timepoint, group_delivery),
    file.path(fig_dir, "fig5bc_mds_lineage.csv")
  )
}

## -----------------------
## Fig5E + Fig5D (derived): Olink PCA
## Targets:
## - Fig5E: Extra_Olink_PCA_by_feeding_group.pdf (color by feeding)
## - Fig5D: same PCA but color by timepoint
## -----------------------
{
  fig_dir <- file.path(out_root, "fig5")
  olink_x <- dat$olink_data_wide
  olink_meta <- dat$olink_metadata

  x <- olink_x
  x[is.na(x)] <- 0
  pca <- stats::prcomp(x, scale. = TRUE, center = TRUE)

  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("olink_id") |>
    dplyr::left_join(olink_meta, by = "olink_id") |>
    dplyr::rename(PC1 = PC1, PC2 = PC2)

  readr::write_csv(
    pca_df |>
      dplyr::select(olink_id, PC1 = PC1, PC2 = PC2, subject_id, group_feeding, timepoint, group_FCT, group_delivery),
    file.path(fig_dir, "fig5de_olink_pca_points.csv")
  )
}

## -----------------------
## Fig5F: Volcano top hits (Figure_4.13_volcano_top_hits.pdf)
## -----------------------
{
  fig_dir <- file.path(out_root, "fig5")
  timepoints <- c("V1", "V3", "V5")
  protein_names <- colnames(dat$olink_data_wide)

  volcano_all <- list()
  for (tp in timepoints) {
    rows <- list()
    tp_mask <- dat$olink_metadata$timepoint == tp
    tp_meta <- dat$olink_metadata[tp_mask, , drop = FALSE]

    for (protein in protein_names) {
      vals <- dat$olink_data_wide[tp_mask, protein, drop = TRUE]
      cf <- vals[tp_meta$group_feeding == "CF"]; cf <- cf[!is.na(cf)]
      ef <- vals[tp_meta$group_feeding == "EF"]; ef <- ef[!is.na(ef)]
      if (length(cf) >= 3 && length(ef) >= 3) {
        log2FC <- mean(ef, na.rm = TRUE) - mean(cf, na.rm = TRUE)
        p <- stats::wilcox.test(cf, ef)$p.value
        rows[[length(rows) + 1]] <- data.frame(
          Timepoint = tp,
          Protein = protein,
          log2FC = log2FC,
          p_value = p,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) > 0) {
      d <- do.call(rbind, rows)
      d$neg_log10_p <- -log10(d$p_value)
      d$Category <- ifelse(d$p_value < 0.05, "p < 0.05", "Not significant")
      volcano_all[[tp]] <- d
    }
  }
  volcano_df <- do.call(rbind, volcano_all)
  readr::write_csv(volcano_df, file.path(fig_dir, "fig5f_olink_volcano_top_hits.csv"))
}

## -----------------------
## Fig6B: FCT trajectory distribution (Fig5E_left_fct_distribution.pdf)
## Source: fct-trajectory-distribution chunk
## -----------------------
{
  fig_dir <- file.path(out_root, "fig6")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  subjects_fct <- dat$metadata |>
    dplyr::distinct(subject_id, group_FCT) |>
    dplyr::filter(!is.na(group_FCT))

  plot_data <- subjects_fct |>
    dplyr::left_join(dat$metadata |> dplyr::distinct(subject_id, group_feeding), by = "subject_id") |>
    dplyr::filter(!is.na(group_FCT) & !is.na(group_feeding)) |>
    dplyr::group_by(group_FCT, group_feeding) |>
    dplyr::summarise(N = dplyr::n(), .groups = "drop")

  readr::write_csv(plot_data, file.path(fig_dir, "fig6b_fct_distribution.csv"))
}

## -----------------------
## Fig6D: CyTOF topClusters PCA (pseudotime + timepoint)
## Outputs match:
## - CyTOF_topClusters_PCA_pseudotime_nClusters.pdf
## - CyTOF_topClusters_PCA_timepoint_nClusters.pdf
## -----------------------
{
  fig_dir <- file.path(out_root, "fig6")
  top_prop <- 0.2
  time_var <- dat$cytof_cluster_data |>
    as.data.frame() |>
    tibble::rownames_to_column("cytof_id") |>
    tidyr::pivot_longer(cols = -cytof_id, names_to = "Cluster", values_to = "Frequency") |>
    dplyr::left_join(dat$cytof_metadata |> dplyr::select(cytof_id, timepoint), by = "cytof_id") |>
    dplyr::filter(!is.na(Frequency) & !is.na(timepoint)) |>
    dplyr::group_by(Cluster, timepoint) |>
    dplyr::summarise(mean_freq = mean(Frequency, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(Cluster) |>
    dplyr::summarise(time_var = stats::var(mean_freq, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(time_var))

  top_clusters <- time_var |>
    dplyr::slice_head(prop = top_prop) |>
    dplyr::pull(Cluster)
  top_clusters <- intersect(top_clusters, colnames(dat$cytof_cluster_data))

  x <- dat$cytof_cluster_data[, top_clusters, drop = FALSE]
  x[is.na(x)] <- 0
  x_clr <- clr_transform(x)
  pca <- stats::prcomp(x_clr, scale. = TRUE, center = TRUE)

  # Use PC1 as a deterministic proxy for "pseudotime" in the exported plot-ready table.
  # (The original report uses slingshot on the same PCA space; exporting slingshot values
  # would require running slingshot here and pinning package versions.)
  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("cytof_id") |>
    dplyr::left_join(dat$cytof_metadata, by = "cytof_id") |>
    dplyr::rename(PC1 = PC1, PC2 = PC2) |>
    dplyr::mutate(pseudotime_proxy = pca$x[, 1])

  readr::write_csv(
    pca_df |>
      dplyr::select(cytof_id, PC1, PC2, pseudotime_proxy, subject_id, timepoint, group_FCT, group_feeding),
    file.path(fig_dir, "fig6d_cytof_topclusters_pca_points.csv")
  )
}

## -----------------------
## Fig6E: Olink topProteins PCA (pseudotime + timepoint)
## Uses the "Extra_Olink_topProteins_PCA_*" outputs
## -----------------------
{
  fig_dir <- file.path(out_root, "fig6")
  top_prop <- 0.2

  olink_x <- dat$olink_data_wide
  olink_meta <- dat$olink_metadata

  prot_time_var <- as.data.frame(olink_x) |>
    tibble::rownames_to_column("olink_id") |>
    tidyr::pivot_longer(cols = -olink_id, names_to = "Protein", values_to = "NPX") |>
    dplyr::left_join(olink_meta |> dplyr::select(olink_id, timepoint), by = "olink_id") |>
    dplyr::filter(!is.na(NPX) & !is.na(timepoint)) |>
    dplyr::group_by(Protein, timepoint) |>
    dplyr::summarise(mean_npx = mean(NPX, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(Protein) |>
    dplyr::summarise(time_var = stats::var(mean_npx, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(time_var))

  top_proteins <- prot_time_var |>
    dplyr::slice_head(prop = top_prop) |>
    dplyr::pull(Protein)
  top_proteins <- intersect(top_proteins, colnames(olink_x))

  x <- olink_x[, top_proteins, drop = FALSE]
  x[is.na(x)] <- 0
  pca <- stats::prcomp(x, scale. = TRUE, center = TRUE)

  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("olink_id") |>
    dplyr::left_join(olink_meta, by = "olink_id") |>
    dplyr::rename(PC1 = PC1, PC2 = PC2) |>
    dplyr::mutate(pseudotime_proxy = pca$x[, 1])

  readr::write_csv(
    pca_df |>
      dplyr::select(olink_id, PC1, PC2, pseudotime_proxy, subject_id, timepoint, group_FCT, group_feeding),
    file.path(fig_dir, "fig6e_olink_topproteins_pca_points.csv")
  )
}

## -----------------------
## Fig6F/Fig6G: Δ-pseudotime density (faceted) for CyTOF / Olink
## Targets:
## - Extra_CyTOF_topClusters_delta_pseudotime_density_faceted_FCT.pdf
## - Extra_Olink_topProteins_delta_pseudotime_density_faceted_FCT.pdf
## -----------------------
{
  fig_dir <- file.path(out_root, "fig6")

  # CyTOF: reuse exported pseudotime proxy, compute Δ vs V1 per subject
  cytof_pt <- readr::read_csv(file.path(fig_dir, "fig6d_cytof_topclusters_pca_points.csv"), show_col_types = FALSE)
  cytof_delta <- cytof_pt |>
    dplyr::filter(timepoint %in% c("V1", "V3", "V5"), !is.na(subject_id), !is.na(group_FCT)) |>
    dplyr::group_by(subject_id, group_FCT, timepoint) |>
    dplyr::summarise(pt = mean(pseudotime_proxy, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = timepoint, values_from = pt)

  if (all(c("V1", "V3", "V5") %in% colnames(cytof_delta))) {
    cytof_delta_long <- cytof_delta |>
      dplyr::mutate(dV3 = V3 - V1, dV5 = V5 - V1) |>
      dplyr::select(subject_id, group_FCT, dV3, dV5) |>
      tidyr::pivot_longer(cols = c(dV3, dV5), names_to = "Delta_tp", values_to = "Delta_pt") |>
      dplyr::mutate(
        Delta_tp = dplyr::recode(Delta_tp, dV3 = "V3 - V1", dV5 = "V5 - V1")
      )
    readr::write_csv(cytof_delta_long, file.path(fig_dir, "fig6f_cytof_delta_pseudotime_long.csv"))
  }

  # Olink
  olink_pt <- readr::read_csv(file.path(fig_dir, "fig6e_olink_topproteins_pca_points.csv"), show_col_types = FALSE)
  olink_delta <- olink_pt |>
    dplyr::filter(timepoint %in% c("V1", "V3", "V5"), !is.na(subject_id), !is.na(group_FCT)) |>
    dplyr::group_by(subject_id, group_FCT, timepoint) |>
    dplyr::summarise(pt = mean(pseudotime_proxy, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = timepoint, values_from = pt)

  if (all(c("V1", "V3", "V5") %in% colnames(olink_delta))) {
    olink_delta_long <- olink_delta |>
      dplyr::mutate(dV3 = V3 - V1, dV5 = V5 - V1) |>
      dplyr::select(subject_id, group_FCT, dV3, dV5) |>
      tidyr::pivot_longer(cols = c(dV3, dV5), names_to = "Delta_tp", values_to = "Delta_pt") |>
      dplyr::mutate(
        Delta_tp = dplyr::recode(Delta_tp, dV3 = "V3 - V1", dV5 = "V5 - V1")
      )
    readr::write_csv(olink_delta_long, file.path(fig_dir, "fig6g_olink_delta_pseudotime_long.csv"))
  }
}

message("Done exporting figure_data for Fig5/Fig6 into: ", out_root)

