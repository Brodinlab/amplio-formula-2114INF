# Same as Fig5_manualgating_mds_age.R (raw Aitchison-distance MDS, colored by
# age/timepoint, NOT plate-corrected), with just the outlier sample excluded
# (subject PHL001-0073, V1, CtrlF, cytof_id 453612193) -- isolates how much of
# the original picture (no visible age gradient, dominated by the plate_4
# cluster) was driven by that one point vs. the plate effect itself. For the
# fully corrected version (outlier removed AND plate regressed out), see
# Fig5_manualgating_mds_age_batchcorrected.R.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "robCompositions"))

root <- get_repo_root()
base <- load_base_tables(root)

OUTLIER_CYTOF_ID <- "453612193"

cytof_manual <- base$cytof_manual |>
  dplyr::filter(as.character(cytof_id) != OUTLIER_CYTOF_ID)

mat <- as.data.frame(cytof_manual)
rownames(mat) <- as.character(mat$cytof_id) # see Fig5_manualgating_mds_age.R for the tibble-rownames caveat
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- mat + 0.001

aitch <- robCompositions::aDist(as.matrix(mat))
aitch[aitch < 0] <- 0
aitch <- (aitch + t(aitch)) / 2
aitch[is.na(aitch) | is.infinite(aitch) | is.nan(aitch)] <- max(aitch[is.finite(aitch)], na.rm = TRUE)

mds <- stats::cmdscale(aitch, k = 2)
df <- as.data.frame(mds) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(base$metadata |> dplyr::select(cytof_id, timepoint), by = "cytof_id") |>
  dplyr::rename(MDS1 = V1, MDS2 = V2) |>
  dplyr::mutate(timepoint = relabel_timepoint(timepoint))

color_timepoint <- setNames(c("#F2AF4AFF", "#C36377FF", "#1D457FFF"), TIMEPOINT_LABELS)

p <- ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::stat_ellipse(level = 0.95, alpha = 0.3) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Manually-gated CyTOF (outlier excluded, NOT plate-corrected): colored by age",
    x = "MDS1", y = "MDS2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "qc", "Fig5_manualgating_mds_age_outlierremoved.pdf"), width = 5.5, height = 4.5)

cat("n samples plotted:", nrow(df), "\n")
