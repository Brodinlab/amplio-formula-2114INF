# MDS of Aitchison distances, plate/batch effect regressed out first.
#
# Fig5_manualgating_mds_age.R showed no visible age structure, but a strong
# plate_4-vs-rest cluster (2026-09-03 investigation) that plate is subject-level
# (86/86 subjects entirely on one plate, none split -- so within-subject
# contrasts are naturally robust, but the raw MDS pools samples across subjects
# and plate, so the batch effect dominates the visible structure there).
#
# Aitchison distance == Euclidean distance in CLR (centered log-ratio) space, so
# "regress out plate" is done by: CLR-transform (clr_transform(), common.R),
# fit value ~ plate per population and take residuals (+ global mean, to keep
# values interpretable), then Euclidean distance + classical MDS on the
# corrected CLR matrix. This is a CLR-space correction, distinct from (but
# consistent in spirit with) the percentage-space correction baked into
# cytof_manual_clean used for the univariate stats scripts -- the two need
# different treatments because Aitchison/compositional geometry only makes
# linear operations valid in log-ratio space.
#
# Outlier excluded here too (subject PHL001-0073, V1, CtrlF, cytof_id
# 453612193 -- see scripts/export/export_cytof_manual_clean.R), per Petter's
# instruction 2026-09-03.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble"))

root <- get_repo_root()
base <- load_base_tables(root)

OUTLIER_CYTOF_ID <- "453612193"

mat <- base$cytof_manual |>
  dplyr::filter(as.character(cytof_id) != OUTLIER_CYTOF_ID) |>
  as.data.frame()
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)

clr_mat <- clr_transform(mat) # samples x populations, log-ratio space

meta_lookup <- base$metadata |>
  dplyr::distinct(cytof_id, cytof_plate, timepoint, group_feeding) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))
plate_vec <- meta_lookup$cytof_plate[match(rownames(clr_mat), meta_lookup$cytof_id)]
stopifnot(!any(is.na(plate_vec)))

# Regress out plate per population, keep residuals + global mean (preserves scale)
clr_corrected <- apply(clr_mat, 2, function(col) {
  fit <- stats::lm(col ~ factor(plate_vec))
  stats::residuals(fit) + mean(col)
})
rownames(clr_corrected) <- rownames(clr_mat)

euclid <- stats::dist(clr_corrected)
mds <- stats::cmdscale(as.matrix(euclid), k = 2)

df <- as.data.frame(mds) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(meta_lookup, by = "cytof_id") |>
  dplyr::rename(MDS1 = V1, MDS2 = V2) |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("V1", "V3", "V5")))

color_timepoint <- c("V1" = "#F2AF4AFF", "V3" = "#C36377FF", "V5" = "#1D457FFF")

p_age <- ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::stat_ellipse(level = 0.95, alpha = 0.3) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Manually-gated CyTOF (outlier removed, plate regressed out): colored by age",
    x = "MDS1", y = "MDS2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

p_plate_check <- ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = plate_vec)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::labs(title = "Same, colored by plate (sanity check -- should now be mixed)", color = "Plate") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p_age, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_mds_age_clean.pdf"), width = 5.5, height = 4.5)
save_pdf(p_plate_check, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_mds_plate_clean_check.pdf"), width = 6, height = 4.5)

cat("n samples:", nrow(df), "\n")
