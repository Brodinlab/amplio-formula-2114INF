# MDS of Aitchison distances between samples (manually-gated CyTOF populations),
# colored by age/timepoint -- does overall cell composition shift with age?
# Same method as the existing FlowSOM-based Fig5B_mds_panelB_timepoint.R, applied
# to base$cytof_manual instead of base$cytof_lineage, for direct comparison.

suppressPackageStartupMessages({
  # NOTE: see Fig5_manualgating_synf_ctrlf_stats.R -- code/lib/ vs scripts/lib/
  # path mismatch upstream, using the real path here.
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "robCompositions"))

root <- get_repo_root()
base <- load_base_tables(root)

cytof_manual <- base$cytof_manual
stopifnot("cytof_id" %in% colnames(cytof_manual))
# NOTE: cytof_manual is a tibble (from readr::read_csv); tibbles silently
# refuse `rownames<-` (deprecation warning, no-op), which would leave the
# distance matrix rownames as bare integers "1","2",... and lose the
# cytof_id -> row mapping entirely -- converting to a plain data.frame first
# is required for rownames() to actually take effect. The existing
# Fig5B_mds_panelB_timepoint.R / Fig5C_mds_panelA_feeding.R scripts in this
# repo have the same `rownames(mat) <- as.character(mat$cytof_id)` pattern on
# a tibble and are very likely affected the same way on this R/tibble version
# -- flagging, not fixed here (out of scope for this script).
mat <- as.data.frame(cytof_manual)
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- mat + 0.001 # pseudocount, avoids log(0) in the Aitchison distance

aitch <- robCompositions::aDist(as.matrix(mat))
aitch[aitch < 0] <- 0
aitch <- (aitch + t(aitch)) / 2
aitch[is.na(aitch) | is.infinite(aitch) | is.nan(aitch)] <- max(aitch[is.finite(aitch)], na.rm = TRUE)

mds <- stats::cmdscale(aitch, k = 2)
df <- as.data.frame(mds) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(
    base$metadata |> dplyr::select(cytof_id, timepoint),
    by = "cytof_id"
  ) |>
  dplyr::rename(MDS1 = V1, MDS2 = V2) |>
  dplyr::mutate(timepoint = factor(timepoint, levels = c("V1", "V3", "V5")))

color_timepoint <- c("V1" = "#F2AF4AFF", "V3" = "#C36377FF", "V5" = "#1D457FFF")

p <- ggplot2::ggplot(df, ggplot2::aes(x = MDS1, y = MDS2, color = timepoint)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::stat_ellipse(level = 0.95, alpha = 0.3) +
  ggplot2::scale_color_manual(values = color_timepoint, name = "Timepoint") +
  ggplot2::labs(
    title = "Manually-gated CyTOF: Aitchison-distance MDS, colored by age",
    x = "MDS1", y = "MDS2"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"))

save_pdf(p, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_mds_age.pdf"), width = 5.5, height = 4.5)

cat("n samples plotted:", nrow(df), "(", sum(is.na(df$timepoint)), "missing timepoint)\n")
