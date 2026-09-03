# Diagnostic: same raw (uncorrected) Aitchison-distance MDS as
# Fig5_manualgating_mds_age.R, colored by cytof_plate and cytof_shipment
# instead of timepoint -- this is what identified the plate_4-vs-rest batch
# cluster on 2026-09-03, which motivated Fig5_manualgating_mds_age_batchcorrected.R.

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
  source("scripts/lib/load_base_tables.R")
})

load_required_packages(c("ggplot2", "dplyr", "tibble", "robCompositions"))

root <- get_repo_root()
base <- load_base_tables(root)

mat <- as.data.frame(base$cytof_manual)
rownames(mat) <- as.character(mat$cytof_id)
mat$cytof_id <- NULL
mat[is.na(mat)] <- 0
mat <- mat + 0.001

aitch <- robCompositions::aDist(as.matrix(mat))
aitch[aitch < 0] <- 0
aitch <- (aitch + t(aitch)) / 2
mds <- stats::cmdscale(aitch, k = 2)

df <- as.data.frame(mds) |>
  tibble::rownames_to_column("cytof_id") |>
  dplyr::left_join(
    base$metadata |> dplyr::select(cytof_id, cytof_plate, cytof_shipment),
    by = "cytof_id"
  ) |>
  dplyr::rename(MDS1 = V1, MDS2 = V2)

p_plate <- ggplot2::ggplot(df, ggplot2::aes(MDS1, MDS2, color = cytof_plate)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::labs(title = "Manually-gated CyTOF MDS (uncorrected), colored by plate", color = "Plate") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

p_shipment <- ggplot2::ggplot(df, ggplot2::aes(MDS1, MDS2, color = cytof_shipment)) +
  ggplot2::geom_point(size = 2, alpha = 0.7) +
  ggplot2::labs(title = "Manually-gated CyTOF MDS (uncorrected), colored by shipment", color = "Shipment") +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 9, face = "bold"))

save_pdf(p_plate, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_mds_plate_check.pdf"), width = 6, height = 4.5)
save_pdf(p_shipment, file.path(root, "output", "figures", "manuscript", "Fig5_manualgating_mds_shipment_check.pdf"), width = 6, height = 4.5)
