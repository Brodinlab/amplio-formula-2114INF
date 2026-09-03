# Harmonizes AMPLIO's manually-gated CyTOF data (cytof_manual_gating_frequency_clean.csv:
# outlier excluded, plate regressed out per population -- see
# export_cytof_manual_clean.R) to the lineage-level schema used by the lab's
# internal cross-cohort longitudinal cell-frequency database (PAWS
# expData.lineageFreq / expData.cytof: Sweden/Iceland/Uganda/Kenya CyTOF +
# spectral flow, see Brodin wiki datasets/cell-frequency-data.md).
#
# Column mapping (AMPLIO fine-grained gate -> lineage-schema column):
#   B.cells                                              -> B_cells
#   CD4T                                                 -> CD4T
#   CD8T                                                 -> CD8T
#   TCRgd.cells                                          -> gdT
#   Tregs                                                -> Tregs
#   NK                                                   -> NK_cells
#   Neutrophils                                          -> Neutrophils
#   Eosinophils                                          -> Eosinophils
#   Basophils                                            -> Basophils
#   Classical + Intermediate + Nonclassical monocytes    -> Monocytes (summed)
#   mDCs                                                 -> mDC
#   pDCs                                                 -> pDC
# Not measured in the AMPLIO panel (left NA): Platelets, MSC.
# Not mapped (AMPLIO-specific finer hierarchy, no lineage-schema equivalent):
#   Neu_Eosino, NonNeutro_Eosino, NonTB, CD56bright/dim NK subsets, NonNK,
#   CD14neg.CD16neg, CD123pos, CD4T/CD8T memory subsets, HLADRpos.T, T.cells
#   (redundant with CD4T+CD8T+TCRgd.cells+HLADRpos.T) -- these remain
#   available in the untouched cytof_manual_gating_frequency_clean.csv for
#   any analysis that wants the full 32-population resolution.
#
# Output: data/tables/cross_cohort_lineage_harmonized.csv
# One row per AMPLIO CyTOF sample, in the same column layout as the PAWS
# expData.lineageFreq table, for appending to (not replacing) that
# cross-cohort table. Written for projects/regional-immune-development
# Paper 2 (see Brodin wiki, 2026-09-03).

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
})

load_required_packages(c("dplyr", "readr"))

root <- get_repo_root()

cytof <- readr::read_csv(file.path(root, "data", "tables", "cytof_manual_gating_frequency_clean.csv"), show_col_types = FALSE)
meta <- readr::read_csv(file.path(root, "data", "tables", "metadata.csv"), show_col_types = FALSE) |>
  dplyr::mutate(cytof_id = as.character(cytof_id))

cytof <- cytof |> dplyr::mutate(cytof_id = as.character(cytof_id))

harmonized <- cytof |>
  dplyr::transmute(
    sample_id = cytof_id,
    B_cells = B.cells,
    CD4T = CD4T,
    CD8T = CD8T,
    gdT = TCRgd.cells,
    Tregs = Tregs,
    NK_cells = NK,
    Neutrophils = Neutrophils,
    Eosinophils = Eosinophils,
    Basophils = Basophils,
    Monocytes = Classical.monocytes + Intermediate.monocytes + Nonclassical.monocytes,
    mDC = mDCs,
    pDC = pDCs,
    Platelets = NA_real_,
    MSC = NA_real_
  ) |>
  dplyr::inner_join(
    meta |> dplyr::select(cytof_id, subject_id, timepoint, group_feeding, group_delivery),
    by = c("sample_id" = "cytof_id")
  ) |>
  dplyr::mutate(
    baby_id = subject_id,
    baby_age = TIMEPOINT_DAYS[as.character(timepoint)],
    cohort = "Philippines (AMPLIO)",
    Location = "Philippines",
    group = group_feeding,
    mode_delivery = group_delivery
  ) |>
  dplyr::select(
    sample_id, baby_id, cohort, Location, baby_age,
    B_cells, CD4T, CD8T, gdT, Tregs, NK_cells, Neutrophils, Eosinophils,
    Basophils, Monocytes, mDC, pDC, Platelets, MSC,
    group, mode_delivery
  )

stopifnot(nrow(harmonized) == nrow(cytof))

out_path <- file.path(root, "data", "tables", "cross_cohort_lineage_harmonized.csv")
readr::write_csv(harmonized, out_path)

cat("Wrote", nrow(harmonized), "harmonized rows (", dplyr::n_distinct(harmonized$baby_id), "subjects) to", out_path, "\n")
