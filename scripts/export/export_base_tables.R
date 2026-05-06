suppressPackageStartupMessages({
  source("code/lib/common.R")
})

load_required_packages(c("dplyr", "readr", "tidyr", "tibble", "stringr", "OlinkAnalyze"))

root <- get_repo_root()
out_dir <- file.path(root, "data", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1) Metadata (analysis-ready) ----
meta_path <- file.path(root, "data/processed/metadata_add_FCT.csv")
metadata <- readr::read_delim(meta_path, delim = ";", show_col_types = FALSE)
readr::write_csv(metadata, file.path(out_dir, "metadata.csv"))

## ---- 2) CyTOF (analysis-ready FlowSOM-derived frequencies) ----
cytof_cluster <- readr::read_csv(
  file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_cluster_abundance_subsample.csv"),
  show_col_types = FALSE
)
cytof_lineage <- readr::read_csv(
  file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_lineage_abundance_subsample.csv"),
  show_col_types = FALSE
)
readr::write_csv(cytof_cluster, file.path(out_dir, "cytof_cluster_frequency.csv"))
readr::write_csv(cytof_lineage, file.path(out_dir, "cytof_lineage_frequency.csv"))

## Optional: marker MFIs (already processed)
cytof_mfi_z <- readr::read_csv(
  file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_res_MFI_zscore.csv"),
  show_col_types = FALSE
)
readr::write_csv(cytof_mfi_z, file.path(out_dir, "cytof_cluster_mfi_zscore.csv"))

## ---- 3) Olink NPX (analysis-ready wide table + sample metadata) ----
olink_raw <- OlinkAnalyze::read_NPX(
  file.path(
    root,
    "data/raw/olink/Amplio_XI-4096_Olink explore 384 INF-1/XI-4096_NPX_2025-02-27/XI-4096_NPX_2025-02-27.csv"
  )
)

olink_sample <- olink_raw |>
  dplyr::filter(Sample_Type == "SAMPLE") |>
  dplyr::filter(QC_Warning != "EXCLUDED") |>
  dplyr::filter(!stringr::str_detect(SampleID, "HD"))

olink_npx_wide <- olink_sample |>
  dplyr::select(SampleID, Assay, NPX) |>
  dplyr::mutate(Assay = make.names(Assay)) |>
  tidyr::pivot_wider(names_from = Assay, values_from = NPX)

olink_metadata <- metadata |>
  dplyr::mutate(SampleID = olink_id) |>
  dplyr::select(subject_id, SampleID, group_delivery, group_feeding, group_FCT, timepoint)

readr::write_csv(olink_npx_wide, file.path(out_dir, "olink_npx_wide.csv"))
readr::write_csv(olink_metadata, file.path(out_dir, "olink_metadata.csv"))

message("Exported base tables to: ", out_dir)

