load_amplio_data <- function(root) {
  load_required_packages(c("OlinkAnalyze", "dplyr"))

  # Metadata + CyTOF processed tables
  metadata <- utils::read.csv2(file.path(root, "data/processed/metadata_add_FCT.csv"), sep = ";")

  cytof_cluster_data <- utils::read.csv(
    file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_cluster_abundance_subsample.csv")
  )
  cytof_lineage_data <- utils::read.csv(
    file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_lineage_abundance_subsample.csv")
  )
  cytof_mfi_data <- utils::read.csv(
    file.path(root, "data/processed/cytof/AMPLIO_EXP-24-DV0748_cytof_flowsom_res_MFI_zscore.csv")
  )

  # Olink raw NPX
  olink_data <- OlinkAnalyze::read_NPX(
    file.path(
      root,
      "data/raw/olink/Amplio_XI-4096_Olink explore 384 INF-1/XI-4096_NPX_2025-02-27/XI-4096_NPX_2025-02-27.csv"
    )
  )

  # Olink QC + wide matrix + metadata join (mirrors `scripts/amplio_merge.Rmd`)
  load_required_packages(c("stringr", "tidyr", "tibble"))
  olink_data_sample <- olink_data |>
    dplyr::filter(Sample_Type == "SAMPLE") |>
    dplyr::filter(QC_Warning != "EXCLUDED") |>
    dplyr::filter(!stringr::str_detect(SampleID, "HD"))

  olink_data_wide <- olink_data_sample |>
    dplyr::select(SampleID, Assay, NPX) |>
    dplyr::mutate(Assay = make.names(Assay)) |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX) |>
    tibble::column_to_rownames("SampleID")

  olink_metadata <- data.frame(SampleID = rownames(olink_data_wide)) |>
    dplyr::left_join(
      metadata |>
        dplyr::mutate(SampleID = olink_id) |>
        dplyr::select(subject_id, SampleID, group_delivery, group_feeding, group_FCT, timepoint),
      by = "SampleID"
    )
  colnames(olink_metadata)[1] <- "olink_id"
  olink_metadata$group_feeding <- factor(olink_metadata$group_feeding, levels = c("CF", "EF"))
  olink_metadata$group_delivery <- factor(olink_metadata$group_delivery, levels = c("VAGINAL", "CAESAREAN"))
  olink_metadata$group_FCT <- factor(olink_metadata$group_FCT, levels = c("T1", "T2", "T3", "T4"))
  olink_metadata$timepoint <- factor(olink_metadata$timepoint, levels = c("V1", "V3", "V5"))

  # Vaccine
  vaccine_data <- utils::read.csv(
    file.path(root, "data/raw/vaccine/241217_shared_data_to_brodin_lab/Combined_data.csv"),
    row.names = 1
  )

  # Reformat CyTOF data as in `scripts/amplio_merge.Rmd`
  rownames(cytof_cluster_data) <- cytof_cluster_data$cytof_id
  cytof_cluster_data <- cytof_cluster_data[, setdiff(colnames(cytof_cluster_data), "cytof_id"), drop = FALSE]

  rownames(cytof_lineage_data) <- cytof_lineage_data$cytof_id
  cytof_lineage_data <- cytof_lineage_data[, setdiff(colnames(cytof_lineage_data), "cytof_id"), drop = FALSE]

  cytof_cluster_data[is.na(cytof_cluster_data)] <- 0
  cytof_lineage_data[is.na(cytof_lineage_data)] <- 0

  cytof_metadata <- metadata |>
    dplyr::select(
      subject_id, group_feeding, group_delivery, group_FCT, timepoint,
      sampling_date, cytof_id, cytof_plate, cytof_shipment
    ) |>
    dplyr::filter(!is.na(cytof_id)) |>
    dplyr::mutate(cytof_id = as.character(cytof_id))

  cytof_metadata <- cytof_metadata[match(rownames(cytof_cluster_data), cytof_metadata$cytof_id), ]
  cytof_metadata$group_feeding <- factor(cytof_metadata$group_feeding, levels = c("CF", "EF"))
  cytof_metadata$group_delivery <- factor(cytof_metadata$group_delivery, levels = c("VAGINAL", "CAESAREAN"))
  cytof_metadata$group_FCT <- factor(cytof_metadata$group_FCT, levels = c("T1", "T2", "T3", "T4"))
  cytof_metadata$timepoint <- factor(cytof_metadata$timepoint, levels = c("V1", "V3", "V5"))

  list(
    metadata = metadata,
    cytof_cluster_data = cytof_cluster_data,
    cytof_lineage_data = cytof_lineage_data,
    cytof_mfi_data = cytof_mfi_data,
    cytof_metadata = cytof_metadata,
    olink_data = olink_data,
    olink_data_wide = olink_data_wide,
    olink_metadata = olink_metadata,
    vaccine_data = vaccine_data
  )
}

