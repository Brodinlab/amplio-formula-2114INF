load_base_tables <- function(root) {
  load_required_packages(c("readr", "dplyr"))

  base_dir <- file.path(root, "data", "tables")

  metadata <- readr::read_csv(file.path(base_dir, "metadata.csv"), show_col_types = FALSE)
  cytof_cluster <- readr::read_csv(file.path(base_dir, "cytof_cluster_frequency.csv"), show_col_types = FALSE)
  cytof_lineage <- readr::read_csv(file.path(base_dir, "cytof_lineage_frequency.csv"), show_col_types = FALSE)
  cytof_manual <- readr::read_csv(file.path(base_dir, "cytof_manual_gating_frequency.csv"), show_col_types = FALSE)
  olink_npx_wide <- readr::read_csv(file.path(base_dir, "olink_npx_wide.csv"), show_col_types = FALSE)
  olink_metadata <- readr::read_csv(file.path(base_dir, "olink_metadata.csv"), show_col_types = FALSE)

  # Normalize ID types for reliable joins
  if ("cytof_id" %in% colnames(metadata)) metadata$cytof_id <- as.character(metadata$cytof_id)
  if ("cytof_id" %in% colnames(cytof_cluster)) cytof_cluster$cytof_id <- as.character(cytof_cluster$cytof_id)
  if ("cytof_id" %in% colnames(cytof_lineage)) cytof_lineage$cytof_id <- as.character(cytof_lineage$cytof_id)
  if ("cytof_id" %in% colnames(cytof_manual)) cytof_manual$cytof_id <- as.character(cytof_manual$cytof_id)

  if ("SampleID" %in% colnames(olink_npx_wide)) olink_npx_wide$SampleID <- as.character(olink_npx_wide$SampleID)
  if ("SampleID" %in% colnames(olink_metadata)) olink_metadata$SampleID <- as.character(olink_metadata$SampleID)

  list(
    metadata = metadata,
    cytof_cluster = cytof_cluster,
    cytof_lineage = cytof_lineage,
    cytof_manual = cytof_manual,
    olink_npx_wide = olink_npx_wide,
    olink_metadata = olink_metadata
  )
}

