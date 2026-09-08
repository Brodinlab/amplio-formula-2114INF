# Builds the analysis-ready vaccine-specific IgG table from the raw Luminex
# bead-array data shared by the assay provider (never modified in place):
#
#   Raw source (OneDrive, not part of this git checkout):
#   "Brodinlab KI/company collabs/Nestlé/AMPLIO Formula RCT 2114INF/datasets/
#    vaccine responses/241217_shared_data_to_brodin_lab.zip"
#   -> Combined_data.csv (long format: Sample.ID, Assay, MFI, z_score, bead_count)
#
# Matches the manuscript's own Methods ("Vaccine-specific antibodies (Luminex)"):
# 36 vaccine-related antigens (IgG, Luminex Flexmap 3D), MFI z-score normalized
# per antigen by the assay provider (see their Assay Report 241217.docx --
# z_score = scale(MFI) within each antigen, already computed, used as-is here
# rather than re-normalizing).
#
# The raw panel has 44 antigen columns total: 36 real antigens + 8 isotype/
# bead controls (Bare bead, Goat IgG, mouse IgG, anti-human IgG, each run in
# both NHS/EDC and Mix&Go coupling batches) -- controls are dropped here.
# Per Petter's decision (2026-09-06): all 36 non-control antigens are
# retained, matching the manuscript's literal "36 vaccine-related antigens"
# count -- this includes antigens outside the infant vaccination schedule
# (CMV, EBV, RSV, HPV, S. pneumoniae) alongside the actual schedule vaccines
# (diphtheria, tetanus, pertussis, HBV, poliovirus, rotavirus) and
# maternally-transferred measles/mumps/rubella titers. NOTE: no Haemophilus
# influenzae type b (Hib) antigen is present in this panel despite being
# named in the manuscript's antigen-category list -- flagged as a
# manuscript-text/data discrepancy, not resolved here.
#
# Sample.ID (e.g. "PHL001-0001-BPV1") matches metadata$vaccine_id directly.
#
# Output: data/tables/vaccine_igg_zscore.csv (long format: subject_id,
#         timepoint, group_feeding, group_delivery, antigen, z_score,
#         MFI, bead_count)

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
})

load_required_packages(c("dplyr", "readr", "stringr"))

root <- get_repo_root()

zip_path <- "/Users/petterbrodin/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Brodinlab KI/company collabs/Nestlé/AMPLIO Formula RCT 2114INF/datasets/vaccine responses/241217_shared_data_to_brodin_lab.zip"
tmp_dir <- tempfile("vaccine_")
dir.create(tmp_dir)
utils::unzip(zip_path, files = "241217_shared_data_to_brodin_lab/Combined_data.csv", exdir = tmp_dir)

combined <- readr::read_csv(
  file.path(tmp_dir, "241217_shared_data_to_brodin_lab", "Combined_data.csv"),
  show_col_types = FALSE
) |>
  dplyr::rename(sample_id = `Sample.ID`, antigen = Assay)

CONTROL_ANTIGENS <- c(
  "Bare bead", "Goat IgG", "mouse IgG", "anti-human IgG",
  "Bare bead Mix&Go", "Goat IgG Mix&Go", "mouse IgG Mix&Go", "anti-human IgG Mix&Go"
)
combined <- combined |>
  dplyr::filter(!antigen %in% CONTROL_ANTIGENS, !sample_id %in% c("Pool", "Blank"))

metadata <- readr::read_csv(file.path(root, "data", "tables", "metadata.csv"), show_col_types = FALSE)

vaccine <- combined |>
  dplyr::inner_join(
    metadata |> dplyr::distinct(subject_id, timepoint, group_feeding, group_delivery, vaccine_id),
    by = c("sample_id" = "vaccine_id")
  ) |>
  dplyr::transmute(
    subject_id, timepoint, group_feeding, group_delivery,
    antigen = stringr::str_squish(antigen), # a few antigen names have trailing/embedded whitespace from the vendor sheet
    MFI, z_score, bead_count
  )

n_antigens <- dplyr::n_distinct(vaccine$antigen)
stopifnot(n_antigens == 36) # sanity check against the manuscript's stated panel size

readr::write_csv(vaccine, file.path(root, "data", "tables", "vaccine_igg_zscore.csv"))
cat("Wrote", nrow(vaccine), "rows,", n_antigens, "antigens,", dplyr::n_distinct(vaccine$subject_id), "subjects\n")
cat("Samples per timepoint:\n")
print(vaccine |> dplyr::distinct(subject_id, timepoint) |> dplyr::count(timepoint))

unlink(tmp_dir, recursive = TRUE)
