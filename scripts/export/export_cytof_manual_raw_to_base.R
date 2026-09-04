# Builds data/tables/cytof_manual_gating_frequency.csv from Kanth's manually-
# gated FlowJo hierarchy export, via the flowjo-gating-freqs skill
# (github.com/Brodinlab/brodinlab-skills) flattened relative-frequency table.
#
# Source (external, not committed -- raw data must never be copied into the
# repo per lab convention): OneDrive
#   "Brodinlab KI/company collabs/Nestlé/AMPLIO Formula RCT 2114INF/datasets/Cytof/"
# The flattened intermediate (`flowjo_tool.py flatten ... --confirm-root
# 'CD45+'`) is saved alongside the raw FlowJo export in that same OneDrive
# folder as <raw-basename>_relative_frequencies.csv, not in this repo, per
# the same raw-data convention.
#
# Version history:
#   - v1 (preliminary): analysis_20260831.csv -> analysis_20260831_relative_frequencies.csv
#   - v1.1 (Kanth QC pass, 2026-09-03): analysis_20260831_v1.1.csv ->
#     analysis_20260831_v1.1_relative_frequencies.csv (same 251 samples, same
#     gating tree, revised cell counts only)
#   - v1.1 Kanth re-gate (2026-09-04): Analysis_v1.1_Amplio_Kanth.csv ->
#     Analysis_v1.1_Amplio_Kanth_relative_frequencies.csv. **Structural gating
#     tree change**, not just revised counts: the previous
#     Neu_Eosino -> {Neutrophils, Eosinophils} / Non-Neutrophils_Eosinophils
#     branch was replaced with Neutrophils / Non-neutrophils ->
#     {Eosinophils, Non-Eosinophils} directly under CD45+. This targets a
#     previously-flagged issue (see Fig5_manualgating_nlr.R) that neutrophil
#     fractions were biologically implausible (median ~12% of CD45+). The new
#     gating raises median Neutrophils to ~17% -- a real improvement but NOT
#     a full resolution (61% of samples still <20%, min still 0%); treat any
#     NLR-based finding as still provisional. `cytof_manual_gating_column_
#     mapping.csv` was updated accordingly: Neu_Eosino and Non-Neutrophils_
#     Eosinophils columns removed (no longer exist in the new tree),
#     Non-neutrophils/Non-Eosinophils added. Deep T-cell subset populations
#     (CD4/CD8 memory subsets, Tregs, TCRgd) are essentially unchanged
#     (mean |diff| <1 percentage point) -- only the myeloid/granulocyte/NK
#     side of the tree was re-gated.
#
# Output: data/tables/cytof_manual_gating_frequency.csv

suppressPackageStartupMessages({
  source("scripts/lib/common.R")
})

load_required_packages(c("dplyr", "readr"))

root <- get_repo_root()
tables_dir <- file.path(root, "data", "tables")

# Path to the flattened relative-frequency table -- update this to point at
# the current version's file each time this script is re-run for a new
# gating pass (see version history above).
flattened_path <- file.path(
  "/Users/petterbrodin/Library/CloudStorage/OneDrive-KarolinskaInstitutet",
  "Brodinlab KI/company collabs/Nestlé/AMPLIO Formula RCT 2114INF/datasets/Cytof",
  "Analysis_v1.1_Amplio_Kanth_relative_frequencies.csv"
)
stopifnot(
  "Flattened relative-frequency file not found -- run flowjo-gating-freqs `flatten` on the raw FlowJo export first" =
    file.exists(flattened_path)
)

mapping <- readr::read_csv(file.path(tables_dir, "cytof_manual_gating_column_mapping.csv"), show_col_types = FALSE)

flat <- readr::read_csv(flattened_path, show_col_types = FALSE)

missing_from_flat <- setdiff(mapping$original_flowjo_leaf_name, colnames(flat))
if (length(missing_from_flat) > 0) {
  stop(
    "Column mapping references population(s) not present in the flattened file: ",
    paste(missing_from_flat, collapse = ", "),
    " -- update cytof_manual_gating_column_mapping.csv for this gating version first."
  )
}

base <- flat |>
  dplyr::transmute(cytof_id = as.character(sample_id))

for (i in seq_len(nrow(mapping))) {
  base[[mapping$r_safe_column_name[i]]] <- flat[[mapping$original_flowjo_leaf_name[i]]]
}

readr::write_csv(base, file.path(tables_dir, "cytof_manual_gating_frequency.csv"))
cat(
  "Wrote", nrow(base), "sample(s) x", ncol(base) - 1, "population column(s) to",
  "data/tables/cytof_manual_gating_frequency.csv, from", basename(flattened_path), "\n"
)
