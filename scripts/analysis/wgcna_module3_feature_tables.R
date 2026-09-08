# Splits the fine-grained WGCNA Module 3 ("brown") feature list into two
# separate tables -- proteins (Olink) and cell populations (CyTOF) -- each
# sorted by module membership (kME) descending, per Petter's request,
# 2026-09-08. Reuses the already-computed assignment table
# (wgcna_finegrained_module_assignments.csv); no recomputation needed.
#
# Output:
#   output/tables/wgcna_module3_proteins.csv
#   output/tables/wgcna_module3_cell_populations.csv

suppressPackageStartupMessages(source("scripts/lib/common.R"))
load_required_packages(c("dplyr", "readr"))

root <- get_repo_root()

assignment <- readr::read_csv(
  file.path(root, "output", "tables", "wgcna_finegrained_module_assignments.csv"),
  show_col_types = FALSE
)
module_sizes <- assignment |> dplyr::filter(module != "grey") |> dplyr::count(module, sort = TRUE)
brown_number <- which(module_sizes$module == "brown")

module_features <- assignment |>
  dplyr::filter(module == "brown") |>
  dplyr::arrange(dplyr::desc(kME)) |>
  dplyr::select(feature_name, modality, kME)

proteins <- module_features |> dplyr::filter(modality == "Olink") |> dplyr::select(-modality) |> dplyr::rename(protein = feature_name)
cell_populations <- module_features |> dplyr::filter(modality == "CyTOF") |> dplyr::select(-modality) |> dplyr::rename(population = feature_name)

readr::write_csv(proteins, file.path(root, "output", "tables", "wgcna_module3_proteins.csv"))
readr::write_csv(cell_populations, file.path(root, "output", "tables", "wgcna_module3_cell_populations.csv"))

cat("Module", brown_number, "(\"brown\") --", nrow(proteins), "proteins,", nrow(cell_populations), "cell populations.\n\n")
cat("Proteins (by kME):\n")
print(proteins, n = Inf)
cat("\nCell populations (by kME):\n")
print(cell_populations, n = Inf)
cat("\nSaved:\n  output/tables/wgcna_module3_proteins.csv\n  output/tables/wgcna_module3_cell_populations.csv\n")
