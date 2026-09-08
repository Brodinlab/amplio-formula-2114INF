# Table panel for Figure 5 ("Module 3 features"): fine-grained WGCNA
# Module 3 ("brown") feature list, split into cell populations and
# proteins, sorted by module membership (kME) -- per Petter's request,
# 2026-09-08, to add the module composition as a table alongside the
# module's other Figure 5 panels (wgcna_finegrained_brown_splitviolin.R,
# wgcna_module3_feature_heatmap.R).
#
# Layout matches Petter's final assembled Fig5: one title ("Module 3
# features"), then a single row of side-by-side tables -- cell populations
# (narrow, 3 rows) followed by the 37 proteins split into 3 columns of
# ~13/13/11 rows each (rather than one long 37-row list) so the panel reads
# left-to-right instead of requiring a long vertical scan.
#
# Reuses the already-generated tables (wgcna_module3_feature_tables.R) --
# no recomputation.
#
# Output:
#   output/figures/manuscript/Fig5_wgcna_module3_table.pdf

suppressPackageStartupMessages(source("scripts/lib/common.R"))
load_required_packages(c("dplyr", "readr", "gridExtra", "grid"))

root <- get_repo_root()

proteins <- readr::read_csv(file.path(root, "output", "tables", "wgcna_module3_proteins.csv"), show_col_types = FALSE) |>
  dplyr::mutate(kME = round(kME, 2)) |>
  dplyr::rename(Protein = protein)
cell_populations <- readr::read_csv(file.path(root, "output", "tables", "wgcna_module3_cell_populations.csv"), show_col_types = FALSE) |>
  dplyr::mutate(kME = round(kME, 2)) |>
  dplyr::rename(`Cell pop.` = population)

table_theme <- gridExtra::ttheme_minimal(
  core = list(fg_params = list(fontsize = 8, hjust = 0, x = 0.05)),
  colhead = list(fg_params = list(fontsize = 8, fontface = "bold", hjust = 0, x = 0.05)),
  padding = grid::unit(c(3, 2), "mm")
)

# Split the kME-sorted protein list into 3 side-by-side column blocks
# (ceiling(n/3) rows each, last block shorter) -- preserves sort order
# within and across columns (col 1 = highest kME, ... col 3 = lowest).
n_proteins <- nrow(proteins)
chunk_size <- ceiling(n_proteins / 3)
chunk_id <- ceiling(seq_len(n_proteins) / chunk_size)
protein_chunks <- split(proteins, chunk_id)

# grid.arrange centers grobs of unequal height within their layout cell --
# pad every column to the same row count with blank rows so headers land on
# the same line across all 4 tables (top-aligned) instead of the shorter
# cell-population table drifting down to the vertical middle.
target_rows <- max(nrow(cell_populations), vapply(protein_chunks, nrow, integer(1)))
pad_table <- function(df, n) {
  if (nrow(df) >= n) return(df)
  blank <- df[rep(1, n - nrow(df)), , drop = FALSE]
  blank[] <- ""
  dplyr::bind_rows(df |> dplyr::mutate(dplyr::across(dplyr::everything(), as.character)), blank)
}
cell_populations_padded <- pad_table(cell_populations, target_rows)
protein_chunks_padded <- lapply(protein_chunks, pad_table, n = target_rows)

cell_grob <- gridExtra::tableGrob(cell_populations_padded, rows = NULL, theme = table_theme)
protein_grobs <- lapply(protein_chunks_padded, gridExtra::tableGrob, rows = NULL, theme = table_theme)

title_grob <- grid::textGrob("Module 3 features", x = 0, hjust = 0, gp = grid::gpar(fontsize = 12, fontface = "bold"))
gap <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))

# grid.arrange's "grobheight" unit does not reliably resolve through its own
# nested gtable layout -- measure each grob's own height directly off its
# $heights slot (resolves correctly against an open device's font metrics)
# inside a throwaway device, then use those as fixed inch values.
grDevices::pdf(NULL)
row_h <- max(vapply(c(list(cell_grob), protein_grobs), function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE), numeric(1)))
grDevices::dev.off()

title_h <- 0.3 # in
gap_h <- 0.15 # in
pdf_height <- title_h + gap_h + row_h + 0.3 # + outer margin
pdf_width <- 9.5

pdf(file.path(root, "output", "figures", "manuscript", "Fig5_wgcna_module3_table.pdf"), width = pdf_width, height = pdf_height)
gridExtra::grid.arrange(
  title_grob, gap,
  gridExtra::arrangeGrob(
    cell_grob, protein_grobs[[1]], protein_grobs[[2]], protein_grobs[[3]],
    ncol = 4, widths = grid::unit(c(1.3, 2, 2, 2), "in")
  ),
  ncol = 1, heights = grid::unit(c(title_h, gap_h, row_h), "in")
)
dev.off()

cat("Saved: output/figures/manuscript/Fig5_wgcna_module3_table.pdf\n")
