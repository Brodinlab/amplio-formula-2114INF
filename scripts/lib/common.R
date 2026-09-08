# Common utilities shared by figure scripts

# Visit-code -> nominal age mapping for the manually-gated CyTOF sub-study
# (V1 = enrollment/baseline, 14-35 days; V3 = 2 months; V5 = 4 months -- see
# the manuscript Methods). Used for x-axis positioning (TIMEPOINT_DAYS) and,
# per Petter's request 2026-09-03, for display labels in every figure instead
# of the internal V1/V3/V5 codes (TIMEPOINT_LABELS).
TIMEPOINT_DAYS <- c(V1 = 0, V3 = 60, V5 = 120)
TIMEPOINT_LABELS <- c(V1 = "Baseline", V3 = "2 months", V5 = "4 months")

# Relabel a vector of V1/V3/V5 codes to the display labels above, as an
# ordered factor (order follows TIMEPOINT_LABELS, i.e. chronological).
relabel_timepoint <- function(x) {
  factor(TIMEPOINT_LABELS[as.character(x)], levels = TIMEPOINT_LABELS)
}

# Manual-gating CyTOF outlier exclusion, final set (Petter's decision,
# 2026-09-04): samples with Euclidean distance from centroid > 10 in
# plate-corrected CLR (Aitchison) space, computed on Kanth's re-gated v1.1
# base table -- see Fig5_manualgating_outlier_check.R for the detection
# method and Fig5_manualgating_pca_biplot_dist10_excluded.R for the
# threshold choice. Supersedes the earlier 2-sample list found ad hoc on the
# pre-v1.1 gating (453612193, 453610960 only -- 453611359 is new to v1.1).
# Single source of truth: every script that excludes outliers from the
# manual-gating table should reference this constant, not a local copy.
MANUAL_GATING_OUTLIER_IDS <- c("453612193", "453611359", "453610960")

# Manual-gating populations excluded from analysis (Petter's decision,
# 2026-09-04): intermediate/complement nodes in the gating hierarchy, not
# informative end-populations in their own right --
#   NonNK               = CD45+ minus NK (parent split, not a population)
#   NonTB               = CD45+ minus T/B cells (parent split)
#   CD14neg.CD16neg     = monocyte-gate leftover after classical/intermediate/
#                         nonclassical are removed, not a defined cell type
#   NonEosinophils      = parent split under CD45+ (sibling of NonNeutrophils
#                         below, same rationale)
#   NonNeutrophils      = parent split under CD45+ (sibling of NonEosinophils
#                         above) -- added per Petter's confirmation, 2026-09-06
# Single source of truth: every script analysing the manual-gating table's
# populations should drop these via this constant, not a local copy. Applies
# to the population set, orthogonal to MANUAL_GATING_OUTLIER_IDS (samples).
MANUAL_GATING_EXCLUDED_POPULATIONS <- c("NonNK", "NonTB", "CD14neg.CD16neg", "NonEosinophils", "NonNeutrophils")

# Olink NPX outlier exclusion (Petter's decision, 2026-09-06): none. Unlike
# the manually-gated CyTOF check (3 samples clearly separated from the rest,
# distance>10 cliff), the Olink per-protein-z-scored distance-from-centroid
# check (Fig5_olink_outlier_check.R) showed a smooth decline with no clear
# separation -- the 8 samples crossing the median+3*MAD line are only mildly
# offset from the bulk, not a distinct outlier cluster. All samples are kept.
# Single source of truth: every script analysing Olink NPX data should
# reference this constant (empty for now), not assume no exclusion locally --
# so a future revision only needs to change it here.
OLINK_OUTLIER_IDS <- character(0)

# Olink NPX cytokine classification (Petter's request, 2026-09-06): the Olink
# Explore Inflammation panel bundles many non-cytokine analytes (receptors,
# enzymes, structural/ECM proteins, growth factors) alongside true cytokines
# -- e.g. the top age-associated proteins overall include LAMA4, COLEC12,
# CKAP4, AGRN, ARNT, HSD11B1, none of which are cytokines. No official Olink
# panel protein-class annotation file was found in this repo or on OneDrive
# (data/raw/olink is not present in this checkout), so this is a manually
# curated, immunology-standard cytokine list -- NOT sourced from an
# Olink-provided annotation. Categories included:
#   - Interleukins (IL-prefixed ligands; IL-receptor subunits like IL10RA,
#     IL2RB, IL4R, etc. are EXCLUDED as receptors, not cytokines themselves;
#     IL1RN/IL-1Ra is kept as it is itself a secreted regulatory cytokine)
#   - Interferons (IFNG; IFNGR1/IFNLR1 receptors excluded)
#   - TNF-superfamily LIGANDS (TNF, TNFSF*, FASLG; TNFRSF* receptors and
#     CD40/LTBR excluded as receptors)
#   - Chemokines (CCL*/CXCL*)
#   - Colony-stimulating factors (CSF1, CSF3)
#   - A small set of other classic secreted immune cytokines (EPO, OSM,
#     CRLF1, FLT3LG)
# Deliberately EXCLUDED as a distinct category, not cytokines by the
# conventional definition, even though Olink's Inflammation panel includes
# them: general growth factors (FGF*, VEGF*, PDGF*, EGF, TGFA, TGFB1, GDNF-
# family ligands ARTN/NRTN/PSPN, PGF, angiopoietins). This is a judgment call
# at the margins -- flag to Petter if any inclusion/exclusion looks wrong.
# Single source of truth: every script classifying Olink proteins as
# cytokines should reference this constant, not a local copy.
OLINK_CYTOKINE_PROTEINS <- c(
  # Interleukins
  "IL2", "IL4", "IL1B", "IL33", "IL13", "IL10", "IL24", "IL20", "IL11", "IL17F",
  "IL5", "IL17A", "IL17C", "IL17D", "IL1A", "IL7", "IL15", "IL6", "IL32", "IL16",
  "IL18", "IL12B", "IL1RN",
  # Interferons
  "IFNG",
  # TNF-superfamily ligands
  "TNF", "TNFSF11", "TNFSF10", "TNFSF12", "TNFSF13", "FASLG",
  # Chemokines
  "CXCL14", "CXCL12", "CXCL9", "CXCL10", "CXCL8", "CXCL6", "CXCL17", "CXCL3", "CXCL1",
  "CCL7", "CCL26", "CCL28", "CCL3", "CCL13", "CCL11", "CCL25", "CCL21", "CCL23",
  "CCL4", "CCL17", "CCL22", "CCL24", "CCL20",
  # Colony-stimulating factors
  "CSF3", "CSF1",
  # Other classic secreted immune cytokines
  "EPO", "OSM", "CRLF1", "FLT3LG"
)

# Vaccine Luminex IgG panel: antigen category classification (Petter's
# request, 2026-09-06), distinguishing which of the 36 antigens (see
# export_vaccine.R) reflect the infant's OWN active vaccine response versus
# maternally-transferred antibody, since these answer different biological
# questions and are easily conflated in a flat antigen list --
#   INFANT_SCHEDULE_VACCINE_ANTIGENS: pathogens covered by the Philippines
#     EPI schedule doses an infant would have received by 4 months
#     (diphtheria, tetanus, pertussis, hepatitis B, poliovirus, rotavirus --
#     typically given at 6/10/14 weeks) -- IgG here can plausibly reflect
#     the infant's own vaccine response.
#   MATERNAL_MMR_ANTIGENS: measles/mumps/rubella antigens -- MMR is not
#     given until ~9-12 months in this schedule, so any IgG detected at V1/
#     V3/V5 (0-4 months) is maternally-transferred (transplacental) antibody
#     and its postnatal decay, NOT an infant response.
#   NON_VACCINE_PATHOGEN_ANTIGENS: CMV, EBV, RSV, HPV, S. pneumoniae --
#     present in this Luminex panel but not part of the infant's
#     immunization schedule at all; reflect natural exposure (mostly
#     maternal/environmental) rather than vaccination of either kind.
# All three sets are disjoint and together cover all 36 non-control antigens
# exactly (stopifnot-checked in export_vaccine.R against the manuscript's
# stated panel size).
INFANT_SCHEDULE_VACCINE_ANTIGENS <- c(
  "Diphtheria mutated toxin", "Diphtheria Toxoid",
  "Clostridium tetani Tetanus Toxoid", "Tetanus Toxoid, Recombinant Heavy Chain Fragment C",
  "B. pertussis toxin (mutant)", "Bordetella pertussis Filamentous Hemagglutinin (FHA)",
  "Bordetella pertussis Filamentous Hemagglutinin (FHA) - Bulk antigen",
  "Bordetella pertussis Filamentous Hemagglutinin (FHA) - Nativeantigen",
  "B. pertussis Pertactin Protein[His]",
  "B. Pertussis whole-cell (strain tahoma I)", "B. Pertussis whole-cell (strain tahoma I) Mix&Go",
  "HBV Surface Antigen (subtype adw)",
  "Recombinant Poliovirus type 1 Capsid protein (strain Sabin)",
  "Recombinant Poliovirus type 2 VP3-VP1 capsid Protein [His]",
  "Recombinant Poliovirus type 3 VP3-VP1 capsid protein",
  "Rotavirus VP7 Protein", "RotaVirus (Strain SA-11)"
)
MATERNAL_MMR_ANTIGENS <- c(
  "Measles Virus Nucleoprotein (HEK293)", "Native Measles virus",
  "Mumps Virus Nucleoprotein Recombinant", "Mumps virus nucleoprotein, inactivated pathogen.",
  "Native Mumps virus", "Mumps virus nucleoprotein",
  "Rubella E1", "Rubella virus E1, C-terminal SHFc-tag", "Rubella Spike Ectodomain (E1-E2)",
  "Rubella virus nucleoprotein, C-terminal His-tag", "Rubella Virus Grade 4, natural antigen."
)
NON_VACCINE_PATHOGEN_ANTIGENS <- c(
  "Cytomegalovirus glycoprotein B (gB)", "Respiratory Syncytial Virus A Glycoprotein G",
  "HPV type 16 L1 Protein (full length)", "HPV type 18 L1 Protein (full length)",
  "Recombinant Human Papilloma Virus type 33 L1 protein (VLP)", "Recombinant HPV type 6 L1 protein (VLP)",
  "Epstein Barr Virus gp125", "S. pneumoniae Cell Wall Polysaccharide Antigen"
)

classify_vaccine_antigen <- function(antigen) {
  dplyr::case_when(
    antigen %in% INFANT_SCHEDULE_VACCINE_ANTIGENS ~ "Infant-schedule vaccine",
    antigen %in% MATERNAL_MMR_ANTIGENS ~ "Maternal/MMR-family",
    antigen %in% NON_VACCINE_PATHOGEN_ANTIGENS ~ "Non-vaccine pathogen",
    TRUE ~ NA_character_
  )
}

get_repo_root <- function() {
  # Prefer running from repo root; fall back to script location if possible.
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  # If user runs from code/ or code/figures/, walk up.
  if (basename(wd) %in% c("figures", "lib")) wd <- dirname(wd)
  if (basename(wd) == "code") wd <- dirname(wd)

  wd
}

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nInstall them (e.g. install.packages(...)) and re-run."
    )
  }
}

save_pdf <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "pdf"
  )
}

add_significance <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  "ns"
}

cohens_d_with_ci <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) {
    return(list(d = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, d_str = "N/A"))
  }
  mean1 <- mean(x, na.rm = TRUE)
  mean2 <- mean(y, na.rm = TRUE)
  var1 <- stats::var(x, na.rm = TRUE)
  var2 <- stats::var(y, na.rm = TRUE)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
  if (pooled_sd == 0 || is.na(pooled_sd)) {
    return(list(d = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, d_str = "N/A"))
  }
  d <- (mean1 - mean2) / pooled_sd
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2)))
  ci_lower <- d - 1.96 * se_d
  ci_upper <- d + 1.96 * se_d
  d_str <- paste0(round(d, 2), " [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]")
  list(d = d, ci_lower = ci_lower, ci_upper = ci_upper, d_str = d_str)
}

# Mood's median test: tests whether two samples come from populations with
# the same median (distinct from Wilcoxon/Mann-Whitney, which tests for a
# general stochastic/rank shift and only reduces to a median comparison
# under a location-shift assumption). Builds the standard 2x2
# above/at-or-below-grand-median x group contingency table; uses Fisher's
# exact test instead of chi-squared whenever an expected cell count < 5.
moods_median_test <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  grand_median <- stats::median(c(x, y))
  tab <- matrix(
    c(sum(x > grand_median), sum(x <= grand_median), sum(y > grand_median), sum(y <= grand_median)),
    nrow = 2
  )
  expected_ok <- all(suppressWarnings(stats::chisq.test(tab)$expected) >= 5)
  test <- if (expected_ok) stats::chisq.test(tab, correct = TRUE) else stats::fisher.test(tab)
  list(p_value = test$p.value, method = if (expected_ok) "chisq" else "fisher", grand_median = grand_median)
}

clr_transform <- function(x) {
  # Add small pseudocount to avoid zeros
  x <- x + 0.001
  geom_mean <- exp(rowMeans(log(x)))
  log(x / geom_mean)
}

# Per-population plate correction on the log scale (residuals + global mean,
# then exponentiate): positive by construction regardless of correction
# magnitude, unlike an additive correction on the raw percentage scale (see
# export_cytof_manual_clean.R header for the negative-percentage bug this
# fixed, 2026-09-03). Shared here so any script needing this exact
# correction (e.g. reconstructing a differently-population-filtered version
# of the "clean" table) doesn't reimplement it and risk drift.
# df must contain pop_cols and have plate_vec aligned row-for-row.
plate_correct_log_scale <- function(df, pop_cols, plate_vec, pseudocount = 0.001) {
  for (pop in pop_cols) {
    log_val <- log(df[[pop]] + pseudocount)
    fit <- stats::lm(log_val ~ factor(plate_vec))
    log_corrected <- stats::residuals(fit) + mean(log_val, na.rm = TRUE)
    df[[pop]] <- pmax(exp(log_corrected) - pseudocount, 0)
  }
  df
}

