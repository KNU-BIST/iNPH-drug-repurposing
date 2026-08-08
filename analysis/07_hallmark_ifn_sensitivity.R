# =============================================================================
# Hallmark IFN-alpha pathway sensitivity analysis
# Script: 07_hallmark_ifn_sensitivity.R
# =============================================================================
#
# Project
# -------
# Systems biology prioritizes drug candidates for interferon-activated
# monocytes in idiopathic normal pressure hydrocephalus
#
# Author
# ------
# Md. Masudul Haque
#
# Purpose
# -------
# Test whether prioritized compound connectivity is reproduced using the
# non-overlapping component of the MSigDB Hallmark IFN-alpha-response gene
# set.
#
# Query definition
# ----------------
# The versioned 97-gene Hallmark set is compared with the complete 300-gene
# iNPH signature. Genes with direct overlap are excluded. The remaining
# genes are submitted as an UP-only LINCS2 query.
#
# Inputs
# ------
# data/processed/inph_ifn_signature.rds
# ExperimentHub resource EH7297
# data/raw/compoundinfo_beta.txt
#
# Outputs
# -------
# data/processed/hallmark_ifn_sensitivity_result.rds
# data/processed/hallmark_ifn_alpha_pathway_sensitivity_analysis.csv
#
# Reproducibility notes
# ---------------------
# - The 97-gene Hallmark IFN-alpha-response set is embedded below with its
#   Entrez identifiers.
# - Overlap is calculated against the complete finalized 300-gene iNPH
#   signature (150 UP + 150 DOWN), not against the UP set alone.
# - Overlapping Hallmark genes are excluded before the UP-only LINCS2 query.
# - This is a pathway-level sensitivity analysis, not a positive control,
#   disease-specific validation, or drug-class enrichment analysis.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 . Environment and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "signatureSearch",
  "ExperimentHub",
  "readr",
  "stringr",
  "purrr",
  "rhdf5",
  "BiocParallel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run renv::restore() before executing this script."
  )
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(signatureSearch)
  library(ExperimentHub)
  library(readr)
  library(stringr)
})

PROCESSED_DIR <- here::here("data", "processed")
RAW_DIR <- here::here("data", "raw")

SIGNATURE_PATH <- file.path(
  PROCESSED_DIR,
  "inph_ifn_signature.rds"
)

COMPOUNDINFO_BETA_PATH <- file.path(
  RAW_DIR,
  "compoundinfo_beta.txt"
)

required_inputs <- c(
  SIGNATURE_PATH,
  COMPOUNDINFO_BETA_PATH
)

input_info <- file.info(required_inputs)

bad_inputs <- required_inputs[
  !file.exists(required_inputs) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(bad_inputs) > 0L) {
  stop(
    "Required input file(s) missing or empty: ",
    paste(bad_inputs, collapse = ", ")
  )
}

cat(
  "Phase 1 signature MD5: ",
  unname(tools::md5sum(SIGNATURE_PATH)),
  "\n",
  sep = ""
)

cat(
  "compoundinfo_beta.txt MD5: ",
  unname(tools::md5sum(COMPOUNDINFO_BETA_PATH)),
  "\n",
  sep = ""
)

# Match the primary LINCS2 script's serial execution setting.
BiocParallel::register(
  BiocParallel::SerialParam()
)

# Resolve the ExperimentHub cache with a project-local fallback.
ext_root <- Sys.getenv("INPH_EXTERNAL_DATA")

eh_cache <- if (nzchar(ext_root)) {
  file.path(ext_root, "experiment_hub_cache")
} else {
  file.path(RAW_DIR, "experiment_hub_cache")
}

dir.create(
  eh_cache,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!dir.exists(eh_cache)) {
  warning(
    "Could not access ExperimentHub cache path '",
    eh_cache,
    "'. Falling back to the project-local cache."
  )

  eh_cache <- file.path(
    RAW_DIR,
    "experiment_hub_cache"
  )

  dir.create(
    eh_cache,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

if (!dir.exists(eh_cache)) {
  stop(
    "Unable to create or access the ExperimentHub cache directory: ",
    eh_cache
  )
}

ExperimentHub::setExperimentHubOption(
  "CACHE",
  eh_cache
)

cat(
  "ExperimentHub cache resolved to: ",
  eh_cache,
  "\n\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Section 1 . Versioned MSigDB Hallmark IFN-alpha-response set
# -----------------------------------------------------------------------------
hallmark_ifn_alpha <- tibble::tribble(
  ~symbol,     ~entrez,
  "ADAR",      "103",     "B2M",       "567",     "BATF2",     "116071",
  "BST2",      "684",     "C1S",       "716",     "CASP1",     "834",
  "CASP8",     "841",     "CCRL2",     "9034",    "CD47",      "961",
  "CD74",      "972",     "CMPK2",     "129607",  "CNP",       "1267",
  "CSF1",      "1435",    "CXCL10",    "3627",    "CXCL11",    "6373",
  "DDX60",     "55601",   "DHX58",     "79132",   "EIF2AK2",   "5610",
  "ELF1",      "1997",    "EPSTI1",    "94240",   "MVB12A",    "93343",
  "TENT5A",    "55603",   "CMTR1",     "23070",   "GBP2",      "2634",
  "GBP4",      "115361",  "GMPR",      "2766",    "HERC6",     "55008",
  "HLA-C",     "3107",    "IFI27",     "3429",    "IFI30",     "10437",
  "IFI35",     "3430",    "IFI44",     "10561",   "IFI44L",    "10964",
  "IFIH1",     "64135",   "IFIT2",     "3433",    "IFIT3",     "3437",
  "IFITM1",    "8519",    "IFITM2",    "10581",   "IFITM3",    "10410",
  "IL15",      "3600",    "IL4R",      "3566",    "IL7",       "3574",
  "IRF1",      "3659",    "IRF2",      "3660",    "IRF7",      "3665",
  "IRF9",      "10379",   "ISG15",     "9636",    "ISG20",     "3669",
  "LAMP3",     "27074",   "LAP3",      "51056",   "LGALS3BP",  "3959",
  "LPAR6",     "10161",   "LY6E",      "4061",    "MOV10",     "4343",
  "MX1",       "4599",    "NCOA7",     "135112",  "NMI",       "9111",
  "NUB1",      "51667",   "OAS1",      "4938",    "OASL",      "8638",
  "OGFR",      "11054",   "PARP12",    "64761",   "PARP14",    "54625",
  "PARP9",     "83666",   "PLSCR1",    "5359",    "PNPT1",     "87178",
  "HELZ2",     "85441",   "PROCR",     "10544",   "PSMA3",     "5684",
  "PSMB8",     "5696",    "PSMB9",     "5698",    "PSME1",     "5720",
  "PSME2",     "5721",    "RIPK2",     "8767",    "RNF31",     "55072",
  "RSAD2",     "91543",   "RTP4",      "64108",   "SAMD9",     "54809",
  "SAMD9L",    "219285",  "SELL",      "6402",    "SLC25A28",  "81894",
  "SP110",     "3431",    "STAT2",     "6773",    "TAP1",      "6890",
  "TDRD7",     "23424",   "TMEM140",   "55281",   "TRAFD1",    "10906",
  "TRIM14",    "9830",    "TRIM21",    "6737",    "TRIM25",    "7706",
  "TRIM26",    "7726",    "TRIM5",     "85363",   "TXNIP",     "10628",
  "UBA7",      "7318",    "UBE2L6",    "9246",    "USP18",     "11274",
  "WARS1",     "7453"
)
stopifnot(nrow(hallmark_ifn_alpha) == 97)
stopifnot(!anyDuplicated(hallmark_ifn_alpha$symbol))
stopifnot(!anyDuplicated(hallmark_ifn_alpha$entrez))

cat("Hallmark IFN-alpha set loaded: 97 genes (Liberzon et al. 2015, M5911)\n\n")

# -----------------------------------------------------------------------------
# Section 2 . Live overlap check against the finalized 300-gene iNPH signature
# -----------------------------------------------------------------------------

inph_signature <- readRDS(
  SIGNATURE_PATH
)

if (
  !is.list(inph_signature) ||
  !all(c("up", "down") %in% names(inph_signature))
) {
  stop(
    "inph_ifn_signature.rds must contain a list with 'up' and 'down' elements."
  )
}

if (
  length(inph_signature$up) != 150L ||
  length(inph_signature$down) != 150L ||
  anyDuplicated(inph_signature$up) ||
  anyDuplicated(inph_signature$down) ||
  length(intersect(inph_signature$up, inph_signature$down)) != 0L
) {
  stop(
    "Phase 1 signature failed the expected disjoint 150-UP/150-DOWN validation."
  )
}

overlap_up <- intersect(
  hallmark_ifn_alpha$symbol,
  inph_signature$up
)

overlap_down <- intersect(
  hallmark_ifn_alpha$symbol,
  inph_signature$down
)

overlap_all <- union(
  overlap_up,
  overlap_down
)

nonoverlap_genes <- setdiff(
  hallmark_ifn_alpha$symbol,
  overlap_all
)

cat("Hallmark overlap with the finalized iNPH signature:\n")
cat("  Overlap with UP set   :", length(overlap_up), "\n")
cat("  Overlap with DOWN set :", length(overlap_down), "\n")
cat("  Unique total overlap  :", length(overlap_all), "of 97\n")
cat("  Non-overlapping genes :", length(nonoverlap_genes), "of 97\n\n")

if (
  length(overlap_all) != 43L ||
  length(nonoverlap_genes) != 54L
) {
  stop(
    "Hallmark/iNPH overlap differs from the finalized analysis: expected ",
    "43 overlapping and 54 non-overlapping genes."
  )
}

query_genes <- hallmark_ifn_alpha |>
  dplyr::filter(
    symbol %in% nonoverlap_genes
  )

if (
  nrow(query_genes) != 54L ||
  anyDuplicated(query_genes$symbol) ||
  anyDuplicated(query_genes$entrez)
) {
  stop(
    "Hallmark sensitivity query should contain exactly 54 unique non-overlapping genes."
  )
}


# -----------------------------------------------------------------------------
# Section 3 . Access LINCS2 and run the UP-only query once
# -----------------------------------------------------------------------------

eh <- ExperimentHub::ExperimentHub(
  ask = FALSE
)

lincs2_path <- eh[["EH7297"]]

if (
  length(lincs2_path) != 1L ||
  !file.exists(lincs2_path) ||
  is.na(file.info(lincs2_path)$size) ||
  file.info(lincs2_path)$size <= 0
) {
  stop(
    "ExperimentHub resource EH7297 was not materialized as a valid local file."
  )
}

if (!rhdf5::H5Fis_hdf5(lincs2_path)) {
  stop(
    "ExperimentHub resource EH7297 is not a valid HDF5 file: ",
    lincs2_path
  )
}

cat(
  "LINCS2 database path: ",
  lincs2_path,
  "\n",
  sep = ""
)

lincs_genes <- rhdf5::h5read(
  lincs2_path,
  "rownames"
)[, 1]

n_query_in_lincs <- sum(
  query_genes$entrez %in% lincs_genes
)

cat(
  "Non-overlapping Hallmark query genes represented in LINCS2: ",
  n_query_in_lincs,
  " of ",
  nrow(query_genes),
  "\n",
  sep = ""
)

if (n_query_in_lincs != 49L) {
  stop(
    "LINCS2 coverage of the 54-gene non-overlapping Hallmark query differs ",
    "from the finalized analysis: expected 49 represented Entrez IDs."
  )
}

qsig_hallmark <- signatureSearch::qSig(
  query = list(
    upset = query_genes$entrez,
    downset = character(0)
  ),
  gess_method = "LINCS",
  refdb = lincs2_path
)

cat(
  "Running gess_lincs() - UP-only query, silent, 5-30 min...\n"
)

# tau = TRUE is retained to match the primary LINCS2 workflow. The sensitivity
# interpretation below uses NCS and WTCS FDR; Tau is not interpreted.
hallmark_result <- signatureSearch::gess_lincs(
  qsig_hallmark,
  sortby = "NCS",
  tau = TRUE,
  workers = 1
)

cat("Scoring complete.\n")

# For this UP-only query, signatureSearch returns a populated S4 result slot
# even where the result() generic may not dispatch cleanly on some versions.
if (!"result" %in% methods::slotNames(hallmark_result)) {
  stop(
    "Hallmark gess_lincs result does not contain the expected 'result' slot."
  )
}

hallmark_table <- methods::slot(
  hallmark_result,
  "result"
)

required_result_columns <- c(
  "pert",
  "type",
  "cell",
  "NCS",
  "WTCS_FDR"
)

missing_result_columns <- setdiff(
  required_result_columns,
  colnames(hallmark_table)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Hallmark LINCS2 result is missing required column(s): ",
    paste(missing_result_columns, collapse = ", ")
  )
}

if (nrow(hallmark_table) != 136460L) {
  stop(
    "Hallmark LINCS2 signature count mismatch: expected 136,460, observed ",
    nrow(hallmark_table),
    "."
  )
}

cat(
  "Total signatures scored: ",
  nrow(hallmark_table),
  "\n\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 4 . Compound identification - compoundinfo_beta.txt (broader
# coverage than GSE70138_Broad_LINCS_pert_info.txt.gz)
# -----------------------------------------------------------------------------
compoundinfo_beta <- readr::read_tsv(
  COMPOUNDINFO_BETA_PATH,
  show_col_types = FALSE
)

required_compound_columns <- c(
  "pert_id",
  "cmap_name"
)

missing_compound_columns <- setdiff(
  required_compound_columns,
  colnames(compoundinfo_beta)
)

if (length(missing_compound_columns) > 0L) {
  stop(
    "compoundinfo_beta.txt is missing required column(s): ",
    paste(missing_compound_columns, collapse = ", ")
  )
}

lookup_brd_ids <- function(drug_name) {
  compoundinfo_beta |>
    dplyr::filter(tolower(cmap_name) == tolower(drug_name)) |>
    dplyr::pull(pert_id) |>
    unique()
}

thp1_lookup <- function(brd_id) {
  hallmark_table |> dplyr::filter(pert == brd_id, toupper(cell) == "THP1")
}

# -----------------------------------------------------------------------------
# Section 5 . Final table
# ONE ROW PER COMPOUND. Where a compound has multiple LINCS2 entries (e.g.
# tofacitinib, methotrexate), the entry with real THP1 data is reported if
# one exists; otherwise the compound is reported once as not identified/
# not profiled - never fanned out into multiple rows.
# -----------------------------------------------------------------------------
comparator_drugs <- tibble::tribble(
  ~Compound,             ~Role,                                        ~selected_brd_id,
  "rosuvastatin",        "Primary lead (iNPH query)",                  "BRD-K82941592",
  "roflumilast",         "Primary lead (iNPH query)",                  "BRD-K03194791",
  "baricitinib",         "Secondary candidate (iNPH query, post hoc)", "BRD-K53581288",
  "tofacitinib",         "JAK inhibitor (approved)",                   "BRD-K36280065",
  "ruxolitinib",         "JAK inhibitor (approved)",                   "BRD-K53972329",
  "upadacitinib",        "JAK inhibitor (approved)",                   NA_character_,
  "filgotinib",          "JAK inhibitor (approved)",                   "BRD-K16803204",
  "methotrexate",        "Conventional DMARD",                         NA_character_,
  "sulfasalazine",       "Conventional DMARD",                         "BRD-K10670311",
  "leflunomide",         "Conventional DMARD",                         "BRD-K78692225",
  "hydroxychloroquine",  "Conventional DMARD",                         "BRD-A99117172"
)

NOT_IN_LIBRARY <- "Not in LINCS2 screening library (checked compoundinfo_beta.txt)"
NOT_PROFILED   <- "Profiled in LINCS2, no THP1 data"

build_row <- function(Compound, Role, selected_brd_id) {
  if (is.na(selected_brd_id) || selected_brd_id == "") {
    ids <- lookup_brd_ids(Compound)
    status <- if (length(ids) == 0L) NOT_IN_LIBRARY else NOT_PROFILED
    return(tibble::tibble(Compound = Compound, Role = Role,
                          THP1_Profiling_Status = status,
                          NCS = NA_real_, FDR = NA_real_,
                          Direction = NA_character_, Interpretation = NA_character_))
  }

  hit <- thp1_lookup(selected_brd_id)

  if (nrow(hit) == 0L) {
    return(tibble::tibble(Compound = Compound, Role = Role,
                          THP1_Profiling_Status = NOT_PROFILED,
                          NCS = NA_real_, FDR = NA_real_,
                          Direction = NA_character_, Interpretation = NA_character_))
  }

  if (nrow(hit) > 1L) {
    stop("Expected exactly one THP1 result for ", Compound, " (", selected_brd_id,
         "), found ", nrow(hit), " - check for duplicate entries.")
  }

  ncs_raw <- hit$NCS[[1]]
  fdr <- hit$WTCS_FDR[[1]]

  if (
    is.na(ncs_raw) ||
    is.na(fdr)
  ) {
    stop(
      "THP1 result for ",
      Compound,
      " (",
      selected_brd_id,
      ") contains missing NCS or WTCS FDR."
    )
  }

  direction <- if (ncs_raw < 0) {
    "Reversal (negative)"
  } else {
    "Activating (positive)"
  }

  interpretation <- dplyr::case_when(
    ncs_raw < -0.5 & fdr < 0.05 ~ "Significant reversal",
    ncs_raw < 0 ~ "Negative, not significant",
    ncs_raw >= 0 & fdr < 0.05 ~ "Significant positive connectivity",
    TRUE ~ "Positive connectivity, not significant"
  )

  tibble::tibble(
    Compound = Compound,
    Role = Role,
    THP1_Profiling_Status = "Profiled",
    NCS = round(ncs_raw, 4),
    FDR = fdr,
    Direction = direction,
    Interpretation = interpretation
  )
}

final_table <- purrr::pmap_dfr(
  comparator_drugs,
  build_row
) |>
  dplyr::mutate(
    Compound = stringr::str_to_sentence(Compound)
  )

if (
  nrow(final_table) != 11L ||
  anyDuplicated(final_table$Compound)
) {
  stop(
    "Hallmark sensitivity table must contain exactly one row for each of 11 compounds."
  )
}

assert_close <- function(actual, expected, tolerance, label) {
  if (
    length(actual) != 1L ||
    is.na(actual) ||
    abs(actual - expected) > tolerance
  ) {
    stop(
      label,
      " mismatch: expected ",
      expected,
      " +/- ",
      tolerance,
      ", observed ",
      paste(actual, collapse = ", ")
    )
  }
  invisible(TRUE)
}

get_compound_row <- function(compound_name) {
  row <- final_table |>
    dplyr::filter(
      Compound == compound_name
    )

  if (nrow(row) != 1L) {
    stop(
      "Expected exactly one Hallmark sensitivity row for ",
      compound_name,
      "."
    )
  }

  row
}

expected_profiled <- tibble::tribble(
  ~Compound,              ~NCS,     ~FDR,       ~Interpretation,
  "Rosuvastatin",         -2.2427,  0.0001043,  "Significant reversal",
  "Roflumilast",          -1.2500,  0.0001043,  "Significant reversal",
  "Baricitinib",          -0.7772,  0.08445,    "Negative, not significant",
  "Tofacitinib",          -1.0656,  0.0002255,  "Significant reversal",
  "Ruxolitinib",          -1.2454,  0.0001043,  "Significant reversal",
  "Leflunomide",           1.2074,  0.0001043,  "Significant positive connectivity",
  "Hydroxychloroquine",    0.9130,  0.002171,   "Significant positive connectivity"
)

for (i in seq_len(nrow(expected_profiled))) {
  expected_row <- expected_profiled[i, ]
  observed_row <- get_compound_row(
    expected_row$Compound
  )

  if (observed_row$THP1_Profiling_Status != "Profiled") {
    stop(
      expected_row$Compound,
      " should have a THP1 profile in the finalized Hallmark sensitivity analysis."
    )
  }

  assert_close(
    observed_row$NCS,
    expected_row$NCS,
    0.002,
    paste0(expected_row$Compound, " Hallmark THP1 NCS")
  )

  fdr_tolerance <- if (expected_row$FDR > 0.01) {
    0.001
  } else {
    0.00001
  }

  assert_close(
    observed_row$FDR,
    expected_row$FDR,
    fdr_tolerance,
    paste0(expected_row$Compound, " Hallmark THP1 FDR")
  )

  if (observed_row$Interpretation != expected_row$Interpretation) {
    stop(
      expected_row$Compound,
      " interpretation mismatch: expected '",
      expected_row$Interpretation,
      "', observed '",
      observed_row$Interpretation,
      "'."
    )
  }
}

expected_unprofiled_status <- c(
  "Upadacitinib" = NOT_IN_LIBRARY,
  "Filgotinib" = NOT_PROFILED,
  "Methotrexate" = NOT_PROFILED,
  "Sulfasalazine" = NOT_PROFILED
)

for (compound_name in names(expected_unprofiled_status)) {
  observed_row <- get_compound_row(
    compound_name
  )

  if (
    observed_row$THP1_Profiling_Status !=
      expected_unprofiled_status[[compound_name]]
  ) {
    stop(
      compound_name,
      " profiling-status mismatch."
    )
  }

  if (
    !is.na(observed_row$NCS) ||
    !is.na(observed_row$FDR)
  ) {
    stop(
      compound_name,
      " should not have a THP1 NCS/FDR value in the finalized table."
    )
  }
}

cat(
  "=== Supplementary Table 4: Hallmark IFN-alpha sensitivity analysis ===\n"
)
cat(
  "(orthogonal interferon-pathway sensitivity check; not a positive control or validation)\n\n"
)
print(
  final_table,
  n = Inf
)
cat("\n")

# -----------------------------------------------------------------------------
# Section 6 . Save and verify outputs
# -----------------------------------------------------------------------------

RESULT_OUTPUT_PATH <- file.path(
  PROCESSED_DIR,
  "hallmark_ifn_sensitivity_result.rds"
)

TABLE_OUTPUT_PATH <- file.path(
  PROCESSED_DIR,
  "hallmark_ifn_alpha_pathway_sensitivity_analysis.csv"
)

saveRDS(
  hallmark_result,
  RESULT_OUTPUT_PATH
)

readr::write_csv(
  final_table,
  TABLE_OUTPUT_PATH
)

output_paths <- c(
  RESULT_OUTPUT_PATH,
  TABLE_OUTPUT_PATH
)

output_info <- file.info(
  output_paths
)

bad_outputs <- output_paths[
  !file.exists(output_paths) |
    is.na(output_info$size) |
    output_info$size <= 0
]

if (length(bad_outputs) > 0L) {
  stop(
    "One or more Hallmark sensitivity outputs were not created correctly: ",
    paste(bad_outputs, collapse = ", ")
  )
}

cat("============================================================\n")
cat("  HALLMARK IFN-ALPHA SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(
  "  Hallmark set                         : ",
  nrow(hallmark_ifn_alpha),
  " genes\n",
  sep = ""
)
cat(
  "  Overlap with iNPH UP                : ",
  length(overlap_up),
  "\n",
  sep = ""
)
cat(
  "  Overlap with iNPH DOWN              : ",
  length(overlap_down),
  "\n",
  sep = ""
)
cat(
  "  Unique overlap excluded             : ",
  length(overlap_all),
  "\n",
  sep = ""
)
cat(
  "  Non-overlapping genes in UP query   : ",
  length(nonoverlap_genes),
  "\n",
  sep = ""
)
cat(
  "  Query Entrez IDs represented in LINCS2: ",
  n_query_in_lincs,
  "\n",
  sep = ""
)
cat(
  "  LINCS2 signatures scored            : ",
  nrow(hallmark_table),
  "\n",
  sep = ""
)
cat("============================================================\n")

cat(
  "\nCitation: Liberzon et al. 2015 (PMID 26771021), MSigDB ",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE (M5911).\n",
  sep = ""
)

cat(
  "Query definition: ",
  length(nonoverlap_genes),
  " of 97 Hallmark genes were non-overlapping with the complete 300-gene ",
  "iNPH signature; ",
  length(overlap_all),
  " overlapping genes were excluded (",
  length(overlap_up),
  " overlapping UP, ",
  length(overlap_down),
  " overlapping DOWN).\n",
  sep = ""
)

cat(
  "LINCS2 coverage: ",
  n_query_in_lincs,
  " of ",
  length(nonoverlap_genes),
  " non-overlapping Hallmark Entrez IDs were represented in the LINCS2 ",
  "reference; the query was UP-only.\n",
  sep = ""
)

cat(
  "Connectivity methodology: same NCS/WTCS framework and tau=TRUE execution ",
  "setting as 03_run_lincs2_connectivity.R; Tau is not interpreted.\n"
)

cat(
  "Interpretation: pathway-level sensitivity analysis; not a positive control, ",
  "disease-specific validation, or class-enrichment analysis.\n"
)

cat("\nOutputs saved and verified:\n")
for (path in output_paths) {
  cat(
    "  ",
    path,
    "\n",
    sep = ""
  )
}
