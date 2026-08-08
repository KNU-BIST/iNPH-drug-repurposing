# =============================================================================
# Phase 3 — LINCS2 transcriptional connectivity analysis
# Script: 03_run_lincs2_connectivity.R
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
# Query the 150-UP/150-DOWN iNPH IFN-activated monocyte signature against
# the LINCS2 L1000 reference database and prioritize compounds using
# THP1-specific normalized connectivity scores.
#
# Inputs
# ------
# data/processed/inph_ifn_signature.rds
# ExperimentHub resource EH7297
# data/raw/GSE70138_Broad_LINCS_pert_info.txt.gz
# data/raw/compoundinfo_beta.txt
#
# Outputs
# -------
# data/processed/lincs2_gess_result.rds
# data/processed/lincs2_compound_summary.rds
# data/processed/lincs2_significant_thp1_reversers.rds
# data/processed/lincs2_final_candidates.rds
# data/processed/lincs2_top200_thp1_reversal_candidates.csv
# data/processed/lincs2_significant_thp1_reversers.csv
#
# Ranking and filtering
# ---------------------
# 1. Rank compounds by full-precision THP1 NCS among the 597 compounds with
#    non-zero THP1 connectivity.
# 2. Retain significant strong reversers: NCS < -0.5 and WTCS FDR < 0.05.
# 3. Apply the curated regulatory-status review to the prespecified top 30
#    THP1 compounds plus post hoc baricitinib.
# =============================================================================

# -----------------------------------------------------------------------------
# Section 0 . Environment and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "ggplot2",
  "readr",
  "signatureSearch",
  "signatureSearchData",
  "ExperimentHub",
  "org.Hs.eg.db",
  "AnnotationDbi",
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
  library(ggplot2)
  library(signatureSearch)
  library(ExperimentHub)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

PROCESSED_DIR <- here::here("data", "processed")
RAW_DIR <- here::here("data", "raw")
DIAGNOSTIC_FIG_DIR <- here::here("figures", "diagnostic")

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIAGNOSTIC_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

SIGNATURE_PATH <- file.path(
  PROCESSED_DIR,
  "inph_ifn_signature.rds"
)

GSE_PERT_INFO_PATH <- file.path(
  RAW_DIR,
  "GSE70138_Broad_LINCS_pert_info.txt.gz"
)

COMPOUNDINFO_BETA_PATH <- file.path(
  RAW_DIR,
  "compoundinfo_beta.txt"
)

required_inputs <- c(
  SIGNATURE_PATH,
  GSE_PERT_INFO_PATH,
  COMPOUNDINFO_BETA_PATH
)

input_info <- file.info(required_inputs)

missing_or_empty <- required_inputs[
  !file.exists(required_inputs) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(missing_or_empty) > 0L) {
  stop(
    "Required input file(s) missing or empty: ",
    paste(missing_or_empty, collapse = ", ")
  )
}

cat(
  "Phase 1 signature MD5: ",
  unname(tools::md5sum(SIGNATURE_PATH)),
  "\n",
  sep = ""
)

cat(
  "GSE70138 metadata MD5: ",
  unname(tools::md5sum(GSE_PERT_INFO_PATH)),
  "\n",
  sep = ""
)

cat(
  "compoundinfo_beta.txt MD5: ",
  unname(tools::md5sum(COMPOUNDINFO_BETA_PATH)),
  "\n",
  sep = ""
)

# Use a serial BiocParallel backend for reproducible single-worker execution.
# gess_lincs() is also called with workers = 1 below. Some versions of
# signatureSearch may still emit an informational Windows backend warning;
# the manuscript uses NCS, WTCS, and WTCS FDR rather than Tau.
BiocParallel::register(
  BiocParallel::SerialParam()
)

# Resolve the ExperimentHub cache location with a project-local fallback.
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
    "Could not reach ExperimentHub cache path '",
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
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 1 . Load the Phase 1 signature
# -----------------------------------------------------------------------------
inph_signature <- readRDS(SIGNATURE_PATH)

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
    "Phase 1 signature failed the expected disjoint 150-UP/150-DOWN check."
  )
}

cat(
  "Phase 1 signature: UP ",
  length(inph_signature$up),
  " | DOWN ",
  length(inph_signature$down),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Section 2 . Access the LINCS2 database
# -----------------------------------------------------------------------------
cat("\nAccessing LINCS2 database via ExperimentHub...\n")

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

h5info <- rhdf5::h5ls(
  lincs2_path
)

print(
  h5info[
    h5info$otype == "H5I_DATASET",
    c("name", "dim")
  ]
)

lincs_genes <- rhdf5::h5read(
  lincs2_path,
  "rownames"
)[, 1]

cat(
  "LINCS2 genes: ",
  length(lincs_genes),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Section 3 . Gene ID conversion (HGNC symbol -> Entrez)
# -----------------------------------------------------------------------------
symbol_to_entrez <- function(symbols) {
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = symbols, column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  mapped <- mapped[!is.na(mapped)]
  cat("  Mapped:", length(mapped), "of", length(symbols), "genes\n")
  mapped
}
cat("Converting UP genes:\n");
up_entrez   <- symbol_to_entrez(inph_signature$up)
cat("Converting DOWN genes:\n");
down_entrez <- symbol_to_entrez(inph_signature$down)

n_up_mapped <- length(up_entrez)
n_down_mapped <- length(down_entrez)
n_up_lincs <- sum(up_entrez %in% lincs_genes)
n_down_lincs <- sum(down_entrez %in% lincs_genes)

cat("\nCoverage after conversion:\n")
cat(
  "  UP mapped      : ",
  n_up_mapped,
  " of 150\n",
  sep = ""
)
cat(
  "  DOWN mapped    : ",
  n_down_mapped,
  " of 150\n",
  sep = ""
)
cat(
  "  UP in LINCS2   : ",
  n_up_lincs,
  " of ",
  n_up_mapped,
  "\n",
  sep = ""
)
cat(
  "  DOWN in LINCS2 : ",
  n_down_lincs,
  " of ",
  n_down_mapped,
  "\n",
  sep = ""
)

if (
  n_up_mapped != 143L ||
  n_down_mapped != 137L ||
  n_up_lincs != 114L ||
  n_down_lincs != 101L
) {
  stop(
    "Gene-mapping coverage differs from the finalized analysis ",
    "(expected mapped UP/DOWN 143/137 and LINCS-covered UP/DOWN 114/101)."
  )
}

up_entrez_vec <- as.character(
  up_entrez
)

down_entrez_vec <- as.character(
  down_entrez
)

# -----------------------------------------------------------------------------
# Section 4 . Run LINCS WTCS connectivity scoring - ONCE
# -----------------------------------------------------------------------------
qsig_lincs2 <- qSig(
  query = list(upset = up_entrez_vec, downset = down_entrez_vec),
  gess_method = "LINCS", refdb = lincs2_path
)
cat("\nRunning gess_lincs() - silent, 5-30 min...\n")

# tau = TRUE is retained to reproduce the historical result object. The
# manuscript conclusions use NCS, WTCS, and WTCS FDR; Tau is not interpreted.
lincs2_result <- signatureSearch::gess_lincs(
  qsig_lincs2,
  sortby = "NCS",
  tau = TRUE,
  workers = 1
)

cat("Scoring complete.\n")

result_table <- signatureSearch::result(
  lincs2_result
)

required_result_columns <- c(
  "pert",
  "type",
  "cell",
  "NCS",
  "WTCS",
  "WTCS_FDR"
)

missing_result_columns <- setdiff(
  required_result_columns,
  colnames(result_table)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "LINCS2 result is missing required column(s): ",
    paste(missing_result_columns, collapse = ", ")
  )
}

if (nrow(result_table) != 136460L) {
  stop(
    "LINCS2 signature count mismatch: expected 136,460, observed ",
    nrow(result_table),
    "."
  )
}

cat(
  "Total signatures scored: ",
  nrow(result_table),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Section 5 . Compound-level accounting + ranking (BRD-only, THP1-primary)
# -----------------------------------------------------------------------------
cell_lines_present <- sort(
  unique(result_table$cell)
)

if (length(cell_lines_present) != 58L) {
  stop(
    "LINCS2 cell-line count mismatch: expected 58, observed ",
    length(cell_lines_present),
    "."
  )
}

cat(
  "\nCell lines represented: ",
  length(cell_lines_present),
  "\n",
  sep = ""
)

# Verify compound-level accounting.
all_compounds <- result_table |>
  dplyr::filter(type == "trt_cp") |>
  dplyr::distinct(pert)

brd_only <- all_compounds |>
  dplyr::filter(
    grepl("^BRD-", pert)
  )

n_non_brd <- nrow(all_compounds) -
  nrow(brd_only)

if (
  nrow(all_compounds) != 30456L ||
  n_non_brd != 45L ||
  nrow(brd_only) != 30411L
) {
  stop(
    "LINCS2 compound accounting differs from the finalized analysis: ",
    "expected 30,456 trt_cp entries, 45 non-BRD entries, and 30,411 ",
    "BRD compounds."
  )
}

cat(
  "\nTotal trt_cp compounds: ",
  nrow(all_compounds),
  "\n",
  sep = ""
)

cat(
  "Non-BRD entries excluded: ",
  n_non_brd,
  "\n",
  sep = ""
)

cat(
  "BRD-only compounds retained: ",
  nrow(brd_only),
  "\n",
  sep = ""
)

# Verify the finalized THP1 non-zero and strong-reversal counts.
thp1_hits <- result_table |>
  dplyr::filter(
    type == "trt_cp",
    toupper(cell) == "THP1",
    grepl("^BRD-", pert),
    !is.na(NCS),
    NCS != 0
  )

n_strong_thp1 <- sum(
  thp1_hits$NCS < -0.5
)

if (
  nrow(thp1_hits) != 597L ||
  n_strong_thp1 != 288L
) {
  stop(
    "THP1 accounting differs from the finalized analysis: expected ",
    "597 non-zero BRD compounds and 288 with NCS < -0.5."
  )
}

cat(
  "\nTHP1 non-zero BRD compounds: ",
  nrow(thp1_hits),
  "\n",
  sep = ""
)

cat(
  "THP1 NCS < -0.5: ",
  n_strong_thp1,
  " (",
  sprintf("%.1f", 100 * n_strong_thp1 / nrow(thp1_hits)),
  "%)\n",
  sep = ""
)

# --- Build the ranked compound table. BRD-only filter applied HERE, not just
#     reported afterward, so the 45 non-drug entries can never appear ranked. ---
thp1_rows <- result_table |>
  dplyr::filter(type == "trt_cp", toupper(cell) == "THP1", grepl("^BRD-", pert), !is.na(NCS))

duplicate_thp1 <- thp1_rows |>
  dplyr::count(pert, name = "n_thp1_rows") |>
  dplyr::filter(n_thp1_rows != 1L)

if (nrow(duplicate_thp1) > 0L) {
  stop("Expected one THP1 result per compound, but found duplicates for: ",
       paste(duplicate_thp1$pert, collapse = ", "))
}

compound_summary <- result_table |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert)) |>
  dplyr::group_by(pert) |>
  dplyr::summarise(
    n_cell_lines = dplyr::n_distinct(cell),
    thp1_NCS  = { value <- NCS[toupper(cell) == "THP1"];      if (length(value) == 1L) value else NA_real_ },
    thp1_WTCS = { value <- WTCS[toupper(cell) == "THP1"];     if (length(value) == 1L) value else NA_real_ },
    thp1_FDR  = { value <- WTCS_FDR[toupper(cell) == "THP1"]; if (length(value) == 1L) value else NA_real_ },
    median_NCS = median(NCS, na.rm = TRUE),
    min_NCS    = min(NCS, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  dplyr::arrange(thp1_NCS)

cat("\nTop 30 compounds ranked by THP1-specific NCS:\n")
print(compound_summary |> head(30), n = 30)


# -----------------------------------------------------------------------------
# Section 6 . Cell-line sensitivity - THP-1 and HL-60
# -----------------------------------------------------------------------------
monocyte_lines <- c("THP1", "HL60")
mono_results <- result_table |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert),
                toupper(cell) %in% monocyte_lines) |>
  dplyr::arrange(NCS)
cat("\nTHP-1 / HL-60 specific top hits:\n")
print(
  head(
    mono_results |>
      dplyr::select(
        pert,
        cell,
        WTCS,
        NCS,
        WTCS_FDR
      ),
    30
  )
)


# -----------------------------------------------------------------------------
# Section 7 . Compound identification and reviewed-candidate filtering
# -----------------------------------------------------------------------------
# Compound names are resolved from two LINCS metadata sources. Both tables can
# contain repeated pert_id values, so each is deduplicated before joining.

compoundinfo_gse <- readr::read_tsv(
  GSE_PERT_INFO_PATH,
  show_col_types = FALSE
)

compoundinfo_beta <- readr::read_tsv(
  COMPOUNDINFO_BETA_PATH,
  show_col_types = FALSE
)

required_gse_columns <- c(
  "pert_id",
  "pert_iname",
  "canonical_smiles",
  "inchi_key"
)

required_beta_columns <- c(
  "pert_id",
  "cmap_name"
)

missing_gse_columns <- setdiff(
  required_gse_columns,
  colnames(compoundinfo_gse)
)

missing_beta_columns <- setdiff(
  required_beta_columns,
  colnames(compoundinfo_beta)
)

if (length(missing_gse_columns) > 0L) {
  stop(
    "GSE70138 metadata is missing required column(s): ",
    paste(missing_gse_columns, collapse = ", ")
  )
}

if (length(missing_beta_columns) > 0L) {
  stop(
    "compoundinfo_beta.txt is missing required column(s): ",
    paste(missing_beta_columns, collapse = ", ")
  )
}

compoundinfo_gse_dedup <- compoundinfo_gse |>
  dplyr::select(
    pert_id,
    pert_iname,
    canonical_smiles,
    inchi_key
  ) |>
  dplyr::distinct(
    pert_id,
    .keep_all = TRUE
  )

compoundinfo_beta_dedup <- compoundinfo_beta |>
  dplyr::select(
    pert_id,
    cmap_name
  ) |>
  dplyr::distinct(
    pert_id,
    .keep_all = TRUE
  )

compound_named <- compound_summary |>
  dplyr::left_join(
    compoundinfo_gse_dedup,
    by = c("pert" = "pert_id")
  ) |>
  dplyr::left_join(
    compoundinfo_beta_dedup,
    by = c("pert" = "pert_id")
  ) |>
  dplyr::mutate(
    pert_iname = dplyr::case_when(
      !is.na(pert_iname) &
        pert_iname != "" &
        pert_iname != pert ~ pert_iname,
      !is.na(cmap_name) &
        cmap_name != "" &
        cmap_name != pert ~ cmap_name,
      TRUE ~ "UNIDENTIFIED"
    )
  ) |>
  dplyr::select(
    -cmap_name
  )

if (
  nrow(compound_named) != nrow(compound_summary) ||
  anyDuplicated(compound_named$pert)
) {
  stop(
    "Compound-name annotation changed the compound row count or introduced ",
    "duplicate BRD identifiers."
  )
}

# Define the 597-compound THP1 ranking used for manuscript ranks.
ranked_thp1 <- compound_named |>
  dplyr::filter(
    !is.na(thp1_NCS),
    thp1_NCS != 0
  ) |>
  dplyr::arrange(
    thp1_NCS
  ) |>
  dplyr::mutate(
    rank = dplyr::row_number()
  )

if (
  nrow(ranked_thp1) != 597L ||
  anyDuplicated(ranked_thp1$pert)
) {
  stop(
    "Expected exactly 597 uniquely ranked compounds with non-zero THP1 NCS."
  )
}

# Significant strong reversers are defined independently of regulatory status.
step2_significant <- compound_named |>
  dplyr::filter(
    !is.na(thp1_NCS),
    thp1_NCS < -0.5,
    !is.na(thp1_FDR),
    thp1_FDR < 0.05
  )

if (nrow(step2_significant) != 138L) {
  stop(
    "Significant strong-reverser count mismatch: expected 138, observed ",
    nrow(step2_significant),
    "."
  )
}

cat(
  "\nSignificant strong THP1 reversers ",
  "(NCS < -0.5 and FDR < 0.05): ",
  nrow(step2_significant),
  "\n",
  sep = ""
)

# Curated regulatory-status determinations are limited to the prespecified
# top-30 THP1 review plus post hoc baricitinib. Compounds outside that reviewed
# set remain explicitly unverified.
candidate_verification <- tibble::tribble(
  ~pert,           ~pert_iname_expected,      ~approval_status,                             ~verification_note,
  "BRD-K82941592", "rosuvastatin",            "Approved (FDA)",                             "FDA-approved statin; retained as a primary computational lead.",
  "BRD-K03194791", "roflumilast",             "Approved (FDA)",                             "FDA-approved PDE4 inhibitor; retained as a primary computational lead.",
  "BRD-K43002773", "GDC-0068",                "Investigational",                            "= ipatasertib, AKT inhibitor. Phase II-III oncology trials, no FDA/EMA approval. Verified 2026-07.",
  "BRD-K53581288", "baricitinib",             "Approved (FDA/EMA)",                         "Olumiant; FDA- and EMA-authorized JAK1/2 inhibitor. Retained as a post hoc mechanistically aligned secondary candidate.",
  "BRD-K92303087", "oglemilast",              "Investigational",                            "PDE4 inhibitor (same class as roflumilast). Phase II only (asthma/COPD), never approved. Verified 2026-07.",
  "BRD-K12683773", "methyl-fasudil",          "Investigational",                            "Fasudil analog (ROCK inhibitor). Parent fasudil approved only in Japan (1995), never FDA-approved; no evidence this analog was independently developed. Verified 2026-07.",
  "BRD-K40624912", "ZM-39923",                "Investigational",                            "JAK1/3 inhibitor. Confirmed research-use-only across multiple vendors; never entered clinical development. Verified 2026-07.",
  "BRD-K59962020", "CHEMBL-374350",           "Investigational",                            "Obscure preclinical screening compound; name is its own ChEMBL ID, no trade/common name. No clinical development evidence. Verified 2026-07.",
  "BRD-K68548958", "C-646",                   "Investigational",                            "p300/CBP HAT inhibitor. Confirmed research-only across multiple vendors. Verified 2026-07.",
  "BRD-K92317137", "UNIDENTIFIED",            "Excluded",                                   "No name resolved from the GSE70138 or compoundinfo_beta metadata; excluded from the reviewed approved-candidate set.",
  "BRD-K59019422", "UNIDENTIFIED",            "Excluded",                                   "No name resolved from the GSE70138 or compoundinfo_beta metadata; excluded from the reviewed approved-candidate set.",
  # Additional regulatory-status determinations from the version-controlled candidate-verification record.
  "BRD-A02367930", "ethinyl-estradiol",       "Approved (FDA; ingredient)",                 "Active ingredient in multiple FDA-approved combination contraceptives; LINCS entry is the ingredient, not a standalone product. Source: FDA label 021676s020.",
  "BRD-K93255255", "indoximod",               "Investigational",                            "IDO-pathway inhibitor, Phase 2; no FDA/EMA authorization. Source: NCT04049669.",
  "BRD-K75641298", "metoclopramide",          "Approved (FDA)",                             "FDA-labeled (Reglan and generics). Source: DailyMed.",
  "BRD-K59013864", "AM-24",                   "Investigational",                            "= Bobel-24/2,4,6-triiodophenol; Broad metadata lists Phase 2. No current FDA/EMA authorization.",
  "BRD-K84925803", "T-5224",                  "Investigational / discontinued",             "AP-1 inhibitor, Phase 2; development did not lead to approval. Source: J Med Chem 10.1021/jm5004733.",
  "BRD-K07859598", "GW-5074",                 "Preclinical / research compound",            "c-Raf inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K14993104", "bemegride",               "Historical/obsolete; status unclear",        "Historical CNS stimulant/barbiturate antidote; NCATS lists 1958 approval but no current FDA/EMA product located. Status genuinely unresolved - do not cite either way without further check.",
  "BRD-K01436366", "XMD-1150",                "Preclinical / research compound",            "L1000 research compound; no clinical marketing authorization identified.",
  "BRD-K96615647", "AV-608",                  "Investigational / discontinued",             "NK1 receptor antagonist, Phase 2; development discontinued, no approval.",
  "BRD-A13084692", "troglitazone",            "Withdrawn",                                  "Approved 1997, withdrawn from market 2000 due to hepatotoxicity. Source: PMID 20020270.",
  "BRD-K00610438", "altanserin",              "Investigational",                            "5-HT2A antagonist/PET ligand, Phase 2; no marketing authorization.",
  "BRD-K60038276", "irbesartan",              "Approved (FDA)",                             "FDA-labeled ARB (Avapro and generics). Source: DailyMed.",
  "BRD-A14262390", "guanadrel",               "Discontinued (US)",                          "Formerly marketed as Hylorel; brand and generics discontinued in the US.",
  "BRD-K77987382", "mebendazole",             "Approved (FDA)",                             "FDA-approved anthelmintic (EMVERM), initial approval 1974.",
  "BRD-K82181189", "ischemin",                "Preclinical / research compound",            "CBP bromodomain inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K94830329", "ataluren",                "Not currently authorized (EMA non-renewal)", "Translarna had conditional EU authorization; European Commission did not renew it (28 Mar 2025). No FDA approval. Source: EMA EPAR.",
  "BRD-K29950728", "clomifene",               "Approved (FDA)",                             "FDA-labeled clomiphene citrate for ovulatory dysfunction.",
  "BRD-K70642949", "GSK-2334470",             "Preclinical / research compound",            "Selective PDK1 inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K01674964", "histamine",               "Approved (EMA formulation only)",            "Histamine dihydrochloride is EMA-authorized as Ceplene (+IL-2, AML maintenance); the LINCS 'histamine' entity may not match this specific formulation. Treat as borderline, not a clean pass. Source: EMA EPAR.",
  "BRD-K99749624", "linifanib",               "Investigational / development unsuccessful", "= ABT-869; reached Phase 3 in HCC but failed primary endpoint. No FDA/EMA approval."
)

# The regulatory-status table must correspond exactly to the prespecified
# top 30 THP1 compounds plus baricitinib at rank 39.
review_ids <- unique(
  c(
    ranked_thp1$pert[seq_len(30L)],
    "BRD-K53581288"
  )
)

if (length(review_ids) != 31L) {
  stop(
    "Reviewed-candidate set should contain exactly 31 compounds."
  )
}

if (
  nrow(candidate_verification) != 31L ||
  !setequal(
    candidate_verification$pert,
    review_ids
  )
) {
  stop(
    "The curated candidate-verification table no longer matches the ",
    "prespecified top-30 review plus post hoc baricitinib."
  )
}

# Use the curated expected name only as a fallback for reviewed compounds that
# remain unresolved after the two public metadata joins.
name_fallback <- candidate_verification |>
  dplyr::select(
    pert,
    pert_iname_expected
  )

compound_named <- compound_named |>
  dplyr::left_join(
    name_fallback,
    by = "pert"
  ) |>
  dplyr::mutate(
    pert_iname = dplyr::case_when(
      pert_iname == "UNIDENTIFIED" &
        !is.na(pert_iname_expected) &
        pert_iname_expected != "UNIDENTIFIED" ~ pert_iname_expected,
      TRUE ~ pert_iname
    )
  ) |>
  dplyr::select(
    -pert_iname_expected
  )

# Rebuild the ranked table after final name resolution. Connectivity values and
# ranks are unchanged; only human-readable names can be updated.
ranked_thp1 <- compound_named |>
  dplyr::filter(
    !is.na(thp1_NCS),
    thp1_NCS != 0
  ) |>
  dplyr::arrange(
    thp1_NCS
  ) |>
  dplyr::mutate(
    rank = dplyr::row_number()
  )

step2_significant <- compound_named |>
  dplyr::filter(
    !is.na(thp1_NCS),
    thp1_NCS < -0.5,
    !is.na(thp1_FDR),
    thp1_FDR < 0.05
  ) |>
  dplyr::left_join(
    candidate_verification |>
      dplyr::select(
        pert,
        approval_status,
        verification_note
      ),
    by = "pert"
  ) |>
  dplyr::mutate(
    approval_status = dplyr::coalesce(
      approval_status,
      "UNVERIFIED - outside reviewed set"
    )
  )

if (nrow(step2_significant) != 138L) {
  stop(
    "Significant strong-reverser count changed after annotation."
  )
}

reviewed_significant <- step2_significant |>
  dplyr::filter(
    pert %in% review_ids
  ) |>
  dplyr::left_join(
    ranked_thp1 |>
      dplyr::select(
        pert,
        rank
      ),
    by = "pert"
  ) |>
  dplyr::arrange(
    rank
  )

if (nrow(reviewed_significant) != 31L) {
  stop(
    "Expected all 31 reviewed compounds to satisfy the significant ",
    "strong-reversal criteria in the finalized analysis."
  )
}

cat(
  "\nReviewed top-30 compounds plus post hoc baricitinib:\n"
)

print(
  reviewed_significant |>
    dplyr::select(
      rank,
      pert,
      pert_iname,
      thp1_NCS,
      thp1_FDR,
      approval_status
    ),
  n = 31
)

# FDA-approved candidates are selected only within the 31-compound reviewed set.
# EMA-only/formulation-mismatch cases such as histamine are not treated as clean
# FDA-approved passes.
final_candidates <- reviewed_significant |>
  dplyr::filter(
    grepl(
      "^Approved \\(FDA",
      approval_status
    )
  )

if (nrow(final_candidates) != 8L) {
  stop(
    "Reviewed FDA-approved candidate count mismatch: expected 8, observed ",
    nrow(final_candidates),
    "."
  )
}

cat(
  "\nVerified FDA-approved reversers from the prespecified top-30 review ",
  "plus post hoc baricitinib: ",
  nrow(final_candidates),
  "\n",
  sep = ""
)

print(
  final_candidates |>
    dplyr::select(
      rank,
      pert,
      pert_iname,
      thp1_NCS,
      thp1_FDR,
      approval_status
    )
)

borderline_candidates <- reviewed_significant |>
  dplyr::filter(
    approval_status == "Approved (EMA formulation only)"
  )

if (nrow(borderline_candidates) > 0L) {
  cat(
    "\nBorderline regulatory/formulation case(s), not included in ",
    "final_candidates:\n",
    sep = ""
  )

  print(
    borderline_candidates |>
      dplyr::select(
        rank,
        pert,
        pert_iname,
        thp1_NCS,
        thp1_FDR,
        approval_status
      )
  )
}

# Manuscript tiers: two primary leads, post hoc baricitinib, and other
# FDA-approved reversers.
manuscript_tier <- tibble::tribble(
  ~pert,           ~tier,
  "BRD-K82941592", "Primary lead",
  "BRD-K03194791", "Primary lead",
  "BRD-K53581288", "Post hoc secondary candidate",
  "BRD-A02367930", "Other FDA-approved reverser",
  "BRD-K75641298", "Other FDA-approved reverser",
  "BRD-K60038276", "Other FDA-approved reverser",
  "BRD-K77987382", "Other FDA-approved reverser",
  "BRD-K29950728", "Other FDA-approved reverser"
)

tier_order <- c(
  "Primary lead",
  "Post hoc secondary candidate",
  "Other FDA-approved reverser"
)

final_candidates <- final_candidates |>
  dplyr::left_join(
    manuscript_tier,
    by = "pert"
  ) |>
  dplyr::arrange(
    match(tier, tier_order),
    thp1_NCS
  )

if (
  anyNA(final_candidates$tier) ||
  !setequal(
    final_candidates$pert,
    manuscript_tier$pert
  )
) {
  stop(
    "Final-candidate membership does not match the manuscript tier table."
  )
}

cat(
  "\nFinal candidates with manuscript tier:\n"
)

print(
  final_candidates |>
    dplyr::select(
      rank,
      pert,
      pert_iname,
      thp1_NCS,
      thp1_FDR,
      tier
    )
)

# Critical manuscript rank checks.
critical_ranks <- ranked_thp1 |>
  dplyr::filter(
    pert %in% c(
      "BRD-K82941592",
      "BRD-K03194791",
      "BRD-K53581288"
    )
  ) |>
  dplyr::select(
    pert,
    rank,
    thp1_NCS,
    thp1_FDR
  )

expected_critical_ranks <- c(
  "BRD-K82941592" = 1L,
  "BRD-K03194791" = 14L,
  "BRD-K53581288" = 39L
)

observed_critical_ranks <- setNames(
  critical_ranks$rank,
  critical_ranks$pert
)

if (
  !all(
    observed_critical_ranks[
      names(expected_critical_ranks)
    ] == expected_critical_ranks
  )
) {
  stop(
    "Critical candidate ranks differ from the finalized analysis."
  )
}

cat(
  "\nCritical THP1 ranks:\n"
)

print(
  critical_ranks |>
    dplyr::arrange(
      rank
    )
)

n_unverified <- sum(
  step2_significant$approval_status ==
    "UNVERIFIED - outside reviewed set"
)

cat(
  "\nSignificant strong reversers outside the regulatory-review set: ",
  n_unverified,
  " of ",
  nrow(step2_significant),
  ". No regulatory-status claim is made for these compounds.\n",
  sep = ""
)

# Descriptive check for known JAK/IFN-pathway compounds among significant hits.
jak_inhibitors <- c(
  "tofacitinib",
  "ruxolitinib",
  "baricitinib",
  "upadacitinib",
  "filgotinib",
  "abrocitinib",
  "fedratinib",
  "pacritinib",
  "momelotinib",
  "itacitinib"
)

ifn_modulators <- c(
  "anifrolumab",
  "sifalimumab",
  "rontalizumab",
  "hydroxychloroquine",
  "chloroquine",
  "quinacrine"
)

priority_hits <- step2_significant |>
  dplyr::filter(
    tolower(pert_iname) %in%
      tolower(
        c(
          jak_inhibitors,
          ifn_modulators
        )
      )
  )

cat(
  "\nKnown JAK/IFN-pathway compounds among significant hits:\n"
)

if (nrow(priority_hits) > 0L) {
  print(
    priority_hits
  )
} else {
  cat(
    "  (None identified by compound name among significant hits.)\n"
  )
}


# -----------------------------------------------------------------------------
# Section 8 . Diagnostic visualization of THP1 connectivity results
# -----------------------------------------------------------------------------
p1 <- compound_summary |>
  dplyr::filter(!is.na(thp1_NCS)) |>
  ggplot(aes(x = thp1_NCS)) +
  geom_histogram(bins = 60, fill = "#5e9ea0", color = NA, alpha = 0.8) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "#d4624a", linewidth = 0.7) +
  labs(x = "THP1-specific NCS", y = "Number of compounds",
       title = "LINCS2 THP1 connectivity score distribution",
       subtitle = "Query: iNPH IFN-activated monocyte signature | Negative NCS = reversal") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(
  file.path(DIAGNOSTIC_FIG_DIR, "lincs2_thp1_ncs_distribution.pdf"),
  p1,
  width = 7,
  height = 4.5
)

top30 <- ranked_thp1 |>
  dplyr::slice_head(n = 30) |>
  dplyr::mutate(
    pert = factor(pert, levels = rev(pert)),
    is_strong = thp1_NCS < -0.5
  )

p2 <- ggplot(top30, aes(x = thp1_NCS, y = pert)) +
  geom_segment(aes(x = 0, xend = thp1_NCS, y = pert, yend = pert), color = "grey80", linewidth = 0.5) +
  geom_point(aes(size = n_cell_lines, color = is_strong), alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = "#d4624a", "FALSE" = "#5e9ea0"),
                     labels = c("TRUE" = "THP1 NCS < -0.5", "FALSE" = "THP1 NCS \u2265 -0.5"), name = "") +
  scale_size_continuous(name = "Cell lines", range = c(2, 6)) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "#d4624a", linewidth = 0.5, alpha = 0.7) +
  labs(x = "THP1-specific NCS", y = NULL,
       title = "Top 30 reversal candidates - ranked by THP1 NCS",
       subtitle = "LINCS2 WTCS scoring") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right",
        panel.grid.major.y = element_blank())
ggsave(
  file.path(DIAGNOSTIC_FIG_DIR, "lincs2_top30_thp1_reversal_candidates.pdf"),
  p2,
  width = 8,
  height = 7
)


# -----------------------------------------------------------------------------
# Section 9 . Save outputs
# -----------------------------------------------------------------------------

output_paths <- c(
  lincs2_gess_result =
    file.path(PROCESSED_DIR, "lincs2_gess_result.rds"),
  lincs2_compound_summary =
    file.path(PROCESSED_DIR, "lincs2_compound_summary.rds"),
  lincs2_significant_thp1_reversers =
    file.path(PROCESSED_DIR, "lincs2_significant_thp1_reversers.rds"),
  lincs2_final_candidates =
    file.path(PROCESSED_DIR, "lincs2_final_candidates.rds"),
  lincs2_top200_thp1_reversal_candidates =
    file.path(PROCESSED_DIR, "lincs2_top200_thp1_reversal_candidates.csv"),
  lincs2_significant_thp1_reversers_csv =
    file.path(PROCESSED_DIR, "lincs2_significant_thp1_reversers.csv")
)

saveRDS(
  lincs2_result,
  output_paths[["lincs2_gess_result"]]
)

saveRDS(
  compound_named,
  output_paths[["lincs2_compound_summary"]]
)

saveRDS(
  step2_significant,
  output_paths[["lincs2_significant_thp1_reversers"]]
)

saveRDS(
  final_candidates,
  output_paths[["lincs2_final_candidates"]]
)

readr::write_csv(
  ranked_thp1 |>
    dplyr::slice_head(n = 200),
  output_paths[["lincs2_top200_thp1_reversal_candidates"]]
)

readr::write_csv(
  step2_significant,
  output_paths[["lincs2_significant_thp1_reversers_csv"]]
)

output_info <- file.info(
  unname(output_paths)
)

bad_outputs <- names(output_paths)[
  !file.exists(unname(output_paths)) |
    is.na(output_info$size) |
    output_info$size <= 0
]

if (length(bad_outputs) > 0L) {
  stop(
    "One or more Phase 3 outputs were not created correctly: ",
    paste(bad_outputs, collapse = ", ")
  )
}

legacy_outputs <- c(
  file.path(PROCESSED_DIR, "lincs2_step2_significant.rds"),
  file.path(PROCESSED_DIR, "phase3_top200_reversal_candidates.csv"),
  file.path(PROCESSED_DIR, "phase3_step2_significant_candidates.csv")
)

legacy_present <- legacy_outputs[
  file.exists(legacy_outputs)
]

if (length(legacy_present) > 0L) {
  cat(
    "\nNOTE: legacy Phase 3 output file(s) are still present from earlier runs:\n"
  )
  cat(
    paste0(
      "  - ",
      legacy_present,
      collapse = "\n"
    ),
    "\n"
  )
  cat(
    "They are not written by this revised script. Remove them before the final ",
    "fresh-clone reproducibility test to avoid accidental use by older scripts.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("  PHASE 3 COMPLETE\n")
cat("============================================================\n")
cat("  Total LINCS2 signatures       :", nrow(result_table), "\n")
cat("  Cell lines                    :", length(cell_lines_present), "\n")
cat("  BRD-only compounds            :", nrow(brd_only), "\n")
cat("  Non-zero THP1 compounds       :", nrow(ranked_thp1), "\n")
cat("  THP1 NCS < -0.5              :", n_strong_thp1, "\n")
cat("  Significant strong reversers  :", nrow(step2_significant), "\n")
cat("  Reviewed compounds            :", length(review_ids), "\n")
cat("  Verified FDA-approved reversers:", nrow(final_candidates), "\n")
cat("  Unreviewed significant hits   :", n_unverified, "\n")
cat("============================================================\n")

cat("\nOutputs saved and verified:\n")
for (path in unname(output_paths)) {
  cat("  ", path, "\n", sep = "")
}

cat("\nDiagnostic figures:\n")
cat(
  "  ",
  file.path(DIAGNOSTIC_FIG_DIR, "lincs2_thp1_ncs_distribution.pdf"),
  "\n",
  sep = ""
)
cat(
  "  ",
  file.path(DIAGNOSTIC_FIG_DIR, "lincs2_top30_thp1_reversal_candidates.pdf"),
  "\n",
  sep = ""
)
