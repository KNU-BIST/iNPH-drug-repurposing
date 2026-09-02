# =============================================================================
# Phase 0 — Reconstruct canonical iNPH classical-monocyte analysis input
# Script: 00_prepare_canonical_monocyte_input.R
# =============================================================================
#
# Project
# -------
# Systems biology prioritizes drug candidates for interferon-activated
# monocytes in idiopathic normal pressure hydrocephalus
#
# Purpose
# -------
# Reconstruct the locked 22,197-cell classical-monocyte analysis population
# used in the manuscript from public GSE292141 PBMC 10x count matrices plus a
# version-controlled canonical cell-membership manifest.
#
# This replaces the former requirement for a pre-built binary
# data/processed/duy_monocytes_classical.rds input.
#
# Reproducibility boundary
# ------------------------
# The manifest preserves the exact historical cell membership and order. This
# script does not claim to independently rediscover the historical monocyte
# cluster by re-running the original stochastic clustering. Instead, it
# reconstructs the locked analysis population from public raw counts and
# verifies the analysis-critical identities and metadata before writing the
# canonical processed object.
#
# Historical validation before workflow adoption
# ----------------------------------------------
# The reconstruction was verified against the previously validated legacy
# canonical object and reproduced:
#   - 36,601 genes and 22,197 cells
#   - exact gene identifiers/order and cell identifiers/order
#   - exact patient/tissue/MoCA metadata
#   - exact raw count matrix (0 non-zero differences)
#   - identical Script 01 signature and DE results
#   - identical subcluster assignments, normalized data, PCA, Harmony, and UMAP
#
# Inputs
# ------
# data/canonical/duy_classical_monocyte_cell_manifest.tsv
# data/canonical/duy_classical_monocyte_manifest_provenance.txt
# data/raw/duy/extracted/<GSE292141 10x triplets>
# R/data_acquisition_and_loading.R
# R/pbmc_clustering_and_annotation.R
#
# Output
# ------
# data/processed/duy_monocytes_classical.rds
#
# Historical primary QC
# ---------------------
# nFeature_RNA >= 200
# percent_mt <= 20
# No total-count MAD filter is applied.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment and locked constants
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "Matrix",
  "Seurat",
  "DropletUtils",
  "dplyr",
  "tibble",
  "purrr",
  "logger"
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
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

RAW_DIR <- here::here("data", "raw")
PROCESSED_DIR <- here::here("data", "processed")
CANONICAL_DIR <- here::here("data", "canonical")

DUY_EXTRACTED <- file.path(RAW_DIR, "duy", "extracted")
MANIFEST_PATH <- file.path(
  CANONICAL_DIR,
  "duy_classical_monocyte_cell_manifest.tsv"
)
PROVENANCE_PATH <- file.path(
  CANONICAL_DIR,
  "duy_classical_monocyte_manifest_provenance.txt"
)
OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "duy_monocytes_classical.rds"
)
DATA_HELPER <- here::here("R", "data_acquisition_and_loading.R")
PBMC_HELPER <- here::here("R", "pbmc_clustering_and_annotation.R")

EXPECTED_MANIFEST_MD5 <- "28417a6aa7efd1dd8fdd6f67c09e7be4"
EXPECTED_LEGACY_RDS_MD5 <- "5ab87562a5b244dfee63ca85c9402ac3"
EXPECTED_N_RAW_PBMC <- 124047L
EXPECTED_N_QC_PBMC <- 122724L
EXPECTED_N_GENES <- 36601L
EXPECTED_N_CANONICAL <- 22197L
EXPECTED_N_PATIENTS <- 9L

EXPECTED_PATIENT_COUNTS <- c(
  "1" = 2852L,
  "2" = 1804L,
  "3" = 3797L,
  "4" = 1953L,
  "5" = 2261L,
  "7" = 1416L,
  "8" = 2230L,
  "9" = 1343L,
  "10" = 4541L
)


# -----------------------------------------------------------------------------
# Section 1 · Required-input checks
# -----------------------------------------------------------------------------

required_files <- c(
  MANIFEST_PATH,
  PROVENANCE_PATH,
  DATA_HELPER,
  PBMC_HELPER
)

required_info <- file.info(required_files)

bad_required <- required_files[
  !file.exists(required_files) |
    is.na(required_info$size) |
    required_info$size <= 0
]

if (length(bad_required) > 0L) {
  stop(
    "Required reconstruction file(s) missing or empty:\n  ",
    paste(bad_required, collapse = "\n  ")
  )
}

if (!dir.exists(DUY_EXTRACTED)) {
  stop("Duy extracted-data directory is missing: ", DUY_EXTRACTED)
}

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

if (
  file.exists(OUTPUT_RDS) &&
  Sys.getenv("INPH_OVERWRITE_CANONICAL_RDS", unset = "0") != "1"
) {
  stop(
    "Canonical output already exists: ", OUTPUT_RDS, "\n",
    "Refusing to overwrite it automatically. A strict fresh pipeline should ",
    "begin without this file. If replacement is intentional, set ",
    "INPH_OVERWRITE_CANONICAL_RDS=1 for that run."
  )
}

manifest_md5 <- unname(tools::md5sum(MANIFEST_PATH))

if (!identical(manifest_md5, EXPECTED_MANIFEST_MD5)) {
  stop(
    "Canonical manifest MD5 mismatch.\n",
    "Expected: ", EXPECTED_MANIFEST_MD5, "\n",
    "Observed: ", manifest_md5, "\n",
    "Do not reconstruct from an unverified manifest."
  )
}

cat("Canonical manifest MD5: ", manifest_md5, "\n", sep = "")

source(DATA_HELPER)
source(PBMC_HELPER)


# -----------------------------------------------------------------------------
# Section 2 · Manifest validation
# -----------------------------------------------------------------------------

manifest <- utils::read.delim(
  MANIFEST_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_manifest_cols <- c(
  "canonical_order",
  "cell_id",
  "patient",
  "tissue",
  "moca_group"
)

missing_manifest_cols <- setdiff(
  required_manifest_cols,
  colnames(manifest)
)

if (length(missing_manifest_cols) > 0L) {
  stop(
    "Canonical manifest is missing required column(s): ",
    paste(missing_manifest_cols, collapse = ", ")
  )
}

if (nrow(manifest) != EXPECTED_N_CANONICAL) {
  stop(
    "Canonical manifest cell count mismatch: expected ",
    EXPECTED_N_CANONICAL,
    ", observed ",
    nrow(manifest),
    "."
  )
}

if (!identical(
  as.integer(manifest$canonical_order),
  seq_len(EXPECTED_N_CANONICAL)
)) {
  stop("canonical_order is not exactly 1:22197 in locked historical order.")
}

if (anyNA(manifest$cell_id) || any(manifest$cell_id == "")) {
  stop("Canonical manifest contains missing or empty cell identifiers.")
}

if (anyDuplicated(manifest$cell_id)) {
  stop("Canonical manifest contains duplicated cell identifiers.")
}

if (length(unique(as.character(manifest$patient))) != EXPECTED_N_PATIENTS) {
  stop("Canonical manifest does not contain the expected 9 patients.")
}

if (!all(manifest$tissue == "PBMC")) {
  stop("Canonical manifest contains non-PBMC entries.")
}

if (!all(manifest$moca_group %in% c("High", "Low"))) {
  stop("Canonical manifest contains unexpected MoCA-group values.")
}

manifest_patient_counts <- table(as.character(manifest$patient))

if (!identical(
  as.integer(manifest_patient_counts[names(EXPECTED_PATIENT_COUNTS)]),
  as.integer(EXPECTED_PATIENT_COUNTS)
)) {
  stop("Canonical manifest per-patient cell counts do not match the locked population.")
}

cat("Canonical manifest cells: ", nrow(manifest), "\n", sep = "")


# -----------------------------------------------------------------------------
# Section 3 · Reconstruct analyzed PBMCs from public GSE292141 matrices
# -----------------------------------------------------------------------------

cat("\n--- Reconstructing analyzed GSE292141 PBMC cohort ---\n")

pbmc_raw <- load_duy_samples(
  extracted_dir = DUY_EXTRACTED,
  tissue_filter = "PBMC",
  drop_unknown_moca = TRUE
)

if (!inherits(pbmc_raw, "SingleCellExperiment")) {
  stop("load_duy_samples() did not return a SingleCellExperiment.")
}

if (ncol(pbmc_raw) != EXPECTED_N_RAW_PBMC) {
  stop(
    "Pre-QC PBMC count mismatch: expected 124,047; observed ",
    ncol(pbmc_raw),
    "."
  )
}

if (nrow(pbmc_raw) != EXPECTED_N_GENES) {
  stop(
    "Pre-QC gene count mismatch: expected 36,601; observed ",
    nrow(pbmc_raw),
    "."
  )
}

if (anyDuplicated(colnames(pbmc_raw))) {
  stop("Raw PBMC reconstruction contains duplicated cell identifiers.")
}

pbmc_raw <- fix_gene_symbols(pbmc_raw)

if (nrow(pbmc_raw) != EXPECTED_N_GENES || anyDuplicated(rownames(pbmc_raw))) {
  stop("Gene-symbol harmonization produced an unexpected feature structure.")
}

cat(
  "Pre-QC PBMC: ",
  nrow(pbmc_raw),
  " genes x ",
  ncol(pbmc_raw),
  " cells\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 4 · Apply historical primary QC only
# -----------------------------------------------------------------------------
# Exact finalized criteria:
#   nFeature_RNA >= 200
#   percent_mt <= 20
# No total-count MAD filter is applied.
#
# Seurat is used for the QC definitions so this matches the validated PBMC
# helper exactly. PCA/Harmony/clustering/UMAP are intentionally not run here;
# cell membership is locked by the canonical manifest.

cat("\n--- Applying historical primary QC ---\n")

qc_seu <- Seurat::CreateSeuratObject(
  counts = SingleCellExperiment::counts(pbmc_raw),
  meta.data = as.data.frame(SummarizedExperiment::colData(pbmc_raw))
)

qc_seu[["percent_mt"]] <- Seurat::PercentageFeatureSet(
  qc_seu,
  pattern = "^MT-"
)

qc_seu <- subset(
  qc_seu,
  subset = nFeature_RNA >= 200 & percent_mt <= 20
)

if (ncol(qc_seu) != EXPECTED_N_QC_PBMC) {
  stop(
    "Historical-QC PBMC count mismatch: expected 122,724; observed ",
    ncol(qc_seu),
    "."
  )
}

qc_ids <- colnames(qc_seu)
pbmc_qc <- pbmc_raw[, qc_ids, drop = FALSE]

pbmc_qc$nCount_RNA <- as.numeric(qc_seu$nCount_RNA)
pbmc_qc$nFeature_RNA <- as.numeric(qc_seu$nFeature_RNA)
pbmc_qc$percent_mt <- as.numeric(qc_seu$percent_mt)

cat("Post-QC PBMC cells: ", ncol(pbmc_qc), "\n", sep = "")


# -----------------------------------------------------------------------------
# Section 5 · Verify all locked IDs and critical metadata
# -----------------------------------------------------------------------------

missing_canonical <- setdiff(manifest$cell_id, colnames(pbmc_qc))

if (length(missing_canonical) > 0L) {
  missing_manifest <- manifest[
    manifest$cell_id %in% missing_canonical,
    ,
    drop = FALSE
  ]

  cat("\nMissing canonical cells by patient:\n")
  print(table(missing_manifest$patient))

  stop(
    length(missing_canonical),
    " of 22,197 canonical cell IDs are absent after historical QC."
  )
}

manifest_index <- match(manifest$cell_id, colnames(pbmc_qc))

if (anyNA(manifest_index)) {
  stop("Internal canonical-ID matching error.")
}

patient_match <- identical(
  as.character(pbmc_qc$patient[manifest_index]),
  as.character(manifest$patient)
)

tissue_match <- identical(
  as.character(pbmc_qc$tissue[manifest_index]),
  manifest$tissue
)

moca_match <- identical(
  as.character(pbmc_qc$moca_group[manifest_index]),
  manifest$moca_group
)

if (!patient_match) {
  stop("GEO-derived patient metadata do not match the canonical manifest.")
}

if (!tissue_match) {
  stop("GEO-derived tissue metadata do not match the canonical manifest.")
}

if (!moca_match) {
  stop("GEO-derived MoCA metadata do not match the canonical manifest.")
}


# -----------------------------------------------------------------------------
# Section 6 · Build the canonical SCE in locked historical cell order
# -----------------------------------------------------------------------------

canonical <- pbmc_qc[, manifest$cell_id, drop = FALSE]
canonical$cell_type <- rep("monocyte", ncol(canonical))

S4Vectors::metadata(canonical)$canonical_reconstruction <- list(
  manifest_file = "data/canonical/duy_classical_monocyte_cell_manifest.tsv",
  manifest_md5 = manifest_md5,
  historical_legacy_rds_md5 = EXPECTED_LEGACY_RDS_MD5,
  raw_source = "GSE292141",
  analyzed_pbmc_pre_qc = EXPECTED_N_RAW_PBMC,
  analyzed_pbmc_post_qc = EXPECTED_N_QC_PBMC,
  canonical_cells = EXPECTED_N_CANONICAL,
  qc = c("nFeature_RNA >= 200", "percent_mt <= 20"),
  reproducibility_boundary = paste(
    "Cell membership is locked by the historical canonical manifest;",
    "the original stochastic monocyte-cluster selection is not",
    "independently rediscovered by this reconstruction."
  )
)

if (nrow(canonical) != EXPECTED_N_GENES ||
    ncol(canonical) != EXPECTED_N_CANONICAL) {
  stop("Canonical reconstruction dimensions are incorrect.")
}

if (!identical(colnames(canonical), manifest$cell_id)) {
  stop("Canonical reconstruction does not preserve manifest cell order.")
}

if (!"counts" %in% SummarizedExperiment::assayNames(canonical)) {
  stop("Canonical reconstruction is missing the raw counts assay.")
}

if (length(unique(as.character(canonical$patient))) != EXPECTED_N_PATIENTS) {
  stop("Canonical reconstruction does not contain the expected 9 patients.")
}

if (!all(as.character(canonical$tissue) == "PBMC")) {
  stop("Canonical reconstruction contains non-PBMC cells.")
}

canonical_patient_counts <- table(as.character(canonical$patient))

if (!identical(
  as.integer(canonical_patient_counts[names(EXPECTED_PATIENT_COUNTS)]),
  as.integer(EXPECTED_PATIENT_COUNTS)
)) {
  stop("Canonical reconstruction per-patient counts do not match the locked population.")
}

# Normalize loader-generated sample metadata.
# The raw-data loader stores its temporary extraction directory in `Sample`
# and `metadata()$Samples`. Replace those session-specific paths with stable
# GEO accessions so repeated reconstructions produce deterministic metadata.
stable_sample_ids <- unique(as.character(canonical$sample_id))

canonical$Sample <- as.character(canonical$sample_id)

canonical_metadata <- S4Vectors::metadata(canonical)
canonical_metadata$Samples <- stable_sample_ids
S4Vectors::metadata(canonical) <- canonical_metadata

if (!identical(
  as.character(canonical$Sample),
  as.character(canonical$sample_id)
)) {
  stop("Failed to normalize the canonical Sample metadata.")
}

if (any(grepl("[/\\\\]", canonical$Sample))) {
  stop("Canonical Sample metadata still contains filesystem paths.")
}

if (!identical(
  S4Vectors::metadata(canonical)$Samples,
  stable_sample_ids
)) {
  stop("Failed to normalize top-level Samples metadata.")
}
# -----------------------------------------------------------------------------
# Section 7 · Save and read back the canonical processed object
# -----------------------------------------------------------------------------

saveRDS(canonical, OUTPUT_RDS, compress = TRUE)

if (!file.exists(OUTPUT_RDS) ||
    is.na(file.info(OUTPUT_RDS)$size) ||
    file.info(OUTPUT_RDS)$size <= 0) {
  stop("Failed to create a non-empty canonical processed object.")
}

canonical_check <- readRDS(OUTPUT_RDS)

if (!inherits(canonical_check, "SingleCellExperiment")) {
  stop("Saved canonical object did not read back as a SingleCellExperiment.")
}

if (!identical(
  dim(canonical_check),
  c(EXPECTED_N_GENES, EXPECTED_N_CANONICAL)
)) {
  stop("Saved canonical object failed the dimension read-back check.")
}

if (!identical(colnames(canonical_check), manifest$cell_id)) {
  stop("Saved canonical object failed the cell-order read-back check.")
}

if (!identical(
  as.character(canonical_check$patient),
  as.character(manifest$patient)
)) {
  stop("Saved canonical object failed the patient-metadata read-back check.")
}

if (!identical(as.character(canonical_check$tissue), manifest$tissue)) {
  stop("Saved canonical object failed the tissue-metadata read-back check.")
}

if (!identical(
  as.character(canonical_check$moca_group),
  manifest$moca_group
)) {
  stop("Saved canonical object failed the MoCA-metadata read-back check.")
}

rm(canonical_check)
gc()


# -----------------------------------------------------------------------------
# Section 8 · Closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PHASE 0 COMPLETE — CANONICAL INPUT RECONSTRUCTED\n")
cat("============================================================\n")
cat("  Public source            : GSE292141\n")
cat("  Analyzed PBMC pre-QC     : ", EXPECTED_N_RAW_PBMC, "\n", sep = "")
cat("  Historical-QC PBMC       : ", EXPECTED_N_QC_PBMC, "\n", sep = "")
cat("  Canonical monocytes      : ", EXPECTED_N_CANONICAL, "\n", sep = "")
cat("  Genes                    : ", EXPECTED_N_GENES, "\n", sep = "")
cat("  Patients                 : ", EXPECTED_N_PATIENTS, "\n", sep = "")
cat("  Manifest MD5             : ", manifest_md5, "\n", sep = "")
cat(
  "  Reconstructed RDS MD5   : ",
  unname(tools::md5sum(OUTPUT_RDS)),
  "\n",
  sep = ""
)
cat("  Output                   : data/processed/duy_monocytes_classical.rds\n")
cat("============================================================\n")
cat(
  "Validated contract: exact raw-count analysis population required by Script 01.\n"
)
cat("Next: Phase 1 — derive the iNPH IFN-activated monocyte signature.\n")
cat("============================================================\n")

# =============================================================================
# End of Phase 0.
# =============================================================================
