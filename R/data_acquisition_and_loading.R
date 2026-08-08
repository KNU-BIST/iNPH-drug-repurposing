# =============================================================================
# Data Acquisition and Loading Utilities
# Script: R/data_acquisition_and_loading.R
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
# Provide reproducible acquisition, validation, and loading utilities for the
# public datasets and metadata used by the iNPH drug-repurposing pipeline.
#
# Data sources
# ------------
# 1. GSE292141 (Duy et al., 2025)
#      iNPH PBMC and CSF single-cell RNA-sequencing data.
#
# 2. Perez et al. (2022) lupus atlas
#      SLE and healthy-control comparator data obtained from CZ CELLxGENE.
#
# 3. GSE70138
#      LINCS Phase II perturbagen metadata used to map BRD identifiers to
#      compound names.
#
# 4. compoundinfo_beta.txt (CLUE.io / CMap LINCS 2020)
#      Additional compound metadata used when GSE70138 does not provide a
#      human-readable compound name.
#
# LINCS2 expression signatures are not downloaded by this module. They are
# retrieved directly through Bioconductor ExperimentHub resource EH7297 by
# the connectivity-analysis scripts.
#
# Used by
# -------
# - 01_derive_inph_ifn_signature.R
# - 02_cross_disease_contextualization.R
# - 03_run_lincs2_connectivity.R
# - 05_generate_all_figures.R
# - 07_hallmark_ifn_sensitivity.R
#
# Outputs / side effects
# ----------------------
# - Downloads and stages raw or metadata files under project-defined paths.
# - Returns SingleCellExperiment objects from supported scRNA-seq inputs.
# - Validates required manually acquired metadata before downstream analysis.
#
# Reproducibility notes
# ---------------------
# - Stable public download endpoints are used where available.
# - CLUE.io compoundinfo_beta.txt requires a documented manual download step
#   because no stable direct-download endpoint is used by this project.
# - Functions fail explicitly when expected files are missing, empty, or
#   inconsistent with the required input structure.
# - GSE292141 sample metadata are defined once in load_duy_samples() and reused
#   by download_duy_samples() to avoid disagreement between acquisition and
#   loading steps.
# =============================================================================


# =============================================================================
# SOURCE 1 . GSE292141 (Duy et al. 2025) - iNPH PBMC + CSF atlas
# =============================================================================
# USED BY: 01_derive_inph_ifn_signature.R indirectly through the released
# canonical monocyte object, and 05_generate_all_figures.R for reconstruction
# of the PBMC-level embedding used in Supplementary Figure S1.
#
# GSE292141 raw counts are deposited as separate GSM-level supplementary
# barcodes/features/matrix triplets rather than as one series-level 10x bundle.
# The acquisition function therefore downloads each GSM supplementary set
# individually.

#' Download and extract all 20 Duy et al. 2025 (GSE292141) sample files
#'
#' Downloads each GSM's supplementary files individually (the correct
#' approach for this dataset) and extracts/renames them into the flat
#' directory structure load_duy_samples() expects:
#' `<GSM>_<sample_tag>_barcodes.tsv.gz`, `..._features.tsv.gz`,
#' `..._matrix.mtx.gz`.
#'
#' @param extracted_dir destination directory (matches load_duy_samples()'s
#'   extracted_dir argument).
#' @param sample_meta the sample table from load_duy_samples() - passed in
#'   rather than duplicated, so the two functions can never drift apart.
#' @return extracted_dir, invisibly.
download_duy_samples <- function(extracted_dir, sample_meta) {
  dir.create(extracted_dir, showWarnings = FALSE, recursive = TRUE)
  logger::log_info("Downloading {nrow(sample_meta)} GSM sample sets from GSE292141...")
  logger::log_info("(Per-GSM download, not per-series - see header note on why.)")

  for (i in seq_len(nrow(sample_meta))) {
    gsm <- sample_meta$gsm[i]
    tag <- sample_meta$sample_tag[i]
    prefix <- paste0(gsm, "_", tag)

    expected_files <- c(
      file.path(extracted_dir, paste0(prefix, "_barcodes.tsv.gz")),
      file.path(extracted_dir, paste0(prefix, "_features.tsv.gz")),
      file.path(extracted_dir, paste0(prefix, "_matrix.mtx.gz"))
    )

    file_info <- file.info(expected_files)

    already_present <- all(
      file.exists(expected_files) &
        !is.na(file_info$size) &
        file_info$size > 0
    )

    if (already_present) {
      logger::log_info("  {prefix}: already present, skipping.")
      next
    }

    logger::log_info("  Downloading {gsm} ({i}/{nrow(sample_meta)})...")
    tmp_dl <- file.path(tempdir(), gsm)
    if (dir.exists(tmp_dl)) {
      unlink(tmp_dl, recursive = TRUE, force = TRUE)
    }
    dir.create(tmp_dl, showWarnings = FALSE, recursive = TRUE)

    GEOquery::getGEOSuppFiles(
      GEO           = gsm,
      baseDir       = dirname(tmp_dl),
      makeDirectory = TRUE,
      fetch_files   = TRUE
    )

    downloaded <- list.files(tmp_dl, full.names = TRUE)
    if (length(downloaded) == 0) {
      stop("No supplementary files found for ", gsm,
           " - check the accession is correct and GEO is reachable.")
    }

    # GEO supplementary files for 10x data are typically one .tar bundling
    # the barcodes/features/matrix triplet - untar if so, otherwise assume
    # they're already loose files.
    tar_file <- downloaded[grepl("\\.tar$", downloaded)]
    if (length(tar_file) == 1) {
      untar(tar_file, exdir = tmp_dl)
      file.remove(tar_file)
      downloaded <- list.files(tmp_dl, full.names = TRUE, recursive = TRUE)
    }

    barcode_f <- downloaded[grepl("barcodes", downloaded, ignore.case = TRUE)]
    feature_f <- downloaded[grepl("features|genes", downloaded, ignore.case = TRUE)]
    matrix_f  <- downloaded[grepl("matrix", downloaded, ignore.case = TRUE)]

    if (length(barcode_f) != 1 || length(feature_f) != 1 || length(matrix_f) != 1) {
      stop("Expected exactly one barcodes/features/matrix file each for ", gsm,
           " - found ", length(barcode_f), "/", length(feature_f), "/",
           length(matrix_f), ". Inspect ", tmp_dl, " manually.")
    }

    copy_ok <- c(
      file.copy(barcode_f, file.path(extracted_dir, paste0(prefix, "_barcodes.tsv.gz")), overwrite = TRUE),
      file.copy(feature_f, file.path(extracted_dir, paste0(prefix, "_features.tsv.gz")), overwrite = TRUE),
      file.copy(matrix_f,  file.path(extracted_dir, paste0(prefix, "_matrix.mtx.gz")),   overwrite = TRUE)
    )

    if (!all(copy_ok)) {
      stop("Failed to stage one or more files for ", gsm, ".")
    }

    unlink(tmp_dl, recursive = TRUE, force = TRUE)
  }

  logger::log_info("Duy et al. 2025 download complete: {extracted_dir}")
  invisible(extracted_dir)
}


# =============================================================================
# SOURCE 2 . Perez et al. 2022 - lupus (SLE) comparator atlas
# =============================================================================
# USED BY: 02_cross_disease_contextualization.R for donor-level comparison
#   with SLE and healthy-control classical monocytes.
#
# PROVENANCE: Perez RK et al. (2022) Science. CZ CELLxGENE Discover
# collection 436154da-bcf1-4130-9c8b-120ff9a888f2
# (https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2).
# The direct file endpoint below corresponds to the Perez et al. dataset
# within this collection. CZ CELLxGENE may re-version dataset files over time.
#
# CAVEAT: CZ CELLxGENE periodically re-versions dataset files under new
# UUIDs. If this URL 404s, go to the collection page above, find the current
# "Perez et al. 2022" dataset entry, and update DATASET_URL below - do not
# guess at a replacement UUID.
#
# VERIFICATION AFTER DOWNLOAD: the downstream contextualization dataset should
# yield 307,429 classical monocytes after the study-specific cell-type filter.
# The source atlas contains substantially more cells across all cell types; any
# material change in the filtered classical-monocyte count should be reviewed
# before downstream analysis.

#' Download the Perez et al. 2022 lupus atlas (h5ad, ~12 GB)
#'
#' @param dest destination file path (should end in .h5ad).
#' @return dest, invisibly.
download_perez_lupus <- function(dest) {
  DATASET_URL <- "https://datasets.cellxgene.cziscience.com/8f2e4a29-e397-4f33-bc9a-3181e06b64b0.h5ad"

  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(dest)) {

    if (
      is.na(file.info(dest)$size) ||
      file.info(dest)$size <= 0
    ) {
      stop(
        "Existing Perez dataset file is empty or unreadable: ",
        dest
      )
    }

    if (
      requireNamespace("rhdf5", quietly = TRUE) &&
      !rhdf5::H5Fis_hdf5(dest)
    ) {
      stop(
        "Existing Perez dataset file is not a valid HDF5 file: ",
        dest
      )
    }

    logger::log_info(
      "Perez 2022 atlas already present and validated: {dest}"
    )

    return(invisible(dest))
  }

  logger::log_info("Downloading Perez et al. 2022 lupus atlas (~12 GB, will take a while)...")
  logger::log_info("Source: {DATASET_URL}")
  logger::log_warn("If this URL fails, check https://cellxgene.cziscience.com/",
                   "collections/436154da-bcf1-4130-9c8b-120ff9a888f2 for the ",
                   "current dataset UUID - CZ CELLxGENE re-versions files ",
                   "periodically.")

  utils::download.file(DATASET_URL, destfile = dest, mode = "wb",
                       method = "libcurl")

  if (!file.exists(dest) || is.na(file.info(dest)$size) || file.info(dest)$size == 0) {
    stop("Perez dataset download failed or produced an empty file: ", dest)
  }

  if (requireNamespace("rhdf5", quietly = TRUE) && !rhdf5::H5Fis_hdf5(dest)) {
    stop("Downloaded Perez file is not a valid HDF5 file: ", dest)
  }

  logger::log_info("Download complete: {dest}")
  logger::log_warn(
    "VERIFY before downstream analysis: the source atlas contains ",
    "1,263,676 cells across all cell types, and the finalized ",
    "classical-monocyte filter should yield 307,429 cells. ",
    "See 02_cross_disease_contextualization.R for the study-specific ",
    "cell-count checks.")
  invisible(dest)
}


# =============================================================================
# SOURCE 3 . GSE70138 - LINCS Phase II pert_info (compound metadata)
# =============================================================================
# USED BY: 03_run_lincs2_connectivity.R and 07_hallmark_ifn_sensitivity.R -
#   resolving BRD compound IDs to human-readable drug names.
#
# This IS a standard series-level GEO supplementary file (unlike GSE292141),
# so the generic getGEOSuppFiles() approach works correctly here.

#' Download the GSE70138 compound/perturbagen info file
#'
#' @param dest_dir destination directory.
#' @return path to the downloaded pert_info file, invisibly.
download_lincs_pert_info <- function(dest_dir) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  expected_file <- file.path(dest_dir, "GSE70138_Broad_LINCS_pert_info.txt.gz")

  if (
    file.exists(expected_file) &&
    !is.na(file.info(expected_file)$size) &&
    file.info(expected_file)$size > 0
  ) {
    logger::log_info("GSE70138 pert_info already present: {expected_file}")
    return(invisible(expected_file))
  }

  logger::log_info("Downloading GSE70138 supplementary files (compound metadata)...")
  GEOquery::getGEOSuppFiles(
    GEO           = "GSE70138",
    baseDir       = dest_dir,
    makeDirectory = FALSE,
    fetch_files   = TRUE
  )

  if (
    !file.exists(expected_file) ||
    is.na(file.info(expected_file)$size) ||
    file.info(expected_file)$size == 0
  ) {
    stop(
      "Expected non-empty file not found after download: ",
      expected_file,
      ". Check GSE70138's current supplementary file listing on GEO; ",
      "the filename may have changed."
    )
  }
  invisible(expected_file)
}


# =============================================================================
# SOURCE 4 . compoundinfo_beta.txt (CLUE.io LINCS 2020 beta release)
# =============================================================================
# USED BY: 03_run_lincs2_connectivity.R and 07_hallmark_ifn_sensitivity.R -
#   provides additional compound-name coverage beyond the GSE70138
#   perturbagen metadata.
#
# PROVENANCE: CLUE.io LINCS Data Releases dashboard
# (https://clue.io/releases/data-dashboard), CMap/LINCS 2020 release.
#
# NOT AUTOMATED: this project does not rely on a stable scriptable endpoint
# for this file. The required manual acquisition step is documented below.
#
# MANUAL DOWNLOAD STEPS:
#   1. Go to https://clue.io/releases/data-dashboard
#   2. Find the "LINCS 2020" / "CMap LINCS 2020" release section
#   3. Download compoundinfo_beta.txt (and, if available, the companion
#      compoundinfo_beta.xlsx)
#   4. Place both at: data/raw/compoundinfo_beta.txt (and .xlsx)

#' Check for compoundinfo_beta.txt and explain how to obtain it if missing
#'
#' Does not attempt automated download - see header note. This function
#' exists so the pipeline fails with a clear, actionable message instead of
#' a confusing downstream error if the file is missing.
#'
#' @param expected_path where this pipeline expects the file
#'   (data/raw/compoundinfo_beta.txt).
check_compoundinfo_beta <- function(expected_path) {
  if (
    file.exists(expected_path) &&
    !is.na(file.info(expected_path)$size) &&
    file.info(expected_path)$size > 0
  ) {
    logger::log_info("compoundinfo_beta.txt found: {expected_path}")
    return(invisible(TRUE))
  }

  stop(
    "compoundinfo_beta.txt not found at: ", expected_path, "\n\n",
    "This file has no automated download in this pipeline - CLUE.io's data ",
    "dashboard does not expose a stable direct-download URL. To obtain it:\n",
    "  1. Go to https://clue.io/releases/data-dashboard\n",
    "  2. Find the 'LINCS 2020' / 'CMap LINCS 2020' release section\n",
    "  3. Download compoundinfo_beta.txt\n",
    "  4. Place it at: ", expected_path, "\n"
  )
}


# =============================================================================
# LINCS2 L1000 DATA (EH7297) - NOTE, NOT A FUNCTION
# =============================================================================
# LINCS2 signatures are obtained directly through Bioconductor ExperimentHub
# by 03_run_lincs2_connectivity.R and 07_hallmark_ifn_sensitivity.R:
#
#     eh <- ExperimentHub::ExperimentHub(ask = FALSE)
#     lincs2_path <- eh[["EH7297"]]
#
# No separate LINCS2 download wrapper is required in this utility module.


# =============================================================================
# GENERIC LOADING HELPERS
# These utilities support standardized loading and subsetting of common
# single-cell object formats.
# =============================================================================

#' Load scRNA-seq data into a SingleCellExperiment
#'
#' Detects file type (10x H5, MTX, RDS, h5ad) and dispatches to the right
#' loader. Returns a [SingleCellExperiment].
#'
#' @param path path to the data directory or file.
load_scrnaseq <- function(path) {
  logger::log_info("Loading scRNA-seq data from: {path}")

  if (!file.exists(path)) {
    stop("Path does not exist: ", path)
  }

  if (!dir.exists(path)) {
    extension <- tolower(tools::file_ext(path))

    if (extension == "h5ad") {
      return(zellkonverter::readH5AD(path))
    }

    if (extension == "rds") {
      object <- readRDS(path)
      if (inherits(object, "Seurat")) {
        return(Seurat::as.SingleCellExperiment(object))
      }
      return(object)
    }

    if (extension == "h5") {
      return(DropletUtils::read10xCounts(path))
    }

    stop("Unsupported scRNA-seq file format: ", path)
  }

  h5ad_files <- list.files(path, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) > 0L) {
    return(zellkonverter::readH5AD(h5ad_files[[1]]))
  }

  rds_files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) > 0L) {
    object <- readRDS(rds_files[[1]])
    if (inherits(object, "Seurat")) {
      return(Seurat::as.SingleCellExperiment(object))
    }
    return(object)
  }

  h5_files <- list.files(path, pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) > 0L) {
    return(DropletUtils::read10xCounts(h5_files[[1]]))
  }

  if (file.exists(file.path(path, "matrix.mtx")) ||
      file.exists(file.path(path, "matrix.mtx.gz"))) {
    return(DropletUtils::read10xCounts(path))
  }

  stop("Unrecognized scRNA-seq format in: ", path)
}

#' Subset a SingleCellExperiment to one cell type
#'
#' Looks for a biological cell-type label in `colData` columns named, in
#' priority order: `cell_type`, `celltype`, `predicted.celltype.l2`, and
#' `annotation`.
#'
#' @param sce [SingleCellExperiment].
#' @param cell_type character label to keep (case-insensitive substring match).
subset_cell_type <- function(sce, cell_type = "monocyte") {
  candidate_cols <- c(
    "cell_type",
    "celltype",
    "predicted.celltype.l2",
    "annotation"
  )
  use <- intersect(candidate_cols, colnames(SummarizedExperiment::colData(sce)))
  if (length(use) == 0) {
    stop("No cell type column found in colData. Available: ",
         paste(colnames(SummarizedExperiment::colData(sce)), collapse = ", "))
  }
  col <- use[1]
  labels <- SummarizedExperiment::colData(sce)[[col]]
  keep   <- grepl(cell_type, labels, ignore.case = TRUE)
  logger::log_info("Subsetting to '{cell_type}' using column '{col}': ",
                   "{sum(keep)} / {length(keep)} cells")
  sce[, keep]
}


# =============================================================================
# GSE292141 SAMPLE LOADER
# The verified sample_meta table below is the single source of truth for both
# loading and acquisition, so download_duy_samples() and load_duy_samples()
# use the same GSM, tissue, patient, and MoCA annotations.
# =============================================================================

#' Load Duy GSE292141 multi-sample 10X data into one SingleCellExperiment
#'
#' Reads all 20 sample triplets (barcodes/features/matrix), tags each cell
#' with patient ID, tissue (CSF/PBMC), and MOCA group (High/Low/Unknown),
#' then concatenates into one SCE. Optionally filters to one tissue. If the
#' expected files are not present in extracted_dir, downloads them first via
#' download_duy_samples().
#'
#' @param extracted_dir directory containing the per-sample mtx triplets.
#' @param tissue_filter "CSF", "PBMC", or NULL for both.
#' @param drop_unknown_moca drop the Patient 6 Unknown-MOCA samples.
#' @return SingleCellExperiment with colData columns: sample_id, patient,
#'   tissue, moca_group.
load_duy_samples <- function(
    extracted_dir,
    tissue_filter = NULL,
    drop_unknown_moca = TRUE
) {

  # ---- Hard-coded sample table (from GEO metadata) ------------------------
  # Verified against phenoData of GSE292141 — May 2025 deposit.
  sample_meta_full <- tibble::tribble(
    ~gsm,           ~sample_tag,    ~patient, ~tissue, ~moca_group,
    "GSM8851048",   "1C",            1,        "CSF",   "High",
    "GSM8851049",   "1B",            1,        "PBMC",  "High",
    "GSM8851050",   "2C",            2,        "CSF",   "High",
    "GSM8851051",   "2B",            2,        "PBMC",  "High",
    "GSM8851052",   "3C",            3,        "CSF",   "High",
    "GSM8851053",   "3B",            3,        "PBMC",  "High",
    "GSM8851054",   "4C",            4,        "CSF",   "High",
    "GSM8851055",   "4B",            4,        "PBMC",  "High",
    "GSM8851056",   "5C",            5,        "CSF",   "High",
    "GSM8851057",   "5B",            5,        "PBMC",  "High",
    "GSM8851058",   "6C",            6,        "CSF",   "Unknown",
    "GSM8851059",   "6B",            6,        "PBMC",  "Unknown",
    "GSM8851060",   "101822C",       7,        "CSF",   "Low",
    "GSM8851061",   "101822B",       7,        "PBMC",  "Low",
    "GSM8851062",   "102522C2",      8,        "CSF",   "Low",
    "GSM8851063",   "102522B2",      8,        "PBMC",  "Low",
    "GSM8851064",   "102522C3",      9,        "CSF",   "Low",
    "GSM8851065",   "102522B3",      9,        "PBMC",  "Low",
    "GSM8851066",   "110122C",       10,       "CSF",   "Low",
    "GSM8851067",   "110122B",       10,       "PBMC",  "Low"
  )

  # Download step: verify the exact three expected files for every GSM before
  # applying tissue/MoCA filters. This prevents stale or unrelated files with
  # the same prefix from being counted as a complete sample set.
  sample_set_complete <- vapply(
    seq_len(nrow(sample_meta_full)),
    function(i) {
      prefix <- paste0(
        sample_meta_full$gsm[i],
        "_",
        sample_meta_full$sample_tag[i]
      )

      expected_files <- c(
        file.path(extracted_dir, paste0(prefix, "_barcodes.tsv.gz")),
        file.path(extracted_dir, paste0(prefix, "_features.tsv.gz")),
        file.path(extracted_dir, paste0(prefix, "_matrix.mtx.gz"))
      )

      all(
        file.exists(expected_files) &
          !is.na(file.info(expected_files)$size) &
          file.info(expected_files)$size > 0
      )
    },
    logical(1)
  )

  n_present <- sum(sample_set_complete)

  if (n_present < nrow(sample_meta_full)) {
    logger::log_info(
      "{n_present}/{nrow(sample_meta_full)} GSM sample sets present; ",
      "downloading missing or incomplete sets."
    )

    download_duy_samples(
      extracted_dir,
      sample_meta_full
    )
  }

  # ---- Apply filters early --------------------------------------------------
  sample_meta <- sample_meta_full
  if (!is.null(tissue_filter)) {
    sample_meta <- sample_meta |>
      dplyr::filter(tissue == tissue_filter)
  }
  if (drop_unknown_moca) {
    sample_meta <- sample_meta |>
      dplyr::filter(moca_group != "Unknown")
  }
  sample_meta$moca_group <- factor(sample_meta$moca_group,
                                   levels = c("High", "Low"))

  logger::log_info("Loading {nrow(sample_meta)} samples from Duy GSE292141.")

  # ---- Load each sample ---------------------------------------------------
  load_one <- function(row) {
    prefix <- paste0(row$gsm, "_", row$sample_tag)

    source_files <- c(
      barcodes = file.path(
        extracted_dir,
        paste0(prefix, "_barcodes.tsv.gz")
      ),
      features = file.path(
        extracted_dir,
        paste0(prefix, "_features.tsv.gz")
      ),
      matrix = file.path(
        extracted_dir,
        paste0(prefix, "_matrix.mtx.gz")
      )
    )

    source_ok <- file.exists(source_files) &
      !is.na(file.info(source_files)$size) &
      file.info(source_files)$size > 0

    if (!all(source_ok)) {
      stop(
        "Missing or empty 10x input file(s) for ",
        prefix,
        ": ",
        paste(
          names(source_files)[!source_ok],
          collapse = ", "
        )
      )
    }

    # Stage into a clean temporary directory with standard Cell Ranger names.
    tmp <- file.path(tempdir(), prefix)

    if (dir.exists(tmp)) {
      unlink(tmp, recursive = TRUE, force = TRUE)
    }

    dir.create(
      tmp,
      recursive = TRUE,
      showWarnings = FALSE
    )

    copy_ok <- c(
      file.copy(
        source_files[["barcodes"]],
        file.path(tmp, "barcodes.tsv.gz"),
        overwrite = TRUE
      ),
      file.copy(
        source_files[["features"]],
        file.path(tmp, "features.tsv.gz"),
        overwrite = TRUE
      ),
      file.copy(
        source_files[["matrix"]],
        file.path(tmp, "matrix.mtx.gz"),
        overwrite = TRUE
      )
    )

    if (!all(copy_ok)) {
      unlink(tmp, recursive = TRUE, force = TRUE)
      stop("Failed to stage one or more 10x files for ", prefix, ".")
    }

    sce <- DropletUtils::read10xCounts(
      tmp,
      col.names = TRUE
    )

    unlink(
      tmp,
      recursive = TRUE,
      force = TRUE
    )

    # Tag colData with sample-level metadata.
    sce$sample_id  <- row$gsm
    sce$patient    <- row$patient
    sce$tissue     <- row$tissue
    sce$moca_group <- row$moca_group

    # Make cell barcodes globally unique across samples.
    colnames(sce) <- paste0(
      row$gsm,
      "_",
      colnames(sce)
    )

    sce
  }

  sce_list <- purrr::map(seq_len(nrow(sample_meta)),
                         ~ load_one(sample_meta[.x, ]),
                         .progress = TRUE)

  # ---- Concatenate ---------------------------------------------------------
  logger::log_info("Concatenating {length(sce_list)} samples into one SCE.")
  merged <- do.call(SingleCellExperiment::cbind, sce_list)
  logger::log_info("Final SCE: {nrow(merged)} genes × {ncol(merged)} cells.")
  merged
}
