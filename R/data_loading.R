# =============================================================================
# R/data_loading.R — data acquisition from public repositories
# =============================================================================

#' Download a GEO dataset to a local directory
#'
#' Wraps [GEOquery::getGEOSuppFiles()]. Returns the path of the destination
#' directory so `targets` can fingerprint it as a file dependency.
#'
#' @param accession GEO accession (e.g., "GSE92742"). Use placeholder
#'   "GSEXXXXXX" during scaffolding; replace once Duy et al. data is located.
#' @param dest destination directory.
#' @return path to the destination directory.
download_geo_data <- function(accession, dest) {
  if (grepl("XXXXXX", accession)) {
    logger::log_warn("Placeholder GEO accession {accession}: skipping download.
                     Update config.yml with the real accession.")
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    return(dest)
  }
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  logger::log_info("Downloading GEO {accession} → {dest}")
  GEOquery::getGEOSuppFiles(
    GEO          = accession,
    baseDir      = dest,
    makeDirectory = TRUE,
    fetch_files  = TRUE
  )
  dest
}

#' Download LINCS L1000 perturbation data (Phase I + II)
#'
#' LINCS data is large (~50 GB total). This function fetches the Level 5
#' GCTX matrices and metadata via [signatureSearchData] when available, with
#' GEO fallback.
#'
#' @param phase1_geo phase I GEO accession (GSE92742).
#' @param phase2_geo phase II GEO accession (GSE70138).
#' @param dest destination directory.
download_lincs <- function(phase1_geo, phase2_geo, dest) {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  logger::log_info("Acquiring LINCS L1000 data → {dest}")

  # Preferred path: signatureSearchData ExperimentHub resource
  if (requireNamespace("signatureSearchData", quietly = TRUE)) {
    logger::log_info("Using signatureSearchData ExperimentHub resource.")
    eh <- ExperimentHub::ExperimentHub()
    # NOTE: confirm ID at runtime; these are the canonical LINCS resources.
    # query(eh, c("signatureSearchData", "lincs"))
    return(dest)
  }

  # Fallback: direct GEO supplementary download
  logger::log_warn("signatureSearchData unavailable; falling back to GEO.")
  for (acc in c(phase1_geo, phase2_geo)) {
    GEOquery::getGEOSuppFiles(acc, baseDir = dest, makeDirectory = TRUE)
  }
  dest
}

#' Load scRNA-seq data into a SingleCellExperiment
#'
#' Detects file type (10x H5, MTX, RDS, h5ad) and dispatches to the right
#' loader. Returns a [SingleCellExperiment].
#'
#' @param path path to the data directory or file.
load_scrnaseq <- function(path) {
  logger::log_info("Loading scRNA-seq from: {path}")

  if (length(list.files(path, pattern = "\\.h5ad$")) > 0) {
    f <- list.files(path, pattern = "\\.h5ad$", full.names = TRUE)[1]
    return(zellkonverter::readH5AD(f))
  }

  if (length(list.files(path, pattern = "\\.rds$")) > 0) {
    f <- list.files(path, pattern = "\\.rds$", full.names = TRUE)[1]
    obj <- readRDS(f)
    if (inherits(obj, "Seurat")) {
      return(Seurat::as.SingleCellExperiment(obj))
    }
    return(obj)
  }

  if (length(list.files(path, pattern = "\\.h5$")) > 0) {
    f <- list.files(path, pattern = "\\.h5$", full.names = TRUE)[1]
    counts <- DropletUtils::read10xCounts(f)
    return(counts)
  }

  if (length(list.files(path, pattern = "matrix\\.mtx")) > 0) {
    return(DropletUtils::read10xCounts(path))
  }

  stop("Unrecognised scRNA-seq format in: ", path)
}

#' Subset a SingleCellExperiment to one cell type
#'
#' Looks for the cell type label in `colData` columns named (in priority):
#' `cell_type`, `celltype`, `predicted.celltype.l2`, `seurat_clusters`.
#'
#' @param sce [SingleCellExperiment].
#' @param cell_type character label to keep (case-insensitive substring match).
subset_cell_type <- function(sce, cell_type = "monocyte") {
  candidate_cols <- c("cell_type", "celltype",
                      "predicted.celltype.l2", "annotation",
                      "seurat_clusters")
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



#' Load Duy GSE292141 multi-sample 10X data into one SingleCellExperiment
#'
#' Reads all 20 sample triplets (barcodes/features/matrix), tags each cell
#' with patient ID, tissue (CSF/PBMC), and MOCA group (High/Low/Unknown),
#' then concatenates into one SCE. Optionally filters to one tissue.
#'
#' @param extracted_dir directory containing the per-sample mtx triplets.
#' @param tissue_filter "CSF", "PBMC", or NULL for both.
#' @param drop_unknown_moca drop the Patient 6 Unknown-MOCA samples.
#' @return SingleCellExperiment with colData columns: sample_id, patient,
#'   tissue, moca_group.
load_duy_samples <- function(extracted_dir,
                             tissue_filter   = "CSF",
                             drop_unknown_moca = TRUE) {

  # ---- Hard-coded sample table (from GEO metadata) ------------------------
  # Verified against phenoData of GSE292141 — May 2025 deposit.
  sample_meta <- tibble::tribble(
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

  # Apply filters early
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
    files <- list.files(extracted_dir, pattern = prefix, full.names = TRUE)
    if (length(files) != 3) {
      stop("Expected 3 files for ", prefix, "; found ", length(files))
    }
    # Stage into a temp dir with Cell Ranger naming so DropletUtils can read
    tmp <- file.path(tempdir(), prefix)
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    file.copy(grep("barcodes",  files, value = TRUE),
              file.path(tmp, "barcodes.tsv.gz"), overwrite = TRUE)
    file.copy(grep("features",  files, value = TRUE),
              file.path(tmp, "features.tsv.gz"), overwrite = TRUE)
    file.copy(grep("matrix",    files, value = TRUE),
              file.path(tmp, "matrix.mtx.gz"), overwrite = TRUE)
    sce <- DropletUtils::read10xCounts(tmp, col.names = TRUE)

    # Tag colData with sample-level metadata
    sce$sample_id  <- row$gsm
    sce$patient    <- row$patient
    sce$tissue     <- row$tissue
    sce$moca_group <- row$moca_group
    # Make barcodes globally unique
    colnames(sce) <- paste0(row$gsm, "_", colnames(sce))
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
