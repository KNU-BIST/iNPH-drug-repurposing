# =============================================================================
# Revision analysis 05B — Ventricular-CSF monocyte IFN contextualization
# Install as analysis/05B_csf_monocyte_ifn_contextualization.R
# Build: 05B_CSF_MONOCYTE_IFN_CONTEXTUALIZATION_V1_20260907
# =============================================================================
#
# Purpose
# -------
# Test whether the blood-derived IFN-response-high monocyte program is
# detectable in independently selected ventricular-CSF classical-monocyte
# candidates from the same GSE292141 cohort. This directly addresses the
# reviewer request to strengthen the blood/CSF biological connection.
#
# Interpretation boundary
# -----------------------
# - CSF monocytes are selected without using IFN genes or IFN scores.
# - The frozen AUCell method, gene sets, aucMaxRank, and threshold are reused.
# - The primary paired comparison uses the nine donors in the canonical PBMC
#   analysis; Patient 6 is retained descriptively because MoCA is unknown.
# - Blood-versus-CSF differences are exploratory compartment comparisons.
# - This analysis does not establish iNPH specificity, CNS causality,
#   cognitive association, or therapeutic efficacy.
# - No LINCS query is executed and no existing output is overwritten.
#
# Inputs
# ------
# data/raw/duy/extracted/<all GSE292141 10x triplets>
# data/processed/inph_ifn_signature.rds
# data/processed/isg_reference_panel.rds
# data/processed/phase2_contextualization.rds
# R/data_acquisition_and_loading.R
# R/pbmc_clustering_and_annotation.R
#
# Outputs
# -------
# data/processed/revision_csf_classical_monocytes.rds
# data/processed/revision_csf_ifn_contextualization.rds
# results/revision/csf_ifn_contextualization/*.csv and reporting notes
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment and frozen configuration
# -----------------------------------------------------------------------------

BUILD <- "05B_CSF_MONOCYTE_IFN_CONTEXTUALIZATION_V1_20260907"

required_packages <- c(
  "here",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "Matrix",
  "Seurat",
  "DropletUtils",
  "scuttle",
  "harmony",
  "AUCell",
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
  library(S4Vectors)
  library(Matrix)
})

RAW_DIR <- here::here("data", "raw")
PROCESSED_DIR <- here::here("data", "processed")
REVISION_RESULTS_DIR <- here::here(
  "results",
  "revision",
  "csf_ifn_contextualization"
)

DUY_EXTRACTED <- file.path(RAW_DIR, "duy", "extracted")
DATA_HELPER <- here::here("R", "data_acquisition_and_loading.R")
CLUSTER_HELPER <- here::here("R", "pbmc_clustering_and_annotation.R")
SIGNATURE_PATH <- file.path(PROCESSED_DIR, "inph_ifn_signature.rds")
ISG_REFERENCE_PATH <- file.path(PROCESSED_DIR, "isg_reference_panel.rds")
PHASE2_PATH <- file.path(PROCESSED_DIR, "phase2_contextualization.rds")

OUTPUT_MONOCYTES_RDS <- file.path(
  PROCESSED_DIR,
  "revision_csf_classical_monocytes.rds"
)
OUTPUT_RESULT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_csf_ifn_contextualization.rds"
)
OUTPUT_CELL_SCORES_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_ifn_cell_scores.csv"
)
OUTPUT_DONOR_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_ifn_donor_summary.csv"
)
OUTPUT_PAIRED_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_pbmc_paired_donor_summary.csv"
)
OUTPUT_PAIRED_TESTS_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_pbmc_paired_tests.csv"
)
OUTPUT_CLUSTER_DIAGNOSTICS_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_classical_monocyte_cluster_diagnostics.csv"
)
OUTPUT_ISG_DETECTION_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_isg_detection_summary.csv"
)
OUTPUT_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_ifn_contextualization_summary.csv"
)
OUTPUT_RAW_INVENTORY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "csf_raw_input_inventory.csv"
)
OUTPUT_REPORTING_NOTES <- file.path(
  REVISION_RESULTS_DIR,
  "csf_contextualization_reporting_notes.txt"
)

OUTPUT_PATHS <- c(
  OUTPUT_MONOCYTES_RDS,
  OUTPUT_RESULT_RDS,
  OUTPUT_CELL_SCORES_CSV,
  OUTPUT_DONOR_SUMMARY_CSV,
  OUTPUT_PAIRED_SUMMARY_CSV,
  OUTPUT_PAIRED_TESTS_CSV,
  OUTPUT_CLUSTER_DIAGNOSTICS_CSV,
  OUTPUT_ISG_DETECTION_CSV,
  OUTPUT_SUMMARY_CSV,
  OUTPUT_RAW_INVENTORY_CSV,
  OUTPUT_REPORTING_NOTES
)

EXPECTED_GENES <- 36601L
EXPECTED_ALL_PATIENTS <- as.character(1:10)
EXPECTED_PRIMARY_PATIENTS <- c(
  "1", "2", "3", "4", "5", "7", "8", "9", "10"
)
EXPECTED_PRIMARY_PBMC_CELLS <- 22197L
EXPECTED_PRIMARY_PBMC_DONORS <- 9L
EXPECTED_ISG_SELECTOR <- c(
  "ISG15", "IFI6", "IFI44", "IFI44L", "MX1",
  "MX2", "OAS1", "OAS2", "OAS3", "IFIT1",
  "IFIT3", "RSAD2", "STAT1", "USP18", "LY6E"
)

MIN_GENES_PER_CELL <- 200L
MAX_PCT_MT <- 20
CSF_RESOLUTION <- 0.5
N_PCS <- 30L
CLUSTER_ALGORITHM <- 1L
CLUSTER_SEED <- 0L
PCA_SEED <- 42L
UMAP_SEED <- 42L
MONOCYTE_RELATIVE_THRESHOLD <- 0.75
SCORING_SEED <- 1L

# The same classical panel used in Script 01 and Analysis 01C. It contains no
# frozen IFN-selector gene and is the only panel used to select CSF candidates.
CLASSICAL_MONOCYTE_SELECTION_PANEL <- c(
  "CD14", "S100A8", "S100A9", "S100A12", "VCAN", "FCN1"
)

DIAGNOSTIC_PANELS <- list(
  classical_monocyte = CLASSICAL_MONOCYTE_SELECTION_PANEL,
  nonclassical_monocyte = c(
    "FCGR3A", "MS4A7", "CDKN1C", "LST1", "CSF1R", "RHOC"
  ),
  T_cell = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R"),
  B_cell = c("CD79A", "MS4A1", "CD37", "CD74", "CD19"),
  NK_cell = c("NKG7", "GNLY", "PRF1", "KLRD1", "CTSW"),
  platelet = c("PPBP", "PF4", "NRGN", "GNG11"),
  dendritic = c("FCER1A", "CLEC10A", "GZMB", "CLEC4C"),
  glial_or_astrocytic = c("GFAP", "AQP4", "SLC1A3", "PLP1", "MOG")
)

if (length(intersect(
  CLASSICAL_MONOCYTE_SELECTION_PANEL,
  EXPECTED_ISG_SELECTOR
)) > 0L) {
  stop("The CSF selection panel unexpectedly overlaps the IFN selector.")
}

required_files <- c(
  DATA_HELPER,
  CLUSTER_HELPER,
  SIGNATURE_PATH,
  ISG_REFERENCE_PATH,
  PHASE2_PATH
)
required_info <- file.info(required_files)
bad_required <- required_files[
  !file.exists(required_files) |
    is.na(required_info$size) |
    required_info$size <= 0
]

if (length(bad_required) > 0L) {
  stop(
    "Required CSF-contextualization input(s) missing or empty:\n  ",
    paste(bad_required, collapse = "\n  ")
  )
}

if (!dir.exists(DUY_EXTRACTED)) {
  stop("Duy extracted-data directory is missing: ", DUY_EXTRACTED)
}

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_CSF_CONTEXTUALIZATION", unset = "0"),
  "1"
)
existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "CSF-contextualization output(s) already exist:\n  ",
    paste(existing_outputs, collapse = "\n  "),
    "\nRefusing to overwrite automatically. If replacement is intentional, ",
    "set INPH_OVERWRITE_CSF_CONTEXTUALIZATION=1 for that run."
  )
}

source(DATA_HELPER)
source(CLUSTER_HELPER)


# -----------------------------------------------------------------------------
# Section 1 · Validate the frozen blood-derived state definition
# -----------------------------------------------------------------------------

cat("\n--- Validating the frozen blood-derived IFN state definition ---\n")

signature <- readRDS(SIGNATURE_PATH)
phase2 <- readRDS(PHASE2_PATH)
expected_isgs <- as.character(readRDS(ISG_REFERENCE_PATH))

if (!is.list(signature) ||
    !all(c("up", "down") %in% names(signature)) ||
    length(signature$up) != 150L ||
    length(signature$down) != 150L ||
    anyDuplicated(signature$up) ||
    anyDuplicated(signature$down)) {
  stop("The frozen primary signature is not the expected 150-UP/150-DOWN set.")
}

required_phase2_fields <- c(
  "engine_A_primary",
  "inph_cell_scores",
  "params",
  "excluded_genes",
  "input_provenance"
)
missing_phase2_fields <- setdiff(required_phase2_fields, names(phase2))

if (length(missing_phase2_fields) > 0L) {
  stop(
    "Phase 2 result is missing field(s): ",
    paste(missing_phase2_fields, collapse = ", ")
  )
}

if (!identical(as.character(phase2$excluded_genes), EXPECTED_ISG_SELECTOR)) {
  stop("The frozen Phase 2 IFN selector differs from the expected panel.")
}

AUC_MAXRANK <- as.integer(phase2$params$auc_maxrank)
IFN_THRESHOLD <- as.numeric(phase2$engine_A_primary$thr_shared)

if (length(AUC_MAXRANK) != 1L ||
    is.na(AUC_MAXRANK) ||
    AUC_MAXRANK != 1500L) {
  stop("Expected the frozen AUCell aucMaxRank to equal 1500.")
}

if (length(IFN_THRESHOLD) != 1L || !is.finite(IFN_THRESHOLD)) {
  stop("The frozen AUCell IFN threshold is missing or non-finite.")
}

if (length(expected_isgs) != 35L || anyDuplicated(expected_isgs)) {
  stop("The saved biological-validation ISG panel is not the expected 35 genes.")
}

canonical_scores <- as.data.frame(phase2$inph_cell_scores)
required_canonical_columns <- c(
  "cell_id", "donor", "A_score", "A_isg"
)

if (!all(required_canonical_columns %in% names(canonical_scores)) ||
    nrow(canonical_scores) != EXPECTED_PRIMARY_PBMC_CELLS ||
    anyDuplicated(canonical_scores$cell_id) ||
    length(unique(as.character(canonical_scores$donor))) !=
      EXPECTED_PRIMARY_PBMC_DONORS) {
  stop("The saved canonical PBMC cell-score table is invalid.")
}

if (!setequal(
  unique(as.character(canonical_scores$donor)),
  EXPECTED_PRIMARY_PATIENTS
)) {
  stop("The canonical PBMC scores do not represent the expected nine donors.")
}

signature_up_exact <- setdiff(
  as.character(signature$up),
  EXPECTED_ISG_SELECTOR
)
signature_down_exact <- as.character(signature$down)

selection_overlap_up <- intersect(
  signature_up_exact,
  CLASSICAL_MONOCYTE_SELECTION_PANEL
)
selection_overlap_down <- intersect(
  signature_down_exact,
  CLASSICAL_MONOCYTE_SELECTION_PANEL
)

signature_up_independent <- setdiff(
  signature_up_exact,
  CLASSICAL_MONOCYTE_SELECTION_PANEL
)
signature_down_independent <- setdiff(
  signature_down_exact,
  CLASSICAL_MONOCYTE_SELECTION_PANEL
)

cat(
  "Frozen threshold:", format(IFN_THRESHOLD, digits = 8),
  "| aucMaxRank:", AUC_MAXRANK, "\n"
)
cat(
  "Transfer signature:", length(signature_up_exact), "UP /",
  length(signature_down_exact), "DOWN genes\n"
)
cat(
  "Selection-panel overlap removed in independence sensitivity:",
  length(selection_overlap_up), "UP /",
  length(selection_overlap_down), "DOWN genes\n"
)


# -----------------------------------------------------------------------------
# Section 2 · Load all ten ventricular-CSF samples
# -----------------------------------------------------------------------------

cat("\n--- Loading all ten GSE292141 ventricular-CSF samples ---\n")

csf_raw <- load_duy_samples(
  extracted_dir = DUY_EXTRACTED,
  tissue_filter = "CSF",
  drop_unknown_moca = FALSE
)

if (!inherits(csf_raw, "SingleCellExperiment")) {
  stop("load_duy_samples() did not return a SingleCellExperiment.")
}

if (nrow(csf_raw) != EXPECTED_GENES) {
  stop(
    "Expected ", EXPECTED_GENES,
    " raw features; observed ", nrow(csf_raw), "."
  )
}

if (anyDuplicated(colnames(csf_raw))) {
  stop("The reconstructed CSF object contains duplicate cell identifiers.")
}

if (!all(as.character(csf_raw$tissue) == "CSF")) {
  stop("The CSF reconstruction contains non-CSF cells.")
}

if (!setequal(
  unique(as.character(csf_raw$patient)),
  EXPECTED_ALL_PATIENTS
)) {
  stop("The CSF reconstruction does not contain Patients 1–10.")
}

# Restore the documented Unknown label lost when the loader applies the
# High/Low factor levels while retaining Patient 6.
moca_group <- as.character(csf_raw$moca_group)
moca_group[as.character(csf_raw$patient) == "6"] <- "Unknown"
csf_raw$moca_group <- factor(
  moca_group,
  levels = c("High", "Low", "Unknown")
)

if (anyNA(csf_raw$moca_group)) {
  stop("Missing MoCA-group metadata remain after restoring Patient 6.")
}

# Normalize loader-generated temporary path metadata before any saved output.
if ("Sample" %in% colnames(SummarizedExperiment::colData(csf_raw))) {
  csf_raw$Sample <- as.character(csf_raw$sample_id)
}
csf_metadata <- S4Vectors::metadata(csf_raw)
csf_metadata$Samples <- unique(as.character(csf_raw$sample_id))
S4Vectors::metadata(csf_raw) <- csf_metadata

csf_raw <- fix_gene_symbols(csf_raw)

if (nrow(csf_raw) != EXPECTED_GENES || anyDuplicated(rownames(csf_raw))) {
  stop("CSF gene-symbol harmonization produced an invalid feature structure.")
}

raw_counts_by_patient <- table(
  factor(as.character(csf_raw$patient), levels = EXPECTED_ALL_PATIENTS)
)

cat(
  "Raw ventricular CSF:", nrow(csf_raw), "genes x",
  ncol(csf_raw), "cells\n"
)
print(raw_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 3 · Uniform CSF QC, integration, and clustering
# -----------------------------------------------------------------------------

cat("\n--- Running CSF QC, Harmony integration, and clustering ---\n")

set.seed(PCA_SEED)

csf_clustered <- run_qc_and_clustering(
  sce = csf_raw,
  min_genes_per_cell = MIN_GENES_PER_CELL,
  max_pct_mt = MAX_PCT_MT,
  resolution = CSF_RESOLUTION,
  n_pcs = N_PCS,
  clustering_algorithm = CLUSTER_ALGORITHM,
  clustering_seed = CLUSTER_SEED,
  umap_seed = UMAP_SEED
)

if (!inherits(csf_clustered, "SingleCellExperiment")) {
  stop("CSF clustering did not return a SingleCellExperiment.")
}

required_clustered_metadata <- c(
  "patient", "sample_id", "moca_group", "seurat_clusters"
)
missing_clustered_metadata <- setdiff(
  required_clustered_metadata,
  colnames(SummarizedExperiment::colData(csf_clustered))
)

if (length(missing_clustered_metadata) > 0L) {
  stop(
    "Clustered CSF object is missing metadata field(s): ",
    paste(missing_clustered_metadata, collapse = ", ")
  )
}

qc_counts_by_patient <- table(
  factor(
    as.character(csf_clustered$patient),
    levels = EXPECTED_ALL_PATIENTS
  )
)

if (any(qc_counts_by_patient == 0L)) {
  stop("One or more CSF donors have no cells after QC.")
}

cat("Post-QC ventricular-CSF cells by patient:\n")
print(qc_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 4 · IFN-independent classical-monocyte candidate selection
# -----------------------------------------------------------------------------

cat("\n--- Selecting CSF classical-monocyte candidates independently of IFN ---\n")

csf_annotated <- identify_monocytes(
  sce = csf_clustered,
  markers = CLASSICAL_MONOCYTE_SELECTION_PANEL,
  relative_threshold = MONOCYTE_RELATIVE_THRESHOLD,
  clusters_override = NULL
)

selected_clusters <- as.character(
  S4Vectors::metadata(csf_annotated)$monocyte_clusters
)
cluster_ids <- as.character(csf_annotated$seurat_clusters)
selected_mask <- as.character(csf_annotated$cell_type) == "monocyte"

if (length(selected_clusters) == 0L || sum(selected_mask) == 0L) {
  stop("No CSF classical-monocyte candidate cluster was selected.")
}

if (!"logcounts" %in% SummarizedExperiment::assayNames(csf_annotated)) {
  csf_annotated <- scuttle::logNormCounts(csf_annotated)
}

logcounts_csf <- SingleCellExperiment::logcounts(csf_annotated)
all_clusters <- sort(unique(cluster_ids))

score_cluster_panel <- function(genes) {
  available <- intersect(genes, rownames(logcounts_csf))
  if (length(available) == 0L) {
    return(setNames(rep(NA_real_, length(all_clusters)), all_clusters))
  }
  per_cell <- Matrix::colMeans(
    logcounts_csf[available, , drop = FALSE]
  )
  scores <- tapply(per_cell, cluster_ids, mean)
  as.numeric(scores[all_clusters])
}

cluster_diagnostics <- data.frame(
  cluster = all_clusters,
  n_cells = as.integer(table(cluster_ids)[all_clusters]),
  selected_as_classical_monocyte_candidate =
    all_clusters %in% selected_clusters,
  stringsAsFactors = FALSE
)

for (panel_name in names(DIAGNOSTIC_PANELS)) {
  cluster_diagnostics[[panel_name]] <- score_cluster_panel(
    DIAGNOSTIC_PANELS[[panel_name]]
  )
}

selection_scores <- as.numeric(
  S4Vectors::metadata(csf_annotated)$monocyte_cluster_scores
)
names(selection_scores) <- names(
  S4Vectors::metadata(csf_annotated)$monocyte_cluster_scores
)

if (!all(selected_clusters %in% names(selection_scores)) ||
    any(
      selection_scores[selected_clusters] <
        MONOCYTE_RELATIVE_THRESHOLD * max(selection_scores)
    )) {
  stop("Selected CSF clusters do not reproduce the frozen marker rule.")
}

csf_monocytes <- csf_annotated[, selected_mask, drop = FALSE]
csf_monocytes$csf_annotation <- rep(
  "IFN-independent_classical_monocyte_candidate",
  ncol(csf_monocytes)
)

candidate_counts_by_patient <- table(
  factor(
    as.character(csf_monocytes$patient),
    levels = EXPECTED_ALL_PATIENTS
  )
)

cat("Selected CSF cluster(s):", paste(selected_clusters, collapse = ", "), "\n")
cat("CSF classical-monocyte candidates by patient:\n")
print(candidate_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 5 · Apply the frozen blood-derived AUCell scores and threshold
# -----------------------------------------------------------------------------

cat(
  "\n--- Applying frozen blood-derived AUCell scoring to ",
  ncol(csf_monocytes), " CSF candidates ---\n",
  sep = ""
)

csf_counts <- SummarizedExperiment::assay(csf_monocytes, "counts")

required_scoring_genes <- unique(c(
  signature_up_exact,
  signature_down_exact,
  signature_up_independent,
  signature_down_independent,
  EXPECTED_ISG_SELECTOR
))
missing_scoring_genes <- setdiff(required_scoring_genes, rownames(csf_counts))

if (length(missing_scoring_genes) > 0L) {
  stop(
    "CSF counts are missing scoring gene(s): ",
    paste(missing_scoring_genes, collapse = ", ")
  )
}

scoring_sets <- list(
  transfer_up_exact = signature_up_exact,
  transfer_down_exact = signature_down_exact,
  transfer_up_selection_independent = signature_up_independent,
  transfer_down_selection_independent = signature_down_independent,
  isg_selector = EXPECTED_ISG_SELECTOR
)

set.seed(SCORING_SEED)
csf_rankings <- AUCell::AUCell_buildRankings(
  csf_counts,
  plotStats = FALSE,
  verbose = FALSE
)
csf_auc <- AUCell::AUCell_calcAUC(
  scoring_sets,
  csf_rankings,
  aucMaxRank = AUC_MAXRANK,
  verbose = FALSE
)
csf_auc_matrix <- AUCell::getAUC(csf_auc)

if (!all(names(scoring_sets) %in% rownames(csf_auc_matrix))) {
  stop("The CSF AUCell result is missing one or more frozen scores.")
}

csf_cell_ids <- colnames(csf_counts)
score_row <- function(name) {
  as.numeric(csf_auc_matrix[name, csf_cell_ids])
}

transfer_score_exact <-
  score_row("transfer_up_exact") -
  score_row("transfer_down_exact")
transfer_score_independent <-
  score_row("transfer_up_selection_independent") -
  score_row("transfer_down_selection_independent")
marker_ifn_score <- score_row("isg_selector")

if (any(!is.finite(c(
  transfer_score_exact,
  transfer_score_independent,
  marker_ifn_score
)))) {
  stop("CSF AUCell scoring produced missing or non-finite values.")
}

if (length(csf_cell_ids) != ncol(csf_monocytes) ||
    anyDuplicated(csf_cell_ids)) {
  stop("CSF cell identifiers are invalid after AUCell scoring.")
}

umap_coordinates <- SingleCellExperiment::reducedDim(
  csf_monocytes,
  "UMAP"
)

if (nrow(umap_coordinates) != ncol(csf_monocytes) ||
    ncol(umap_coordinates) != 2L) {
  stop("The selected CSF monocyte UMAP has unexpected dimensions.")
}

csf_cell_scores <- data.frame(
  cell_id = csf_cell_ids,
  patient = as.character(csf_monocytes$patient),
  moca_group = as.character(csf_monocytes$moca_group),
  csf_cluster = as.character(csf_monocytes$seurat_clusters),
  classical_selection_score = as.numeric(csf_monocytes$marker_score),
  transfer_score_exact = transfer_score_exact,
  transfer_score_selection_independent = transfer_score_independent,
  marker_ifn_score = marker_ifn_score,
  marker_ifn_state = ifelse(
    marker_ifn_score >= IFN_THRESHOLD,
    "marker_ifn_high",
    "marker_ifn_low"
  ),
  UMAP_1 = umap_coordinates[, 1],
  UMAP_2 = umap_coordinates[, 2],
  stringsAsFactors = FALSE
)

score_correlation <- stats::cor(
  csf_cell_scores$transfer_score_exact,
  csf_cell_scores$transfer_score_selection_independent,
  method = "spearman"
)

if (!is.finite(score_correlation)) {
  stop("Cannot calculate the CSF transfer-score sensitivity correlation.")
}


# -----------------------------------------------------------------------------
# Section 6 · Donor-level summaries and matched blood/CSF comparison
# -----------------------------------------------------------------------------

summarize_cells_by_donor <- function(
    cell_table,
    donor_column,
    transfer_column,
    marker_column,
    compartment
) {
  donors <- unique(as.character(cell_table[[donor_column]]))
  donors <- donors[order(as.numeric(donors))]

  do.call(rbind, lapply(donors, function(patient_id) {
    keep <- as.character(cell_table[[donor_column]]) == patient_id
    transfer <- as.numeric(cell_table[[transfer_column]][keep])
    marker <- as.numeric(cell_table[[marker_column]][keep])
    high <- marker >= IFN_THRESHOLD

    data.frame(
      compartment = compartment,
      patient = patient_id,
      n_cells = sum(keep),
      marker_ifn_high_n = sum(high),
      marker_ifn_high_pct = 100 * mean(high),
      median_transfer_score = stats::median(transfer),
      mean_transfer_score = mean(transfer),
      median_marker_ifn_score = stats::median(marker),
      mean_marker_ifn_score = mean(marker),
      stringsAsFactors = FALSE
    )
  }))
}

csf_donor_summary <- summarize_cells_by_donor(
  csf_cell_scores,
  donor_column = "patient",
  transfer_column = "transfer_score_exact",
  marker_column = "marker_ifn_score",
  compartment = "ventricular_CSF_classical_monocyte_candidates"
)

canonical_scores$marker_ifn_state <- ifelse(
  canonical_scores$A_isg >= IFN_THRESHOLD,
  "marker_ifn_high",
  "marker_ifn_low"
)

pbmc_donor_summary <- summarize_cells_by_donor(
  canonical_scores,
  donor_column = "donor",
  transfer_column = "A_score",
  marker_column = "A_isg",
  compartment = "canonical_PBMC_classical_monocytes"
)

paired_donor_summary <- merge(
  pbmc_donor_summary,
  csf_donor_summary,
  by = "patient",
  suffixes = c("_PBMC", "_CSF"),
  all = FALSE,
  sort = FALSE
)
paired_donor_summary <- paired_donor_summary[
  match(EXPECTED_PRIMARY_PATIENTS, paired_donor_summary$patient),
  ,
  drop = FALSE
]

if (nrow(paired_donor_summary) != EXPECTED_PRIMARY_PBMC_DONORS ||
    anyNA(paired_donor_summary$patient) ||
    !identical(
      as.character(paired_donor_summary$patient),
      EXPECTED_PRIMARY_PATIENTS
    )) {
  stop("The paired PBMC/CSF summary does not contain the primary nine donors.")
}

paired_donor_summary$transfer_score_difference_CSF_minus_PBMC <-
  paired_donor_summary$median_transfer_score_CSF -
  paired_donor_summary$median_transfer_score_PBMC
paired_donor_summary$marker_score_difference_CSF_minus_PBMC <-
  paired_donor_summary$median_marker_ifn_score_CSF -
  paired_donor_summary$median_marker_ifn_score_PBMC
paired_donor_summary$marker_high_pct_difference_CSF_minus_PBMC <-
  paired_donor_summary$marker_ifn_high_pct_CSF -
  paired_donor_summary$marker_ifn_high_pct_PBMC

paired_test <- function(csf_values, pbmc_values, metric) {
  test <- suppressWarnings(stats::wilcox.test(
    csf_values,
    pbmc_values,
    paired = TRUE,
    exact = FALSE,
    alternative = "two.sided"
  ))

  data.frame(
    metric = metric,
    donors_n = length(csf_values),
    PBMC_median = stats::median(pbmc_values),
    CSF_median = stats::median(csf_values),
    median_paired_difference_CSF_minus_PBMC = stats::median(
      csf_values - pbmc_values
    ),
    wilcoxon_V = unname(test$statistic),
    nominal_unadjusted_p = test$p.value,
    interpretation_scope = paste(
      "exploratory paired compartment comparison;",
      "not a disease-specificity or causal test"
    ),
    stringsAsFactors = FALSE
  )
}

paired_tests <- rbind(
  paired_test(
    paired_donor_summary$median_transfer_score_CSF,
    paired_donor_summary$median_transfer_score_PBMC,
    "median blood-derived transfer score"
  ),
  paired_test(
    paired_donor_summary$median_marker_ifn_score_CSF,
    paired_donor_summary$median_marker_ifn_score_PBMC,
    "median frozen IFN-selector score"
  ),
  paired_test(
    paired_donor_summary$marker_ifn_high_pct_CSF,
    paired_donor_summary$marker_ifn_high_pct_PBMC,
    "percentage above frozen IFN threshold"
  )
)


# -----------------------------------------------------------------------------
# Section 7 · Gene-level detection and neutral result summary
# -----------------------------------------------------------------------------

available_expected_isgs <- intersect(expected_isgs, rownames(csf_counts))
isg_counts <- csf_counts[available_expected_isgs, , drop = FALSE]
donor_vector <- as.character(csf_monocytes$patient)

isg_detection <- do.call(rbind, lapply(
  available_expected_isgs,
  function(gene) {
    gene_counts <- as.numeric(isg_counts[gene, ])
    donor_detected <- vapply(
      EXPECTED_ALL_PATIENTS,
      function(patient_id) {
        any(gene_counts[donor_vector == patient_id] > 0)
      },
      logical(1)
    )

    data.frame(
      gene = gene,
      expressing_cells_n = sum(gene_counts > 0),
      expressing_cells_pct = 100 * mean(gene_counts > 0),
      donors_with_detected_counts_n = sum(donor_detected),
      detected_in_all_ten_donors = all(donor_detected),
      total_raw_counts = sum(gene_counts),
      stringsAsFactors = FALSE
    )
  }
))

csf_high_n <- sum(csf_cell_scores$marker_ifn_state == "marker_ifn_high")
csf_high_pct <- 100 * csf_high_n / nrow(csf_cell_scores)
donors_with_high <- sum(csf_donor_summary$marker_ifn_high_n > 0L)
donors_with_at_least_5_high <- sum(
  csf_donor_summary$marker_ifn_high_n >= 5L
)
donors_with_at_least_10_high <- sum(
  csf_donor_summary$marker_ifn_high_n >= 10L
)

contextualization_summary <- data.frame(
  source_dataset = "GSE292141",
  CSF_source = "ventricular CSF collected during iNPH neurosurgical care",
  raw_CSF_cells_n = ncol(csf_raw),
  post_QC_CSF_cells_n = ncol(csf_clustered),
  selected_CSF_classical_monocyte_candidates_n = ncol(csf_monocytes),
  selected_CSF_clusters = paste(selected_clusters, collapse = " | "),
  CSF_donors_n = length(unique(csf_cell_scores$patient)),
  primary_matched_PBMC_CSF_donors_n = nrow(paired_donor_summary),
  frozen_IFN_threshold = IFN_THRESHOLD,
  CSF_marker_IFN_high_n = csf_high_n,
  CSF_marker_IFN_high_pct = csf_high_pct,
  CSF_donors_with_at_least_1_IFN_high_cell = donors_with_high,
  CSF_donors_with_at_least_5_IFN_high_cells = donors_with_at_least_5_high,
  CSF_donors_with_at_least_10_IFN_high_cells =
    donors_with_at_least_10_high,
  median_CSF_transfer_score = stats::median(
    csf_cell_scores$transfer_score_exact
  ),
  transfer_score_independence_spearman = score_correlation,
  expected_ISGs_detected_in_at_least_1_CSF_cell = sum(
    isg_detection$expressing_cells_n > 0L
  ),
  expected_ISGs_detected_in_all_ten_CSF_donors = sum(
    isg_detection$detected_in_all_ten_donors
  ),
  biological_result_status =
    "REPORT_EXACT_RESULTS_WITHOUT_PREDECLARED_DIRECTIONAL_PASS_FAIL",
  interpretation = paste(
    "Tests cross-compartment detectability of the blood-derived",
    "IFN-response-high program; does not establish iNPH specificity,",
    "CNS causality, cognition association, or treatment efficacy"
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# Section 8 · Input provenance and outputs
# -----------------------------------------------------------------------------

raw_csf_prefixes <- c(
  "GSM8851048_1C",
  "GSM8851050_2C",
  "GSM8851052_3C",
  "GSM8851054_4C",
  "GSM8851056_5C",
  "GSM8851058_6C",
  "GSM8851060_101822C",
  "GSM8851062_102522C2",
  "GSM8851064_102522C3",
  "GSM8851066_110122C"
)
raw_suffixes <- c(
  "_barcodes.tsv.gz",
  "_features.tsv.gz",
  "_matrix.mtx.gz"
)
raw_csf_files <- unlist(lapply(
  raw_csf_prefixes,
  function(prefix) file.path(DUY_EXTRACTED, paste0(prefix, raw_suffixes))
), use.names = FALSE)

if (!all(file.exists(raw_csf_files)) ||
    any(is.na(file.info(raw_csf_files)$size)) ||
    any(file.info(raw_csf_files)$size <= 0)) {
  stop("One or more expected raw CSF 10x files are missing or empty.")
}

raw_input_inventory <- data.frame(
  file = raw_csf_files,
  basename = basename(raw_csf_files),
  size_bytes = as.numeric(file.info(raw_csf_files)$size),
  modified = as.character(file.info(raw_csf_files)$mtime),
  stringsAsFactors = FALSE
)

# Hash only immediate code and processed analytical inputs. The 30 public raw
# files are inventoried by exact GSM filename and size to avoid an unnecessary
# long-running transitive checksum pass.
immediate_inputs <- c(
  DATA_HELPER,
  CLUSTER_HELPER,
  SIGNATURE_PATH,
  ISG_REFERENCE_PATH,
  PHASE2_PATH
)
own_script <- here::here(
  "analysis",
  "05B_csf_monocyte_ifn_contextualization.R"
)
if (file.exists(own_script)) {
  immediate_inputs <- c(immediate_inputs, own_script)
}
immediate_md5_before <- unname(tools::md5sum(immediate_inputs))

if (anyNA(immediate_md5_before)) {
  stop("An immediate CSF-contextualization checksum could not be calculated.")
}

reporting_notes <- c(
  BUILD,
  "REPORTING GUARDRAILS",
  "1. Describe the transferred program as a blood-derived IFN-response-high monocyte program observed in iNPH.",
  "2. Describe CSF cells as independently selected classical-monocyte candidates from ventricular CSF.",
  "3. Report exact cell and donor counts; do not convert the result into a disease-specificity claim.",
  "4. Treat paired PBMC-versus-CSF tests as exploratory compartment comparisons with nominal unadjusted P values.",
  "5. Do not imply that CSF enrichment proves CNS origin, brain-to-blood signaling, cognition effects, or therapeutic efficacy.",
  "6. Patient 6 is descriptive because MoCA is unknown and is not part of the canonical nine-donor PBMC comparison.",
  "7. No High-versus-Low MoCA test is performed in this analysis.",
  paste("Frozen AUCell threshold:", format(IFN_THRESHOLD, digits = 10)),
  paste("Frozen aucMaxRank:", AUC_MAXRANK),
  paste("Selection-panel overlap with transfer UP:",
        paste(selection_overlap_up, collapse = " | ")),
  paste("Selection-panel overlap with transfer DOWN:",
        paste(selection_overlap_down, collapse = " | ")),
  paste("Source paper: Duy et al., PNAS 2025, PMID 40324076,"),
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12087963/"
)

result <- list(
  analysis_id = "revision_05B_csf_monocyte_ifn_contextualization",
  build_id = BUILD,
  status = "PASS_ANALYSIS_COMPLETE_INTERPRET_FROM_EXPORTED_RESULTS",
  configuration = list(
    source_dataset = "GSE292141",
    tissue = "ventricular CSF",
    patients = EXPECTED_ALL_PATIENTS,
    primary_matched_patients = EXPECTED_PRIMARY_PATIENTS,
    min_genes_per_cell = MIN_GENES_PER_CELL,
    max_pct_mt = MAX_PCT_MT,
    resolution = CSF_RESOLUTION,
    n_pcs = N_PCS,
    clustering_algorithm = CLUSTER_ALGORITHM,
    clustering_seed = CLUSTER_SEED,
    pca_seed = PCA_SEED,
    umap_seed = UMAP_SEED,
    selection_panel = CLASSICAL_MONOCYTE_SELECTION_PANEL,
    selection_relative_threshold = MONOCYTE_RELATIVE_THRESHOLD,
    selection_uses_IFN_genes = FALSE,
    auc_maxrank = AUC_MAXRANK,
    scoring_seed = SCORING_SEED,
    IFN_selector = EXPECTED_ISG_SELECTOR,
    frozen_IFN_threshold = IFN_THRESHOLD,
    transfer_up_exact = signature_up_exact,
    transfer_down_exact = signature_down_exact,
    transfer_up_selection_independent = signature_up_independent,
    transfer_down_selection_independent = signature_down_independent,
    selection_overlap_up = selection_overlap_up,
    selection_overlap_down = selection_overlap_down
  ),
  contextualization_summary = contextualization_summary,
  csf_cell_scores = csf_cell_scores,
  csf_donor_summary = csf_donor_summary,
  canonical_pbmc_donor_summary = pbmc_donor_summary,
  paired_donor_summary = paired_donor_summary,
  paired_tests = paired_tests,
  cluster_diagnostics = cluster_diagnostics,
  isg_detection = isg_detection,
  reporting_notes = reporting_notes,
  input_provenance = list(
    immediate_files = data.frame(
      file = immediate_inputs,
      md5 = immediate_md5_before,
      stringsAsFactors = FALSE
    ),
    raw_CSF_files = raw_input_inventory,
    signature_md5 = unname(tools::md5sum(SIGNATURE_PATH)),
    phase2_md5 = unname(tools::md5sum(PHASE2_PATH))
  ),
  session_info = utils::sessionInfo()
)

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REVISION_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(csf_monocytes, OUTPUT_MONOCYTES_RDS, compress = TRUE)
saveRDS(result, OUTPUT_RESULT_RDS, compress = TRUE)
utils::write.csv(
  csf_cell_scores,
  OUTPUT_CELL_SCORES_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  csf_donor_summary,
  OUTPUT_DONOR_SUMMARY_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  paired_donor_summary,
  OUTPUT_PAIRED_SUMMARY_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  paired_tests,
  OUTPUT_PAIRED_TESTS_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  cluster_diagnostics,
  OUTPUT_CLUSTER_DIAGNOSTICS_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  isg_detection,
  OUTPUT_ISG_DETECTION_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  contextualization_summary,
  OUTPUT_SUMMARY_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  raw_input_inventory,
  OUTPUT_RAW_INVENTORY_CSV,
  row.names = FALSE,
  na = ""
)
writeLines(reporting_notes, OUTPUT_REPORTING_NOTES)

if (!all(file.exists(OUTPUT_PATHS)) ||
    any(is.na(file.info(OUTPUT_PATHS)$size)) ||
    any(file.info(OUTPUT_PATHS)$size <= 0)) {
  stop("One or more CSF-contextualization outputs are missing or empty.")
}

result_check <- readRDS(OUTPUT_RESULT_RDS)
monocyte_check <- readRDS(OUTPUT_MONOCYTES_RDS)

if (!identical(result_check, result) ||
    !inherits(monocyte_check, "SingleCellExperiment") ||
    !identical(colnames(monocyte_check), colnames(csf_monocytes)) ||
    !identical(
      as.character(result_check$csf_cell_scores$cell_id),
      colnames(monocyte_check)
    )) {
  stop("A saved CSF-contextualization output failed read-back validation.")
}

if (!identical(
  immediate_md5_before,
  unname(tools::md5sum(immediate_inputs))
)) {
  stop("An immediate analytical input changed during CSF contextualization.")
}

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 05B COMPLETE\n")
cat("============================================================\n")
cat("  Raw ventricular-CSF cells       :", ncol(csf_raw), "\n")
cat("  Post-QC ventricular-CSF cells   :", ncol(csf_clustered), "\n")
cat("  Selected classical candidates   :", ncol(csf_monocytes), "\n")
cat("  Selected CSF cluster(s)         :", paste(selected_clusters, collapse = ", "), "\n")
cat("  CSF donors represented          :", length(unique(csf_cell_scores$patient)), "\n")
cat("  Frozen-threshold IFN-high cells :", csf_high_n, "/", nrow(csf_cell_scores),
    sprintf("(%.2f%%)", csf_high_pct), "\n")
cat("  Donors with >=1 IFN-high cell   :", donors_with_high, "/ 10\n")
cat("  Donors with >=10 IFN-high cells :", donors_with_at_least_10_high, "/ 10\n")
cat("  Primary paired PBMC/CSF donors  :", nrow(paired_donor_summary), "/ 9\n")
cat("  Transfer-score sensitivity rho  :", round(score_correlation, 4), "\n")
cat("  Expected ISGs detected          :",
    sum(isg_detection$expressing_cells_n > 0L), "/ 35\n")
cat("\nExploratory paired compartment tests:\n")
print(paired_tests, row.names = FALSE, digits = 4)
cat("\nPer-donor CSF summary:\n")
print(csf_donor_summary, row.names = FALSE, digits = 4)
cat("\nOutputs saved:\n")
cat(" ", OUTPUT_RESULT_RDS, "\n")
cat(" ", OUTPUT_MONOCYTES_RDS, "\n")
cat(" ", REVISION_RESULTS_DIR, "\n")
cat("============================================================\n")
cat(
  "Interpret exact results as cross-compartment contextualization only; ",
  "do not infer disease specificity, cognition effects, CNS causality, ",
  "or therapeutic efficacy.\n",
  sep = ""
)
