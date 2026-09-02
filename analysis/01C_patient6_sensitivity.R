# =============================================================================
# Revision analysis 1C — Patient 6 ten-patient PBMC sensitivity
# Script: 01C_patient6_sensitivity.R
# =============================================================================
#
# Project
# -------
# Systems biology prioritizes drug candidates for interferon-activated
# monocytes in idiopathic normal pressure hydrocephalus
#
# Purpose
# -------
# Reconstruct a uniform marker-enriched classical-monocyte candidate population
# from all ten GSE292141 PBMC samples, including Patient 6, whose MoCA value is
# unknown. The same PBMC QC, Harmony integration, Louvain clustering, and fixed
# marker-based annotation rule are applied to every patient.
#
# The submitted primary analysis remains the locked nine-patient canonical
# population. This script is a sensitivity analysis only. It validates the
# uniformly reconstructed population against the locked nine-patient manifest
# before declaring it eligible for downstream marker-defined IFN-state work.
#
# Revision-stage choices (frozen before execution)
# ------------------------------------------------
# - Include all ten PBMC samples: drop_unknown_moca = FALSE.
# - QC: nFeature_RNA >= 200 and percent_mt <= 20.
# - PBMC representation: 2,000 variable features, 30 PCs, Harmony by patient.
# - Clustering: Louvain algorithm 1, resolution 0.5, random.seed = 0.
# - Classical-monocyte annotation panel: CD14, S100A8, S100A9, S100A12,
#   VCAN, FCN1 (the same classical panel used in Script 01).
# - Marker-enriched clusters: mean panel score >= 0.75 times the highest
#   cluster score, using the existing identify_monocytes() helper.
# - IFN score: mean normalized expression of IFITM3, IFITM2, ISG15, MX1,
#   OAS1, IRF7, STAT1; no IFN-high threshold is applied in this script.
#
# Operational eligibility criteria for later sensitivity analyses
# ---------------------------------------------------------------
# - At least 80% of locked canonical cells are recovered among Patients 1–5
#   and 7–10.
# - Jaccard overlap with the locked canonical population is at least 0.50.
# - All ten patients contribute at least 30 candidate cells.
#
# These criteria are operational revision-stage quality checks, not externally
# validated biological thresholds. Failure does not alter the primary analysis;
# it marks this Patient 6 branch as unsuitable for downstream signature work.
#
# Dependency boundary
# -------------------
# This script does not define marker-based IFN-high/IFN-low states, perform
# differential expression, or run LINCS. The shared healthy-control-anchored
# IFN threshold is generated downstream from Script 02.
#
# Inputs
# -----
# data/raw/duy/extracted/<all GSE292141 10x triplets>
# data/canonical/duy_classical_monocyte_cell_manifest.tsv
# R/data_acquisition_and_loading.R
# R/pbmc_clustering_and_annotation.R
#
# Outputs
# -------
# data/processed/revision_patient6_uniform_monocytes.rds
# data/processed/revision_patient6_sensitivity.rds
# results/revision/patient6_reconstruction_summary.csv
# results/revision/patient6_canonical_overlap_by_donor.csv
# results/revision/patient6_cluster_diagnostics.csv
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment and frozen configuration
# -----------------------------------------------------------------------------

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
CANONICAL_DIR <- here::here("data", "canonical")
REVISION_RESULTS_DIR <- here::here("results", "revision")

DUY_EXTRACTED <- file.path(RAW_DIR, "duy", "extracted")
MANIFEST_PATH <- file.path(
  CANONICAL_DIR,
  "duy_classical_monocyte_cell_manifest.tsv"
)
DATA_HELPER <- here::here("R", "data_acquisition_and_loading.R")
PBMC_HELPER <- here::here("R", "pbmc_clustering_and_annotation.R")

OUTPUT_MONOCYTES_RDS <- file.path(
  PROCESSED_DIR,
  "revision_patient6_uniform_monocytes.rds"
)
OUTPUT_RESULT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_patient6_sensitivity.rds"
)
OUTPUT_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "patient6_reconstruction_summary.csv"
)
OUTPUT_DONOR_OVERLAP_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "patient6_canonical_overlap_by_donor.csv"
)
OUTPUT_CLUSTER_DIAGNOSTICS_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "patient6_cluster_diagnostics.csv"
)

OUTPUT_PATHS <- c(
  OUTPUT_MONOCYTES_RDS,
  OUTPUT_RESULT_RDS,
  OUTPUT_SUMMARY_CSV,
  OUTPUT_DONOR_OVERLAP_CSV,
  OUTPUT_CLUSTER_DIAGNOSTICS_CSV
)

EXPECTED_GENES <- 36601L
EXPECTED_PATIENT_IDS <- as.character(1:10)
EXPECTED_CANONICAL_CELLS <- 22197L
EXPECTED_CANONICAL_PATIENT_IDS <- c(
  "1", "2", "3", "4", "5", "7", "8", "9", "10"
)

MIN_GENES_PER_CELL <- 200L
MAX_PCT_MT <- 20
PBMC_RESOLUTION <- 0.5
N_PCS <- 30L
CLUSTER_ALGORITHM <- 1L
CLUSTER_SEED <- 0L
UMAP_SEED <- 42L
PCA_SEED <- 42L
MONOCYTE_RELATIVE_THRESHOLD <- 0.75
MIN_CELLS_PER_PATIENT <- 30L
MIN_CANONICAL_RECALL <- 0.80
MIN_CANONICAL_JACCARD <- 0.50

CLASSICAL_MONOCYTE_PANEL <- c(
  "CD14", "S100A8", "S100A9", "S100A12", "VCAN", "FCN1"
)

IFN_PANEL <- c(
  "IFITM3", "IFITM2", "ISG15", "MX1", "OAS1", "IRF7", "STAT1"
)

required_files <- c(MANIFEST_PATH, DATA_HELPER, PBMC_HELPER)
required_info <- file.info(required_files)

bad_required <- required_files[
  !file.exists(required_files) |
    is.na(required_info$size) |
    required_info$size <= 0
]

if (length(bad_required) > 0L) {
  stop(
    "Required Patient 6 input file(s) missing or empty:\n  ",
    paste(bad_required, collapse = "\n  ")
  )
}

if (!dir.exists(DUY_EXTRACTED)) {
  stop("Duy extracted-data directory is missing: ", DUY_EXTRACTED)
}

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_PATIENT6", unset = "0"),
  "1"
)

existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "Patient 6 sensitivity output(s) already exist:\n  ",
    paste(existing_outputs, collapse = "\n  "),
    "\nRefusing to overwrite them automatically. If replacement is ",
    "intentional, set INPH_OVERWRITE_PATIENT6=1 for that run."
  )
}

source(DATA_HELPER)
source(PBMC_HELPER)


# -----------------------------------------------------------------------------
# Section 1 · Load all ten PBMC samples and restore Patient 6 metadata
# -----------------------------------------------------------------------------

cat("\n--- Loading all ten GSE292141 PBMC samples ---\n")

pbmc_raw <- load_duy_samples(
  extracted_dir = DUY_EXTRACTED,
  tissue_filter = "PBMC",
  drop_unknown_moca = FALSE
)

if (!inherits(pbmc_raw, "SingleCellExperiment")) {
  stop("load_duy_samples() did not return a SingleCellExperiment.")
}

if (nrow(pbmc_raw) != EXPECTED_GENES) {
  stop(
    "Expected ", EXPECTED_GENES,
    " raw features; observed ", nrow(pbmc_raw), "."
  )
}

if (anyDuplicated(colnames(pbmc_raw))) {
  stop("Ten-patient PBMC reconstruction contains duplicate cell identifiers.")
}

observed_patient_ids <- unique(as.character(pbmc_raw$patient))

if (!setequal(observed_patient_ids, EXPECTED_PATIENT_IDS)) {
  stop(
    "Expected PBMC samples for Patients 1–10; observed: ",
    paste(observed_patient_ids, collapse = ", ")
  )
}

if (!all(as.character(pbmc_raw$tissue) == "PBMC")) {
  stop("Ten-patient reconstruction contains non-PBMC cells.")
}

# load_duy_samples() currently factors MoCA as High/Low even when unknown
# samples are retained, which converts Patient 6's documented Unknown value to
# NA. Restore that known sample-level label explicitly for this sensitivity run.
moca_group <- as.character(pbmc_raw$moca_group)
moca_group[as.character(pbmc_raw$patient) == "6"] <- "Unknown"
pbmc_raw$moca_group <- factor(
  moca_group,
  levels = c("High", "Low", "Unknown")
)

if (anyNA(pbmc_raw$moca_group)) {
  stop("Missing MoCA-group metadata remain after restoring Patient 6.")
}

if (!all(
  as.character(pbmc_raw$moca_group[as.character(pbmc_raw$patient) == "6"]) ==
    "Unknown"
)) {
  stop("Patient 6 MoCA metadata were not restored as Unknown.")
}

# Normalize loader-generated temporary paths before downstream conversion.
if ("Sample" %in% colnames(SummarizedExperiment::colData(pbmc_raw))) {
  pbmc_raw$Sample <- as.character(pbmc_raw$sample_id)
}

pbmc_raw_metadata <- S4Vectors::metadata(pbmc_raw)
pbmc_raw_metadata$Samples <- unique(as.character(pbmc_raw$sample_id))
S4Vectors::metadata(pbmc_raw) <- pbmc_raw_metadata

pbmc_raw <- fix_gene_symbols(pbmc_raw)

if (nrow(pbmc_raw) != EXPECTED_GENES || anyDuplicated(rownames(pbmc_raw))) {
  stop("Gene-symbol harmonization produced an unexpected feature structure.")
}

raw_counts_by_patient <- table(as.character(pbmc_raw$patient))

cat(
  "Ten-patient PBMC input: ",
  nrow(pbmc_raw), " genes x ", ncol(pbmc_raw), " cells\n",
  sep = ""
)
print(raw_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 2 · Apply the uniform PBMC reconstruction
# -----------------------------------------------------------------------------

set.seed(PCA_SEED)

pbmc_clustered <- run_qc_and_clustering(
  sce = pbmc_raw,
  min_genes_per_cell = MIN_GENES_PER_CELL,
  max_pct_mt = MAX_PCT_MT,
  resolution = PBMC_RESOLUTION,
  n_pcs = N_PCS,
  clustering_algorithm = CLUSTER_ALGORITHM,
  clustering_seed = CLUSTER_SEED,
  umap_seed = UMAP_SEED
)

if (!inherits(pbmc_clustered, "SingleCellExperiment")) {
  stop("run_qc_and_clustering() did not return a SingleCellExperiment.")
}

required_clustered_metadata <- c(
  "patient",
  "sample_id",
  "moca_group",
  "seurat_clusters"
)

missing_clustered_metadata <- setdiff(
  required_clustered_metadata,
  colnames(SummarizedExperiment::colData(pbmc_clustered))
)

if (length(missing_clustered_metadata) > 0L) {
  stop(
    "Uniform PBMC object is missing required metadata field(s): ",
    paste(missing_clustered_metadata, collapse = ", ")
  )
}

qc_counts_by_patient <- table(as.character(pbmc_clustered$patient))

if (!setequal(names(qc_counts_by_patient), EXPECTED_PATIENT_IDS)) {
  stop("One or more patients are absent after uniform PBMC QC.")
}

cat("\nPost-QC PBMC cells by patient:\n")
print(qc_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 3 · Uniform marker-enriched classical-monocyte annotation
# -----------------------------------------------------------------------------

pbmc_annotated <- identify_monocytes(
  sce = pbmc_clustered,
  markers = CLASSICAL_MONOCYTE_PANEL,
  relative_threshold = MONOCYTE_RELATIVE_THRESHOLD,
  clusters_override = NULL
)

selected_clusters <- as.character(
  S4Vectors::metadata(pbmc_annotated)$monocyte_clusters
)
cluster_marker_scores <- S4Vectors::metadata(
  pbmc_annotated
)$monocyte_cluster_scores

if (length(selected_clusters) == 0L) {
  stop("Uniform annotation did not select any monocyte-enriched cluster.")
}

cluster_ids <- as.character(pbmc_annotated$seurat_clusters)
candidate_mask <- as.character(pbmc_annotated$cell_type) == "monocyte"

if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(pbmc_annotated)) {
  stop("Uniform PBMC object is missing the expected UMAP reduction.")
}

umap_coordinates <- SingleCellExperiment::reducedDim(
  pbmc_annotated,
  "UMAP"
)

if (nrow(umap_coordinates) != ncol(pbmc_annotated) ||
    ncol(umap_coordinates) != 2L) {
  stop("Uniform PBMC UMAP dimensions are unexpected.")
}

pbmc_cell_diagnostics <- data.frame(
  cell_id = colnames(pbmc_annotated),
  patient = as.character(pbmc_annotated$patient),
  moca_group = as.character(pbmc_annotated$moca_group),
  pbmc_cluster = cluster_ids,
  classical_marker_score = as.numeric(pbmc_annotated$marker_score),
  selected_as_classical_monocyte_candidate = candidate_mask,
  UMAP_1 = umap_coordinates[, 1],
  UMAP_2 = umap_coordinates[, 2],
  stringsAsFactors = FALSE
)

cluster_diagnostics <- data.frame(
  cluster = names(cluster_marker_scores),
  classical_marker_score = as.numeric(cluster_marker_scores),
  n_cells = as.integer(
    table(cluster_ids)[names(cluster_marker_scores)]
  ),
  selected_as_classical_monocyte_candidate =
    names(cluster_marker_scores) %in% selected_clusters,
  stringsAsFactors = FALSE
)

candidate_monocytes <- pbmc_annotated[, candidate_mask, drop = FALSE]
candidate_monocytes$uniform_annotation <- rep(
  "classical_monocyte_candidate",
  ncol(candidate_monocytes)
)

candidate_counts_by_patient <- table(
  factor(
    as.character(candidate_monocytes$patient),
    levels = EXPECTED_PATIENT_IDS
  )
)

cat("\nUniformly selected candidate cells by patient:\n")
print(candidate_counts_by_patient)


# -----------------------------------------------------------------------------
# Section 4 · Validate the nine-patient reconstruction against the manifest
# -----------------------------------------------------------------------------

manifest <- utils::read.delim(
  MANIFEST_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_manifest_columns <- c("cell_id", "patient", "tissue")
missing_manifest_columns <- setdiff(
  required_manifest_columns,
  colnames(manifest)
)

if (length(missing_manifest_columns) > 0L) {
  stop(
    "Canonical manifest is missing required column(s): ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}

if (nrow(manifest) != EXPECTED_CANONICAL_CELLS) {
  stop(
    "Expected ", EXPECTED_CANONICAL_CELLS,
    " canonical manifest cells; observed ", nrow(manifest), "."
  )
}

if (!setequal(
  unique(as.character(manifest$patient)),
  EXPECTED_CANONICAL_PATIENT_IDS
)) {
  stop("Canonical manifest patient composition is unexpected.")
}

canonical_ids <- manifest$cell_id
candidate_ids <- colnames(candidate_monocytes)
candidate_ids_9 <- candidate_ids[
  as.character(candidate_monocytes$patient) != "6"
]

canonical_overlap_ids <- intersect(canonical_ids, candidate_ids_9)
canonical_union_ids <- union(canonical_ids, candidate_ids_9)

overall_canonical_recall <- length(canonical_overlap_ids) /
  length(canonical_ids)
overall_candidate_precision <- if (length(candidate_ids_9) > 0L) {
  length(canonical_overlap_ids) / length(candidate_ids_9)
} else {
  NA_real_
}
overall_canonical_jaccard <- length(canonical_overlap_ids) /
  length(canonical_union_ids)

donor_overlap <- dplyr::bind_rows(lapply(
  EXPECTED_CANONICAL_PATIENT_IDS,
  function(patient_id) {
    canonical_patient_ids <- manifest$cell_id[
      as.character(manifest$patient) == patient_id
    ]
    candidate_patient_ids <- candidate_ids[
      as.character(candidate_monocytes$patient) == patient_id
    ]
    overlap_n <- length(intersect(
      canonical_patient_ids,
      candidate_patient_ids
    ))
    union_n <- length(union(
      canonical_patient_ids,
      candidate_patient_ids
    ))

    data.frame(
      patient = patient_id,
      canonical_n = length(canonical_patient_ids),
      candidate_n = length(candidate_patient_ids),
      overlap_n = overlap_n,
      union_n = union_n,
      canonical_recall = overlap_n / length(canonical_patient_ids),
      candidate_precision = if (length(candidate_patient_ids) > 0L) {
        overlap_n / length(candidate_patient_ids)
      } else {
        NA_real_
      },
      jaccard = overlap_n / union_n,
      stringsAsFactors = FALSE
    )
  }
))

all_patients_meet_minimum <- all(
  as.integer(candidate_counts_by_patient) >= MIN_CELLS_PER_PATIENT
)

eligible_for_downstream <-
  overall_canonical_recall >= MIN_CANONICAL_RECALL &&
  overall_canonical_jaccard >= MIN_CANONICAL_JACCARD &&
  all_patients_meet_minimum


# -----------------------------------------------------------------------------
# Section 5 · Calculate the seven-gene IFN score without defining IFN-high
# -----------------------------------------------------------------------------

if (!"logcounts" %in% SummarizedExperiment::assayNames(candidate_monocytes)) {
  stop("Uniform candidate-monocyte object is missing normalized logcounts.")
}

missing_ifn_genes <- setdiff(IFN_PANEL, rownames(candidate_monocytes))

if (length(missing_ifn_genes) > 0L) {
  stop(
    "Seven-gene IFN panel is incomplete in the candidate-monocyte object: ",
    paste(missing_ifn_genes, collapse = ", ")
  )
}

candidate_monocytes$interferon_panel_score <- Matrix::colMeans(
  SingleCellExperiment::logcounts(candidate_monocytes)[
    IFN_PANEL,
    ,
    drop = FALSE
  ]
)

patient_summary <- data.frame(
  patient = EXPECTED_PATIENT_IDS,
  moca_group = vapply(
    EXPECTED_PATIENT_IDS,
    function(patient_id) {
      unique_value <- unique(as.character(
        pbmc_clustered$moca_group[
          as.character(pbmc_clustered$patient) == patient_id
        ]
      ))
      if (length(unique_value) != 1L) {
        stop("Non-unique MoCA group for Patient ", patient_id, ".")
      }
      unique_value
    },
    character(1)
  ),
  raw_pbmc_n = as.integer(raw_counts_by_patient[EXPECTED_PATIENT_IDS]),
  post_qc_pbmc_n = as.integer(qc_counts_by_patient[EXPECTED_PATIENT_IDS]),
  candidate_monocyte_n = as.integer(
    candidate_counts_by_patient[EXPECTED_PATIENT_IDS]
  ),
  candidate_pct_post_qc_pbmc = 100 *
    as.integer(candidate_counts_by_patient[EXPECTED_PATIENT_IDS]) /
    as.integer(qc_counts_by_patient[EXPECTED_PATIENT_IDS]),
  mean_interferon_panel_score = vapply(
    EXPECTED_PATIENT_IDS,
    function(patient_id) {
      patient_scores <- candidate_monocytes$interferon_panel_score[
        as.character(candidate_monocytes$patient) == patient_id
      ]
      if (length(patient_scores) > 0L) mean(patient_scores) else NA_real_
    },
    numeric(1)
  ),
  median_interferon_panel_score = vapply(
    EXPECTED_PATIENT_IDS,
    function(patient_id) {
      patient_scores <- candidate_monocytes$interferon_panel_score[
        as.character(candidate_monocytes$patient) == patient_id
      ]
      if (length(patient_scores) > 0L) {
        stats::median(patient_scores)
      } else {
        NA_real_
      }
    },
    numeric(1)
  ),
  meets_30_cell_minimum = as.integer(
    candidate_counts_by_patient[EXPECTED_PATIENT_IDS]
  ) >= MIN_CELLS_PER_PATIENT,
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# Section 6 · Save results before reporting eligibility
# -----------------------------------------------------------------------------

result <- list(
  analysis_id = "revision_01C_patient6_sensitivity",
  configuration = list(
    include_patient6 = TRUE,
    qc = c(
      paste0("nFeature_RNA >= ", MIN_GENES_PER_CELL),
      paste0("percent_mt <= ", MAX_PCT_MT)
    ),
    pbmc_resolution = PBMC_RESOLUTION,
    n_pcs = N_PCS,
    clustering_algorithm = CLUSTER_ALGORITHM,
    clustering_seed = CLUSTER_SEED,
    pca_seed = PCA_SEED,
    umap_seed = UMAP_SEED,
    classical_monocyte_panel = CLASSICAL_MONOCYTE_PANEL,
    monocyte_relative_threshold = MONOCYTE_RELATIVE_THRESHOLD,
    ifn_panel = IFN_PANEL,
    minimum_cells_per_patient = MIN_CELLS_PER_PATIENT,
    minimum_canonical_recall = MIN_CANONICAL_RECALL,
    minimum_canonical_jaccard = MIN_CANONICAL_JACCARD,
    annotation_scope = paste(
      "uniform marker-enriched classical-monocyte candidate population;",
      "sensitivity analysis only"
    )
  ),
  selected_pbmc_clusters = selected_clusters,
  cluster_diagnostics = cluster_diagnostics,
  pbmc_cell_diagnostics = pbmc_cell_diagnostics,
  patient_summary = patient_summary,
  donor_overlap = donor_overlap,
  overall_validation = list(
    canonical_n = length(canonical_ids),
    candidate_n_9_patients = length(candidate_ids_9),
    overlap_n = length(canonical_overlap_ids),
    union_n = length(canonical_union_ids),
    canonical_recall = overall_canonical_recall,
    candidate_precision = overall_candidate_precision,
    jaccard = overall_canonical_jaccard,
    all_ten_patients_meet_30_cell_minimum = all_patients_meet_minimum,
    eligible_for_downstream = eligible_for_downstream
  ),
  input_provenance = list(
    raw_source = "GSE292141",
    canonical_manifest_path =
      "data/canonical/duy_classical_monocyte_cell_manifest.tsv",
    canonical_manifest_md5 = unname(tools::md5sum(MANIFEST_PATH)),
    data_helper_path = "R/data_acquisition_and_loading.R",
    data_helper_md5 = unname(tools::md5sum(DATA_HELPER)),
    pbmc_helper_path = "R/pbmc_clustering_and_annotation.R",
    pbmc_helper_md5 = unname(tools::md5sum(PBMC_HELPER))
  )
)

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REVISION_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(candidate_monocytes, OUTPUT_MONOCYTES_RDS)
saveRDS(result, OUTPUT_RESULT_RDS)

utils::write.csv(
  patient_summary,
  OUTPUT_SUMMARY_CSV,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  donor_overlap,
  OUTPUT_DONOR_OVERLAP_CSV,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  cluster_diagnostics,
  OUTPUT_CLUSTER_DIAGNOSTICS_CSV,
  row.names = FALSE,
  na = ""
)

if (!all(file.exists(OUTPUT_PATHS))) {
  stop("One or more Patient 6 sensitivity outputs were not created.")
}

if (any(is.na(file.info(OUTPUT_PATHS)$size)) ||
    any(file.info(OUTPUT_PATHS)$size <= 0)) {
  stop("One or more Patient 6 sensitivity outputs are empty.")
}

candidate_check <- readRDS(OUTPUT_MONOCYTES_RDS)
result_check <- readRDS(OUTPUT_RESULT_RDS)

if (!inherits(candidate_check, "SingleCellExperiment")) {
  stop("Saved candidate-monocyte object failed class validation.")
}

if (!identical(colnames(candidate_check), colnames(candidate_monocytes))) {
  stop("Saved candidate-monocyte cell order failed read-back validation.")
}

if (!identical(
  result_check$overall_validation$eligible_for_downstream,
  eligible_for_downstream
)) {
  stop("Saved Patient 6 eligibility flag failed read-back validation.")
}


# -----------------------------------------------------------------------------
# Section 7 · Closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 1C COMPLETE\n")
cat("============================================================\n")
cat("  Post-QC PBMC cells                : ", ncol(pbmc_clustered), "\n", sep = "")
cat("  Selected PBMC cluster(s)          : ",
    paste(selected_clusters, collapse = ", "), "\n", sep = "")
cat("  Candidate monocytes, all patients : ",
    ncol(candidate_monocytes), "\n", sep = "")
cat("  Patient 6 candidate cells         : ",
    candidate_counts_by_patient[["6"]], "\n", sep = "")
cat("  Canonical recovery, nine patients : ",
    sprintf("%.1f%%", 100 * overall_canonical_recall), "\n", sep = "")
cat("  Candidate precision               : ",
    sprintf("%.1f%%", 100 * overall_candidate_precision), "\n", sep = "")
cat("  Canonical Jaccard                  : ",
    sprintf("%.3f", overall_canonical_jaccard), "\n", sep = "")
cat("  All patients >=30 candidate cells : ",
    all_patients_meet_minimum, "\n", sep = "")
cat("  Eligible for downstream analysis  : ",
    eligible_for_downstream, "\n", sep = "")
cat("\nPer-patient summary:\n")
print(patient_summary, row.names = FALSE)
cat("\nCanonical overlap by donor:\n")
print(donor_overlap, row.names = FALSE)
cat("\nOutputs saved:\n")
cat("  ", OUTPUT_MONOCYTES_RDS, "\n", sep = "")
cat("  ", OUTPUT_RESULT_RDS, "\n", sep = "")
cat("  ", OUTPUT_SUMMARY_CSV, "\n", sep = "")
cat("  ", OUTPUT_DONOR_OVERLAP_CSV, "\n", sep = "")
cat("  ", OUTPUT_CLUSTER_DIAGNOSTICS_CSV, "\n", sep = "")

if (!eligible_for_downstream) {
  cat("\nDECISION: Patient 6 branch did not meet the frozen operational\n")
  cat("eligibility criteria. Preserve and report these diagnostics, but do\n")
  cat("not use this reconstructed population for downstream signatures.\n")
}

cat("============================================================\n")

# =============================================================================
# End of revision analysis 1C.
# =============================================================================
