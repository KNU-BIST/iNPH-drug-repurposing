# =============================================================================
# Revision analysis 1B — Clustering-resolution sensitivity
# Script: 01B_clustering_resolution_sensitivity.R
# =============================================================================
#
# Project
# -------
# Systems biology prioritizes drug candidates for interferon-activated
# monocytes in idiopathic normal pressure hydrocephalus
#
# Purpose
# -------
# Test whether an interferon-response-high classical-monocyte population
# persists across a revision-stage resolution grid without changing the
# submitted primary analysis. The IFN-high cluster at each resolution is
# selected automatically as the cluster with the highest mean score for the
# same predefined seven-gene interferon panel used in Script 01.
#
# Revision-stage analysis choices (frozen before execution)
# ---------------------------------------------------------
# - Input population: the 22,197-cell, nine-patient canonical monocyte object
#   subclustered by Script 01.
# - Existing neighbor graph: RNA_snn, generated from 30 Harmony dimensions.
# - Resolution grid: 0.3, 0.5, 0.7, 1.0.
# - Clustering: Louvain algorithm 1, random.seed = 0.
# - IFN selector: IFITM3, IFITM2, ISG15, MX1, OAS1, IRF7, STAT1.
# - Automatic selection: highest unrounded mean seven-gene score; a tie is
#   resolved by the smallest numeric cluster label, then lexical label.
# - Submitted reference state: resolution 0.5, cluster 5 (940 cells).
#
# Dependency boundary
# -------------------
# This script evaluates clustering-resolution robustness only. It does not
# derive differential-expression signatures and does not run LINCS. The
# marker-defined IFN-high analysis is intentionally separate because its
# donor-balanced healthy-control-anchored threshold is produced by Script 02.
#
# Input
# -----
# data/processed/seu_monocyte_subclustered.rds
#
# Outputs
# -------
# data/processed/revision_state_resolution_sensitivity.rds
# results/revision/state_resolution_summary.csv
# results/revision/state_resolution_cluster_diagnostics.csv
# results/revision/state_resolution_donor_representation.csv
# results/revision/state_resolution_pairwise_overlap.csv
#
# The RDS is the authoritative machine-readable output. It contains the frozen
# configuration, selected IFN-high cluster at each resolution, all cluster
# diagnostics, donor representation, pairwise overlaps, and per-cell cluster
# assignments needed by later signature-generation scripts.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment and frozen configuration
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "Seurat",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "Matrix",
  "dplyr",
  "tibble"
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
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tibble)
})

PROCESSED_DIR <- here::here("data", "processed")
REVISION_RESULTS_DIR <- here::here("results", "revision")

INPUT_SEURAT <- file.path(
  PROCESSED_DIR,
  "seu_monocyte_subclustered.rds"
)

OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_state_resolution_sensitivity.rds"
)

OUTPUT_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "state_resolution_summary.csv"
)

OUTPUT_DIAGNOSTICS_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "state_resolution_cluster_diagnostics.csv"
)

OUTPUT_DONOR_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "state_resolution_donor_representation.csv"
)

OUTPUT_PAIRWISE_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "state_resolution_pairwise_overlap.csv"
)

RESOLUTION_GRID <- c(0.3, 0.5, 0.7, 1.0)
GRAPH_NAME <- "RNA_snn"
CLUSTER_ALGORITHM <- 1L
CLUSTER_SEED <- 0L
PRIMARY_RESOLUTION <- 0.5
PRIMARY_CLUSTER <- "5"
PRIMARY_CLUSTER_N <- 940L
PRIMARY_REFERENCE_CLUSTER <- "0"
PRIMARY_REFERENCE_N <- 4904L
EXPECTED_CELLS <- 22197L
EXPECTED_PATIENTS <- 9L

PANELS <- list(
  classical         = c("CD14", "S100A8", "S100A9", "S100A12", "VCAN", "FCN1"),
  non_classical     = c("FCGR3A", "MS4A7", "CDKN1C", "LST1", "CSF1R", "RHOC"),
  intermediate      = c("HLA-DRA", "HLA-DRB1", "HLA-DPB1", "CD74", "CIITA"),
  proinflammatory   = c("TNF", "IL1B", "IL6", "NLRP3", "CCL3", "CCL4",
                        "CXCL8", "PTGS2", "IL1A", "CCL20"),
  interferon        = c("IFITM3", "IFITM2", "ISG15", "MX1", "OAS1",
                        "IRF7", "STAT1"),
  anti_inflammatory = c("CD163", "MRC1", "MAF", "MAFB", "SELENOP", "STAB1")
)

IFN_PANEL <- PANELS$interferon

resolution_tag <- function(resolution) {
  sprintf("%02d", as.integer(round(10 * resolution)))
}

resolution_column <- function(resolution) {
  paste0("rev_res", resolution_tag(resolution))
}

natural_label_order <- function(labels) {
  numeric_labels <- suppressWarnings(as.numeric(labels))
  order(
    is.na(numeric_labels),
    numeric_labels,
    labels,
    method = "radix"
  )
}

OUTPUT_PATHS <- c(
  OUTPUT_RDS,
  OUTPUT_SUMMARY_CSV,
  OUTPUT_DIAGNOSTICS_CSV,
  OUTPUT_DONOR_CSV,
  OUTPUT_PAIRWISE_CSV
)

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_REVISION_STATE", unset = "0"),
  "1"
)

existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "Revision state-resolution output(s) already exist:\n  ",
    paste(existing_outputs, collapse = "\n  "),
    "\nRefusing to overwrite the frozen outputs automatically. ",
    "If replacement is intentional, set ",
    "INPH_OVERWRITE_REVISION_STATE=1 for that run."
  )
}


# -----------------------------------------------------------------------------
# Section 1 · Load and validate the submitted primary object
# -----------------------------------------------------------------------------

if (
  !file.exists(INPUT_SEURAT) ||
  is.na(file.info(INPUT_SEURAT)$size) ||
  file.info(INPUT_SEURAT)$size <= 0
) {
  stop(
    "Required Script 01 output is missing or empty: ",
    INPUT_SEURAT,
    "\nRun analysis/01_derive_inph_ifn_signature.R first."
  )
}

seu <- readRDS(INPUT_SEURAT)

if (!inherits(seu, "Seurat")) {
  stop(
    "Script 01 output must be a Seurat object; observed class: ",
    paste(class(seu), collapse = ", ")
  )
}

if (!identical(Seurat::DefaultAssay(seu), "RNA")) {
  stop(
    "Expected RNA to be the default assay; observed ",
    Seurat::DefaultAssay(seu), "."
  )
}

required_metadata <- c("patient", "sub_res05")
missing_metadata <- setdiff(required_metadata, colnames(seu[[]]))

if (length(missing_metadata) > 0L) {
  stop(
    "Script 01 output is missing required metadata field(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

if (ncol(seu) != EXPECTED_CELLS) {
  stop(
    "Expected ", EXPECTED_CELLS,
    " canonical monocytes; observed ", ncol(seu), "."
  )
}

if (length(unique(as.character(seu$patient))) != EXPECTED_PATIENTS) {
  stop(
    "Expected ", EXPECTED_PATIENTS,
    " patients in the canonical object."
  )
}

if (!GRAPH_NAME %in% names(seu@graphs)) {
  stop(
    "Required shared-nearest-neighbor graph is missing: ",
    GRAPH_NAME
  )
}

if (nrow(seu@graphs[[GRAPH_NAME]]) != EXPECTED_CELLS) {
  stop("RNA_snn graph dimensions do not match the canonical cell population.")
}

primary_ids <- as.character(seu$sub_res05)
primary_cluster_cells <- colnames(seu)[primary_ids == PRIMARY_CLUSTER]
primary_reference_cells <- colnames(seu)[primary_ids == PRIMARY_REFERENCE_CLUSTER]

if (length(primary_cluster_cells) != PRIMARY_CLUSTER_N) {
  stop(
    "Submitted cluster 5 should contain ", PRIMARY_CLUSTER_N,
    " cells; observed ", length(primary_cluster_cells), "."
  )
}

if (length(primary_reference_cells) != PRIMARY_REFERENCE_N) {
  stop(
    "Submitted cluster 0 should contain ", PRIMARY_REFERENCE_N,
    " cells; observed ", length(primary_reference_cells), "."
  )
}

cat(
  "Validated Script 01 input: ",
  nrow(seu), " genes x ", ncol(seu), " cells | ",
  length(unique(as.character(seu$patient))), " patients\n",
  sep = ""
)
cat(
  "Submitted reference: resolution 0.5 cluster 5 = ",
  length(primary_cluster_cells), " cells\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 2 · Recluster the existing RNA_snn graph at the frozen grid
# -----------------------------------------------------------------------------

resolution_columns <- setNames(
  vapply(RESOLUTION_GRID, resolution_column, character(1)),
  as.character(RESOLUTION_GRID)
)

for (resolution in RESOLUTION_GRID) {
  cluster_col <- unname(resolution_columns[as.character(resolution)])

  seu <- Seurat::FindClusters(
    object = seu,
    graph.name = GRAPH_NAME,
    resolution = resolution,
    algorithm = CLUSTER_ALGORITHM,
    random.seed = CLUSTER_SEED,
    cluster.name = cluster_col,
    verbose = FALSE
  )

  if (!cluster_col %in% colnames(seu[[]])) {
    stop(
      "FindClusters did not create the expected metadata column: ",
      cluster_col
    )
  }
}

# Re-running the submitted resolution with the submitted graph and seed must
# recover the exact submitted assignments. This is a hard reproducibility check.
rerun_primary_col <- unname(
  resolution_columns[as.character(PRIMARY_RESOLUTION)]
)

if (!identical(
  as.character(seu[[rerun_primary_col, drop = TRUE]]),
  primary_ids
)) {
  mismatch_n <- sum(
    as.character(seu[[rerun_primary_col, drop = TRUE]]) != primary_ids
  )
  stop(
    "Resolution 0.5 assignments did not reproduce sub_res05; mismatched cells: ",
    mismatch_n,
    ". Do not continue until software and graph provenance are resolved."
  )
}

cat("Resolution 0.5 assignment reproduction: PASS\n")


# -----------------------------------------------------------------------------
# Section 3 · Calculate the predefined panel scores
# -----------------------------------------------------------------------------

sce <- Seurat::as.SingleCellExperiment(seu)
logc <- SingleCellExperiment::logcounts(sce)

missing_ifn_genes <- setdiff(IFN_PANEL, rownames(logc))

if (length(missing_ifn_genes) > 0L) {
  stop(
    "The predefined seven-gene IFN panel is incomplete in the input object: ",
    paste(missing_ifn_genes, collapse = ", ")
  )
}

score_panel_per_cell <- function(genes) {
  present_genes <- intersect(genes, rownames(logc))

  if (length(present_genes) == 0L) {
    return(rep(NA_real_, ncol(logc)))
  }

  Matrix::colMeans(logc[present_genes, , drop = FALSE])
}

panel_cell_scores <- lapply(PANELS, score_panel_per_cell)
panel_cell_scores <- as.data.frame(
  panel_cell_scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(panel_cell_scores) <- colnames(seu)

if (!identical(rownames(panel_cell_scores), colnames(seu))) {
  stop("Cell-score order does not match the Seurat object.")
}


# -----------------------------------------------------------------------------
# Section 4 · Select the IFN-high cluster and quantify stability
# -----------------------------------------------------------------------------

cluster_diagnostics_list <- vector("list", length(RESOLUTION_GRID))
selected_summary_list <- vector("list", length(RESOLUTION_GRID))
donor_representation_list <- vector("list", length(RESOLUTION_GRID))
selected_cell_sets <- vector("list", length(RESOLUTION_GRID))

names(cluster_diagnostics_list) <- as.character(RESOLUTION_GRID)
names(selected_summary_list) <- as.character(RESOLUTION_GRID)
names(donor_representation_list) <- as.character(RESOLUTION_GRID)
names(selected_cell_sets) <- as.character(RESOLUTION_GRID)

all_patient_ids <- unique(as.character(seu$patient))
all_patient_ids <- all_patient_ids[natural_label_order(all_patient_ids)]
patient_total_counts <- table(as.character(seu$patient))

for (i in seq_along(RESOLUTION_GRID)) {
  resolution <- RESOLUTION_GRID[[i]]
  resolution_key <- as.character(resolution)
  cluster_col <- unname(resolution_columns[resolution_key])
  cluster_ids <- as.character(seu[[cluster_col, drop = TRUE]])
  cluster_labels <- unique(cluster_ids)
  cluster_labels <- cluster_labels[natural_label_order(cluster_labels)]

  diagnostic <- data.frame(
    resolution = resolution,
    cluster_column = cluster_col,
    cluster = cluster_labels,
    n_cells = as.integer(table(cluster_ids)[cluster_labels]),
    pct_all_monocytes = 100 * as.integer(table(cluster_ids)[cluster_labels]) /
      EXPECTED_CELLS,
    stringsAsFactors = FALSE
  )

  diagnostic$n_donors <- vapply(
    cluster_labels,
    function(cluster_label) {
      length(unique(as.character(seu$patient)[cluster_ids == cluster_label]))
    },
    integer(1)
  )

  diagnostic$n_donors_ge_30_cells <- vapply(
    cluster_labels,
    function(cluster_label) {
      donor_counts <- table(
        as.character(seu$patient)[cluster_ids == cluster_label]
      )
      sum(donor_counts >= 30L)
    },
    integer(1)
  )

  for (panel_name in names(PANELS)) {
    cell_scores <- panel_cell_scores[[panel_name]]
    cluster_scores <- tapply(cell_scores, cluster_ids, mean)
    diagnostic[[paste0(panel_name, "_score")]] <- as.numeric(
      cluster_scores[cluster_labels]
    )
  }

  if (anyNA(diagnostic$interferon_score)) {
    stop(
      "Missing IFN-panel score at resolution ",
      resolution, "."
    )
  }

  cluster_number <- suppressWarnings(as.numeric(diagnostic$cluster))
  selection_order <- order(
    -diagnostic$interferon_score,
    is.na(cluster_number),
    cluster_number,
    diagnostic$cluster,
    method = "radix"
  )

  selected_row <- selection_order[[1L]]
  selected_cluster <- diagnostic$cluster[[selected_row]]
  selected_cells <- colnames(seu)[cluster_ids == selected_cluster]
  selected_cell_sets[[resolution_key]] <- selected_cells

  next_best_score <- if (length(selection_order) >= 2L) {
    diagnostic$interferon_score[selection_order[[2L]]]
  } else {
    NA_real_
  }

  overlap_n <- length(intersect(selected_cells, primary_cluster_cells))
  union_n <- length(union(selected_cells, primary_cluster_cells))

  selected_summary_list[[resolution_key]] <- data.frame(
    resolution = resolution,
    cluster_column = cluster_col,
    n_clusters = length(cluster_labels),
    selected_ifn_cluster = selected_cluster,
    selected_n_cells = length(selected_cells),
    selected_pct_all_monocytes = 100 * length(selected_cells) /
      EXPECTED_CELLS,
    selected_n_donors = diagnostic$n_donors[[selected_row]],
    selected_n_donors_ge_30_cells =
      diagnostic$n_donors_ge_30_cells[[selected_row]],
    selected_ifn_score = diagnostic$interferon_score[[selected_row]],
    next_best_ifn_score = next_best_score,
    ifn_score_margin = diagnostic$interferon_score[[selected_row]] -
      next_best_score,
    overlap_with_primary_c5_n = overlap_n,
    union_with_primary_c5_n = union_n,
    jaccard_with_primary_c5 = overlap_n / union_n,
    selected_cells_in_primary_c5_pct = 100 * overlap_n /
      length(selected_cells),
    primary_c5_recovered_pct = 100 * overlap_n /
      length(primary_cluster_cells),
    stringsAsFactors = FALSE
  )

  selected_patient_counts <- table(
    factor(
      as.character(seu$patient)[cluster_ids == selected_cluster],
      levels = all_patient_ids
    )
  )

  donor_representation_list[[resolution_key]] <- data.frame(
    resolution = resolution,
    selected_ifn_cluster = selected_cluster,
    patient = all_patient_ids,
    selected_n_cells = as.integer(selected_patient_counts[all_patient_ids]),
    patient_total_monocytes = as.integer(
      patient_total_counts[all_patient_ids]
    ),
    selected_pct_of_patient_monocytes = 100 *
      as.integer(selected_patient_counts[all_patient_ids]) /
      as.integer(patient_total_counts[all_patient_ids]),
    meets_30_cell_pseudobulk_minimum =
      as.integer(selected_patient_counts[all_patient_ids]) >= 30L,
    stringsAsFactors = FALSE
  )

  cluster_diagnostics_list[[resolution_key]] <- diagnostic
}

cluster_diagnostics <- dplyr::bind_rows(cluster_diagnostics_list)
resolution_summary <- dplyr::bind_rows(selected_summary_list)
donor_representation <- dplyr::bind_rows(donor_representation_list)

primary_summary <- resolution_summary[
  resolution_summary$resolution == PRIMARY_RESOLUTION,
  ,
  drop = FALSE
]

if (nrow(primary_summary) != 1L) {
  stop("Could not identify the single resolution 0.5 summary row.")
}

if (!identical(primary_summary$selected_ifn_cluster, PRIMARY_CLUSTER)) {
  stop(
    "Automatic selection at resolution 0.5 did not recover submitted cluster 5."
  )
}

if (primary_summary$selected_n_cells != PRIMARY_CLUSTER_N) {
  stop("Automatic resolution 0.5 selection did not recover all 940 cells.")
}

if (abs(primary_summary$jaccard_with_primary_c5 - 1) > .Machine$double.eps^0.5) {
  stop("Resolution 0.5 selected-cell Jaccard with submitted cluster 5 is not 1.")
}


# -----------------------------------------------------------------------------
# Section 5 · Pairwise overlap among automatically selected IFN-high clusters
# -----------------------------------------------------------------------------

pairwise_grid <- expand.grid(
  resolution_1 = RESOLUTION_GRID,
  resolution_2 = RESOLUTION_GRID,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pairwise_grid <- pairwise_grid[
  pairwise_grid$resolution_1 <= pairwise_grid$resolution_2,
  ,
  drop = FALSE
]

pairwise_overlap <- dplyr::bind_rows(lapply(
  seq_len(nrow(pairwise_grid)),
  function(i) {
    resolution_1 <- pairwise_grid$resolution_1[[i]]
    resolution_2 <- pairwise_grid$resolution_2[[i]]
    cells_1 <- selected_cell_sets[[as.character(resolution_1)]]
    cells_2 <- selected_cell_sets[[as.character(resolution_2)]]
    overlap_n <- length(intersect(cells_1, cells_2))
    union_n <- length(union(cells_1, cells_2))

    data.frame(
      resolution_1 = resolution_1,
      resolution_2 = resolution_2,
      selected_n_1 = length(cells_1),
      selected_n_2 = length(cells_2),
      overlap_n = overlap_n,
      union_n = union_n,
      jaccard = overlap_n / union_n,
      pct_resolution_1_recovered = 100 * overlap_n / length(cells_1),
      pct_resolution_2_recovered = 100 * overlap_n / length(cells_2),
      stringsAsFactors = FALSE
    )
  }
))


# -----------------------------------------------------------------------------
# Section 6 · Save the locked resolution-sensitivity outputs
# -----------------------------------------------------------------------------

cell_assignments <- data.frame(
  cell_id = colnames(seu),
  patient = as.character(seu$patient),
  submitted_sub_res05 = primary_ids,
  interferon_panel_score = panel_cell_scores$interferon,
  stringsAsFactors = FALSE
)

for (resolution in RESOLUTION_GRID) {
  cluster_col <- unname(resolution_columns[as.character(resolution)])
  cell_assignments[[cluster_col]] <- as.character(
    seu[[cluster_col, drop = TRUE]]
  )
}

selected_cluster_lookup <- setNames(
  resolution_summary$selected_ifn_cluster,
  as.character(resolution_summary$resolution)
)

for (resolution in RESOLUTION_GRID) {
  cluster_col <- unname(resolution_columns[as.character(resolution)])
  flag_col <- paste0(cluster_col, "_ifn_high")
  cell_assignments[[flag_col]] <-
    cell_assignments[[cluster_col]] ==
    selected_cluster_lookup[[as.character(resolution)]]
}

result <- list(
  analysis_id = "revision_01B_clustering_resolution_sensitivity",
  configuration = list(
    resolution_grid = RESOLUTION_GRID,
    graph_name = GRAPH_NAME,
    clustering_algorithm = CLUSTER_ALGORITHM,
    clustering_seed = CLUSTER_SEED,
    primary_resolution = PRIMARY_RESOLUTION,
    primary_cluster = PRIMARY_CLUSTER,
    primary_cluster_n = PRIMARY_CLUSTER_N,
    primary_reference_cluster = PRIMARY_REFERENCE_CLUSTER,
    primary_reference_n = PRIMARY_REFERENCE_N,
    expected_cells = EXPECTED_CELLS,
    expected_patients = EXPECTED_PATIENTS,
    panels = PANELS,
    selection_rule = paste(
      "highest unrounded mean seven-gene IFN-panel score;",
      "ties resolved by smallest numeric then lexical cluster label"
    )
  ),
  resolution_summary = resolution_summary,
  cluster_diagnostics = cluster_diagnostics,
  donor_representation = donor_representation,
  pairwise_overlap = pairwise_overlap,
  cell_assignments = cell_assignments,
  selected_cell_ids = selected_cell_sets,
  input_provenance = list(
    seurat_input_path = "data/processed/seu_monocyte_subclustered.rds",
    seurat_input_md5 = unname(tools::md5sum(INPUT_SEURAT)),
    primary_assignment_reproduced = TRUE
  )
)

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REVISION_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(result, OUTPUT_RDS)

utils::write.csv(
  resolution_summary,
  OUTPUT_SUMMARY_CSV,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  cluster_diagnostics,
  OUTPUT_DIAGNOSTICS_CSV,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  donor_representation,
  OUTPUT_DONOR_CSV,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  pairwise_overlap,
  OUTPUT_PAIRWISE_CSV,
  row.names = FALSE,
  na = ""
)

if (!all(file.exists(OUTPUT_PATHS))) {
  stop("One or more resolution-sensitivity outputs were not created.")
}

if (any(is.na(file.info(OUTPUT_PATHS)$size)) ||
    any(file.info(OUTPUT_PATHS)$size <= 0)) {
  stop("One or more resolution-sensitivity outputs are empty.")
}

result_check <- readRDS(OUTPUT_RDS)

if (!identical(
  result_check$configuration$resolution_grid,
  RESOLUTION_GRID
)) {
  stop("Saved resolution grid failed read-back validation.")
}

if (!identical(
  result_check$cell_assignments$cell_id,
  colnames(seu)
)) {
  stop("Saved cell assignments failed read-back order validation.")
}


# -----------------------------------------------------------------------------
# Section 7 · Closing summary
# -----------------------------------------------------------------------------

print_summary <- resolution_summary
print_summary$selected_ifn_score <- round(
  print_summary$selected_ifn_score,
  3
)
print_summary$ifn_score_margin <- round(
  print_summary$ifn_score_margin,
  3
)
print_summary$jaccard_with_primary_c5 <- round(
  print_summary$jaccard_with_primary_c5,
  3
)

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 1B COMPLETE\n")
cat("============================================================\n")
print(
  print_summary[, c(
    "resolution",
    "n_clusters",
    "selected_ifn_cluster",
    "selected_n_cells",
    "selected_n_donors",
    "selected_n_donors_ge_30_cells",
    "selected_ifn_score",
    "ifn_score_margin",
    "jaccard_with_primary_c5"
  )],
  row.names = FALSE
)
cat("\nPrimary resolution 0.5 assignment reproduction: PASS\n")
cat("Automatic resolution 0.5 IFN-cluster selection: cluster 5\n")
cat("Outputs saved:\n")
cat("  ", OUTPUT_RDS, "\n", sep = "")
cat("  ", OUTPUT_SUMMARY_CSV, "\n", sep = "")
cat("  ", OUTPUT_DIAGNOSTICS_CSV, "\n", sep = "")
cat("  ", OUTPUT_DONOR_CSV, "\n", sep = "")
cat("  ", OUTPUT_PAIRWISE_CSV, "\n", sep = "")
cat("============================================================\n")

# =============================================================================
# End of revision analysis 1B.
# =============================================================================
