# =============================================================================
# 02C_generate_signature_variants.R
# =============================================================================
#
# Purpose
# -------
# Generate and freeze all non-LOO differential-expression signatures required
# for the reviewer revision before any LINCS query is run.
#
# Analysis axes
# -------------
# 1. Primary C5-vs-C0 signature sizes: 50, 100, 150, and 200 per direction.
# 2. Comparator sensitivity: C5 vs all remaining canonical monocytes.
# 3. Comparator sensitivity: C5 vs marker-defined IFN-low monocytes.
# 4. Clustering-independent state: marker IFN-high vs marker IFN-low.
# 5. Resolution sensitivity: selected IFN cluster vs all remaining monocytes
#    at resolutions 0.3, 0.7, and 1.0. Resolution 0.5 is represented by the
#    C5-vs-all-rest query and is not duplicated.
# 6. Patient 6 sensitivity: marker IFN-high vs marker IFN-low in the uniform
#    ten-patient population.
#
# Primary signature-size queries use their exact fixed sizes. Other sensitivity
# contrasts use a balanced FDR-significant signature of up to 150 genes per
# direction, limited by the direction with fewer significant genes. A minimum
# of 25 genes per direction is required for downstream connectivity analysis.
#
# Every DE model uses donor-level pseudobulk counts, a paired
# ~ patient + state DESeq2 design, apeglm shrinkage, and a minimum of 30 cells
# per patient-state pseudobulk. This script performs no LINCS query.
#
# Inputs
# ------
# data/processed/seu_monocyte_subclustered.rds
# data/processed/inph_ifn_signature.rds
# data/processed/de_ifn_subtype_results.rds
# data/processed/revision_state_resolution_sensitivity.rds
# data/processed/revision_marker_ifn_states.rds
# data/processed/revision_patient6_uniform_monocytes.rds
#
# Outputs
# -------
# data/processed/revision_signature_variants.rds
# results/revision/signature_variant_summary.csv
# results/revision/signature_variant_donor_counts.csv
# results/revision/signature_variant_genes.csv
# results/revision/signature_variant_de_results.csv.gz
# =============================================================================

# -----------------------------------------------------------------------------
# Section 1 · Packages, paths, and overwrite protection
# -----------------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "Matrix",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "DESeq2",
  "apeglm",
  "dplyr",
  "tibble",
  "readr",
  "here"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

PROCESSED_DIR <- here::here("data", "processed")
REVISION_RESULTS_DIR <- here::here("results", "revision")

SEURAT_PATH <- file.path(
  PROCESSED_DIR,
  "seu_monocyte_subclustered.rds"
)
PRIMARY_SIGNATURE_PATH <- file.path(
  PROCESSED_DIR,
  "inph_ifn_signature.rds"
)
PRIMARY_DE_PATH <- file.path(
  PROCESSED_DIR,
  "de_ifn_subtype_results.rds"
)
RESOLUTION_PATH <- file.path(
  PROCESSED_DIR,
  "revision_state_resolution_sensitivity.rds"
)
MARKER_STATE_PATH <- file.path(
  PROCESSED_DIR,
  "revision_marker_ifn_states.rds"
)
UNIFORM_MONOCYTES_PATH <- file.path(
  PROCESSED_DIR,
  "revision_patient6_uniform_monocytes.rds"
)

OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_signature_variants.rds"
)
OUTPUT_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "signature_variant_summary.csv"
)
OUTPUT_DONOR_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "signature_variant_donor_counts.csv"
)
OUTPUT_GENES_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "signature_variant_genes.csv"
)
OUTPUT_DE_CSV_GZ <- file.path(
  REVISION_RESULTS_DIR,
  "signature_variant_de_results.csv.gz"
)

OUTPUT_PATHS <- c(
  OUTPUT_RDS,
  OUTPUT_SUMMARY_CSV,
  OUTPUT_DONOR_CSV,
  OUTPUT_GENES_CSV,
  OUTPUT_DE_CSV_GZ
)
required_inputs <- c(
  SEURAT_PATH,
  PRIMARY_SIGNATURE_PATH,
  PRIMARY_DE_PATH,
  RESOLUTION_PATH,
  MARKER_STATE_PATH,
  UNIFORM_MONOCYTES_PATH
)

if (!all(file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_SIGNATURE_VARIANTS"),
  "1"
)
existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "Signature-variant output(s) already exist: ",
    paste(existing_outputs, collapse = ", "),
    ". Refusing to overwrite automatically. If replacement is intentional, ",
    "set INPH_OVERWRITE_SIGNATURE_VARIANTS=1 for that run."
  )
}

# -----------------------------------------------------------------------------
# Section 2 · Frozen analysis configuration
# -----------------------------------------------------------------------------

MIN_PSEUDOBULK_CELLS <- 30L
FDR_THRESHOLD <- 0.05
PRIMARY_SIGNATURE_SIZES <- c(50L, 100L, 150L, 200L)
VARIANT_SIGNATURE_SIZE <- 150L
MIN_VARIANT_SIGNATURE_SIZE <- 25L
EXPECTED_CANONICAL_CELLS <- 22197L
EXPECTED_UNIFORM_CELLS <- 24703L
EXPECTED_CANONICAL_PATIENTS <- c(
  "1", "2", "3", "4", "5", "7", "8", "9", "10"
)
EXPECTED_UNIFORM_PATIENTS <- as.character(1:10)
EXPECTED_PRIMARY_FDR_GENES <- 1231L

numeric_order <- function(x) {
  x <- unique(as.character(x))
  x[order(as.numeric(x))]
}

jaccard_similarity <- function(a, b) {
  union_n <- length(union(a, b))
  if (union_n == 0L) {
    return(NA_real_)
  }
  length(intersect(a, b)) / union_n
}

extract_signature <- function(
    de_table,
    signature_size,
    allow_balanced_reduction = FALSE,
    minimum_size = signature_size
) {
  passing <- de_table |>
    dplyr::filter(
      !is.na(padj),
      padj < FDR_THRESHOLD,
      !is.na(log2FoldChange)
    )

  available_up <- passing |>
    dplyr::filter(log2FoldChange > 0) |>
    dplyr::arrange(dplyr::desc(log2FoldChange))

  available_down <- passing |>
    dplyr::filter(log2FoldChange < 0) |>
    dplyr::arrange(log2FoldChange)

  actual_size <- min(
    signature_size,
    nrow(available_up),
    nrow(available_down)
  )

  if (!allow_balanced_reduction && actual_size != signature_size) {
    stop(
      "Could not extract ", signature_size,
      " FDR-significant genes in each direction: observed ",
      nrow(available_up), " UP and ", nrow(available_down), " DOWN."
    )
  }

  if (actual_size < minimum_size) {
    stop(
      "Balanced signature would contain only ", actual_size,
      " genes per direction, below the frozen minimum of ",
      minimum_size, "."
    )
  }

  up <- available_up |>
    dplyr::slice_head(n = actual_size) |>
    dplyr::pull(gene)

  down <- available_down |>
    dplyr::slice_head(n = actual_size) |>
    dplyr::pull(gene)

  if (anyDuplicated(up) ||
      anyDuplicated(down) ||
      length(intersect(up, down)) > 0L) {
    stop("Extracted signature is duplicated or not directionally disjoint.")
  }

  list(
    up = up,
    down = down,
    requested_size = signature_size,
    actual_size = actual_size,
    available_fdr_up = nrow(available_up),
    available_fdr_down = nrow(available_down),
    balanced_reduction_used = actual_size < signature_size
  )
}

# -----------------------------------------------------------------------------
# Section 3 · Shared paired-pseudobulk DE function
# -----------------------------------------------------------------------------

run_paired_pseudobulk_de <- function(
    contrast_id,
    counts_mat,
    assignments,
    expected_patients
) {
  required_assignment_columns <- c(
    "cell_id", "patient", "state"
  )
  missing_columns <- setdiff(
    required_assignment_columns,
    colnames(assignments)
  )

  if (length(missing_columns) > 0L) {
    stop(
      contrast_id, " assignments are missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  assignments <- as.data.frame(
    assignments[, required_assignment_columns, drop = FALSE],
    stringsAsFactors = FALSE
  )
  assignments$cell_id <- as.character(assignments$cell_id)
  assignments$patient <- as.character(assignments$patient)
  assignments$state <- as.character(assignments$state)

  if (anyDuplicated(assignments$cell_id)) {
    stop(contrast_id, " assignments contain duplicated cell identifiers.")
  }

  unexpected_states <- setdiff(
    unique(assignments$state),
    c("low", "high")
  )
  if (length(unexpected_states) > 0L) {
    stop(
      contrast_id, " contains unexpected state label(s): ",
      paste(unexpected_states, collapse = ", ")
    )
  }

  missing_cells <- setdiff(assignments$cell_id, colnames(counts_mat))
  if (length(missing_cells) > 0L) {
    stop(
      contrast_id, " has ", length(missing_cells),
      " assigned cells absent from the count matrix."
    )
  }

  patient_state_counts <- table(
    assignments$patient,
    factor(assignments$state, levels = c("low", "high"))
  )
  patient_ids <- numeric_order(assignments$patient)

  if (!identical(patient_ids, expected_patients)) {
    stop(
      contrast_id, " patient set differs from its frozen cohort. Observed: ",
      paste(patient_ids, collapse = ", ")
    )
  }

  patient_state_counts <- patient_state_counts[
    patient_ids,
    c("low", "high"),
    drop = FALSE
  ]

  if (any(patient_state_counts < MIN_PSEUDOBULK_CELLS)) {
    failing <- which(
      patient_state_counts < MIN_PSEUDOBULK_CELLS,
      arr.ind = TRUE
    )
    stop(
      contrast_id,
      " contains patient-state pseudobulk(s) below ",
      MIN_PSEUDOBULK_CELLS,
      " cells."
    )
  }

  group_table <- expand.grid(
    state = c("low", "high"),
    patient = patient_ids,
    stringsAsFactors = FALSE
  )
  group_table <- group_table[
    order(
      match(group_table$patient, patient_ids),
      match(group_table$state, c("low", "high"))
    ),
    ,
    drop = FALSE
  ]
  group_table$group <- paste(
    group_table$patient,
    group_table$state,
    sep = "__"
  )

  assignment_group <- paste(
    assignments$patient,
    assignments$state,
    sep = "__"
  )

  pb_columns <- lapply(group_table$group, function(group_id) {
    cells <- assignments$cell_id[assignment_group == group_id]
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })
  pb_mat <- do.call(cbind, pb_columns)
  colnames(pb_mat) <- group_table$group
  rownames(pb_mat) <- rownames(counts_mat)

  if (ncol(pb_mat) != 2L * length(patient_ids)) {
    stop(contrast_id, " pseudobulk aggregation returned wrong dimensions.")
  }

  if (any(!is.finite(pb_mat)) || any(pb_mat < 0)) {
    stop(contrast_id, " pseudobulk count matrix is invalid.")
  }

  if (any(pb_mat != round(pb_mat))) {
    stop(contrast_id, " pseudobulk count matrix is not integer-valued.")
  }
  storage.mode(pb_mat) <- "integer"

  pb_coldata <- data.frame(
    patient = factor(
      group_table$patient,
      levels = patient_ids
    ),
    state = factor(
      group_table$state,
      levels = c("low", "high")
    ),
    row.names = group_table$group,
    stringsAsFactors = FALSE
  )

  if (!identical(colnames(pb_mat), rownames(pb_coldata))) {
    stop(contrast_id, " pseudobulk counts and metadata are misaligned.")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pb_mat,
    colData = pb_coldata,
    design = ~ patient + state
  )
  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE,
    parallel = FALSE
  )

  expected_coef <- "state_high_vs_low"
  if (!expected_coef %in% DESeq2::resultsNames(dds)) {
    stop(
      contrast_id, " expected coefficient not found: ",
      expected_coef
    )
  }

  shrunken <- DESeq2::lfcShrink(
    dds,
    coef = expected_coef,
    type = "apeglm",
    quiet = TRUE
  )

  de_table <- shrunken |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble() |>
    dplyr::arrange(padj)

  signature_max150 <- extract_signature(
    de_table,
    VARIANT_SIGNATURE_SIZE,
    allow_balanced_reduction = TRUE,
    minimum_size = MIN_VARIANT_SIGNATURE_SIZE
  )

  donor_counts <- data.frame(
    contrast_id = contrast_id,
    patient = patient_ids,
    low_cells = as.integer(patient_state_counts[, "low"]),
    high_cells = as.integer(patient_state_counts[, "high"]),
    stringsAsFactors = FALSE
  )

  list(
    contrast_id = contrast_id,
    de_table = de_table,
    signature_max150 = signature_max150,
    donor_counts = donor_counts,
    n_patients = length(patient_ids),
    n_low_cells = sum(assignments$state == "low"),
    n_high_cells = sum(assignments$state == "high"),
    n_fdr_significant = sum(
      !is.na(de_table$padj) &
        de_table$padj < FDR_THRESHOLD
    )
  )
}

# -----------------------------------------------------------------------------
# Section 4 · Load and align frozen inputs
# -----------------------------------------------------------------------------

seu <- readRDS(SEURAT_PATH)
if (!inherits(seu, "Seurat")) {
  stop("Canonical input must be a Seurat object.")
}
if (ncol(seu) != EXPECTED_CANONICAL_CELLS) {
  stop("Canonical Seurat object does not contain 22,197 cells.")
}

canonical_counts <- Seurat::GetAssayData(
  seu,
  layer = "counts"
)
canonical_cells <- colnames(canonical_counts)
canonical_patient <- as.character(seu$patient)
canonical_cluster <- as.character(seu$sub_res05)

if (!identical(
  numeric_order(canonical_patient),
  EXPECTED_CANONICAL_PATIENTS
)) {
  stop("Canonical Seurat patient set differs from the frozen cohort.")
}

primary_signature <- readRDS(PRIMARY_SIGNATURE_PATH)
primary_de_saved <- readRDS(PRIMARY_DE_PATH)

if (!is.list(primary_signature) ||
    !all(c("up", "down") %in% names(primary_signature)) ||
    length(primary_signature$up) != 150L ||
    length(primary_signature$down) != 150L) {
  stop("Submitted primary signature failed structure validation.")
}

if (!"gene" %in% colnames(primary_de_saved)) {
  stop("Submitted primary DE table lacks a gene column.")
}

resolution_result <- readRDS(RESOLUTION_PATH)
marker_result <- readRDS(MARKER_STATE_PATH)
uniform_monocytes <- readRDS(UNIFORM_MONOCYTES_PATH)

if (!isTRUE(
  resolution_result$input_provenance$primary_assignment_reproduced
)) {
  stop("Resolution analysis did not reproduce the primary assignments.")
}

if (!isTRUE(
  marker_result$downstream_eligibility$canonical_marker_de_eligible
)) {
  stop("Canonical marker-defined state is not eligible for paired DE.")
}
if (!isTRUE(
  marker_result$downstream_eligibility[[
    "uniform_ten_patient_marker_de_eligible"
  ]]
)) {
  stop("Ten-patient marker-defined state is not eligible for paired DE.")
}

resolution_assignments <- as.data.frame(
  resolution_result$cell_assignments,
  stringsAsFactors = FALSE
)
canonical_marker <- as.data.frame(
  marker_result$canonical_assignments,
  stringsAsFactors = FALSE
)
uniform_marker <- as.data.frame(
  marker_result$uniform_ten_patient_assignments,
  stringsAsFactors = FALSE
)

if (!identical(
  resolution_assignments$cell_id,
  canonical_cells
)) {
  stop("Resolution assignments do not match canonical cell order.")
}

canonical_marker <- canonical_marker[
  match(canonical_cells, canonical_marker$cell_id),
  ,
  drop = FALSE
]
if (!identical(canonical_marker$cell_id, canonical_cells)) {
  stop("Canonical marker assignments failed alignment.")
}

if (!inherits(uniform_monocytes, "SingleCellExperiment") ||
    ncol(uniform_monocytes) != EXPECTED_UNIFORM_CELLS) {
  stop("Uniform ten-patient monocyte input failed validation.")
}

uniform_counts <- SummarizedExperiment::assay(
  uniform_monocytes,
  "counts"
)
uniform_cells <- colnames(uniform_counts)
uniform_patient <- as.character(uniform_monocytes$patient)

uniform_marker <- uniform_marker[
  match(uniform_cells, uniform_marker$cell_id),
  ,
  drop = FALSE
]
if (!identical(uniform_marker$cell_id, uniform_cells)) {
  stop("Uniform marker assignments failed alignment.")
}

if (!identical(
  numeric_order(uniform_patient),
  EXPECTED_UNIFORM_PATIENTS
)) {
  stop("Uniform patient set differs from the frozen ten-patient cohort.")
}

# -----------------------------------------------------------------------------
# Section 5 · Construct the eight frozen DE contrasts
# -----------------------------------------------------------------------------

make_assignments <- function(cell_id, patient, state) {
  keep <- !is.na(state)
  data.frame(
    cell_id = as.character(cell_id[keep]),
    patient = as.character(patient[keep]),
    state = as.character(state[keep]),
    stringsAsFactors = FALSE
  )
}

primary_state <- ifelse(
  canonical_cluster == "5",
  "high",
  ifelse(canonical_cluster == "0", "low", NA_character_)
)

c5_allrest_state <- ifelse(
  canonical_cluster == "5",
  "high",
  "low"
)

marker_high <- canonical_marker$marker_ifn_state ==
  "marker_ifn_high"
marker_low <- canonical_marker$marker_ifn_state ==
  "marker_ifn_low"

c5_markerlow_state <- ifelse(
  canonical_cluster == "5",
  "high",
  ifelse(
    canonical_cluster != "5" & marker_low,
    "low",
    NA_character_
  )
)

marker_state <- ifelse(
  marker_high,
  "high",
  "low"
)

required_resolution_flags <- c(
  "rev_res03_ifn_high",
  "rev_res07_ifn_high",
  "rev_res10_ifn_high"
)
missing_resolution_flags <- setdiff(
  required_resolution_flags,
  colnames(resolution_assignments)
)
if (length(missing_resolution_flags) > 0L) {
  stop(
    "Resolution result is missing flag(s): ",
    paste(missing_resolution_flags, collapse = ", ")
  )
}

uniform_state <- ifelse(
  uniform_marker$marker_ifn_state == "marker_ifn_high",
  "high",
  "low"
)

contrast_specs <- list(
  PRIMARY_C5_VS_C0 = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      primary_state
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "primary_reproduction",
    description = "cluster 5 versus cluster 0",
    resolution = 0.5,
    patient6_included = FALSE
  ),
  C5_ALLREST = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      c5_allrest_state
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "comparator_and_resolution",
    description =
      "cluster 5 versus all remaining canonical monocytes",
    resolution = 0.5,
    patient6_included = FALSE
  ),
  C5_MARKERLOW = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      c5_markerlow_state
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "comparator",
    description =
      "cluster 5 versus non-cluster-5 marker IFN-low monocytes",
    resolution = 0.5,
    patient6_included = FALSE
  ),
  MARKER = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      marker_state
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "marker_defined",
    description =
      "marker IFN-high versus marker IFN-low canonical monocytes",
    resolution = NA_real_,
    patient6_included = FALSE
  ),
  RES03_ALLREST = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      ifelse(
        resolution_assignments$rev_res03_ifn_high,
        "high",
        "low"
      )
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "resolution",
    description =
      "resolution 0.3 selected IFN cluster versus all remaining",
    resolution = 0.3,
    patient6_included = FALSE
  ),
  RES07_ALLREST = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      ifelse(
        resolution_assignments$rev_res07_ifn_high,
        "high",
        "low"
      )
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "resolution",
    description =
      "resolution 0.7 selected IFN cluster versus all remaining",
    resolution = 0.7,
    patient6_included = FALSE
  ),
  RES10_ALLREST = list(
    counts = canonical_counts,
    assignments = make_assignments(
      canonical_cells,
      canonical_patient,
      ifelse(
        resolution_assignments$rev_res10_ifn_high,
        "high",
        "low"
      )
    ),
    expected_patients = EXPECTED_CANONICAL_PATIENTS,
    family = "resolution",
    description =
      "resolution 1.0 selected IFN cluster versus all remaining",
    resolution = 1.0,
    patient6_included = FALSE
  ),
  P6_INCLUDED_MARKER = list(
    counts = uniform_counts,
    assignments = make_assignments(
      uniform_cells,
      uniform_patient,
      uniform_state
    ),
    expected_patients = EXPECTED_UNIFORM_PATIENTS,
    family = "patient6",
    description = paste(
      "marker IFN-high versus marker IFN-low in uniform",
      "ten-patient population"
    ),
    resolution = NA_real_,
    patient6_included = TRUE
  )
)

# -----------------------------------------------------------------------------
# Section 6 · Run all DE models and validate the primary implementation
# -----------------------------------------------------------------------------

de_results <- setNames(
  vector("list", length(contrast_specs)),
  names(contrast_specs)
)

cat(
  "Running ", length(contrast_specs),
  " paired pseudobulk DE models ...\n",
  sep = ""
)

for (contrast_id in names(contrast_specs)) {
  spec <- contrast_specs[[contrast_id]]
  cat("  ", contrast_id, " ... ", sep = "")
  start_time <- proc.time()

  de_results[[contrast_id]] <- run_paired_pseudobulk_de(
    contrast_id = contrast_id,
    counts_mat = spec$counts,
    assignments = spec$assignments,
    expected_patients = spec$expected_patients
  )

  elapsed <- round((proc.time() - start_time)[["elapsed"]])
  cat(
    "done (", elapsed, "s; ",
    de_results[[contrast_id]]$n_fdr_significant,
    " FDR-significant genes)\n",
    sep = ""
  )
  invisible(gc(FALSE))
}

primary_recomputed <- de_results$PRIMARY_C5_VS_C0

if (primary_recomputed$n_fdr_significant !=
    EXPECTED_PRIMARY_FDR_GENES) {
  stop(
    "Primary DE reproduction returned ",
    primary_recomputed$n_fdr_significant,
    " FDR-significant genes instead of ",
    EXPECTED_PRIMARY_FDR_GENES, "."
  )
}

if (!identical(
  primary_recomputed$signature_max150$up,
  primary_signature$up
) ||
    !identical(
      primary_recomputed$signature_max150$down,
      primary_signature$down
    )) {
  stop(
    "Shared pseudobulk function did not reproduce the exact submitted ",
    "150-UP/150-DOWN signature."
  )
}

cat("Primary DE/signature implementation check: PASS\n")

# -----------------------------------------------------------------------------
# Section 7 · Build the 11 frozen query signatures
# -----------------------------------------------------------------------------

signatures <- list()
query_metadata <- list()

add_signature <- function(
    query_id,
    family,
    source_contrast,
    maximum_requested_size,
    signature,
    description,
    resolution,
    patient6_included
) {
  signatures[[query_id]] <<- list(
    up = signature$up,
    down = signature$down
  )
  query_metadata[[query_id]] <<- data.frame(
    query_id = query_id,
    family = family,
    source_contrast = source_contrast,
    maximum_requested_size = maximum_requested_size,
    signature_size_per_direction = length(signature$up),
    available_fdr_up = signature$available_fdr_up,
    available_fdr_down = signature$available_fdr_down,
    balanced_reduction_used = signature$balanced_reduction_used,
    description = description,
    resolution = resolution,
    patient6_included = patient6_included,
    stringsAsFactors = FALSE
  )
}

for (signature_size in PRIMARY_SIGNATURE_SIZES) {
  query_id <- paste0("PRIMARY_", signature_size)
  add_signature(
    query_id = query_id,
    family = "signature_size",
    source_contrast = "PRIMARY_C5_VS_C0",
    maximum_requested_size = signature_size,
    signature = extract_signature(
      primary_de_saved,
      signature_size
    ),
    description = paste0(
      "submitted cluster 5 versus cluster 0; ",
      signature_size, " genes per direction"
    ),
    resolution = 0.5,
    patient6_included = FALSE
  )
}

variant_query_map <- c(
  C5_ALLREST = "C5_ALLREST_MAX150",
  C5_MARKERLOW = "C5_MARKERLOW_MAX150",
  MARKER = "MARKER_MAX150",
  RES03_ALLREST = "RES03_ALLREST_MAX150",
  RES07_ALLREST = "RES07_ALLREST_MAX150",
  RES10_ALLREST = "RES10_ALLREST_MAX150",
  P6_INCLUDED_MARKER = "P6_INCLUDED_MARKER_MAX150"
)

for (source_contrast in names(variant_query_map)) {
  query_id <- unname(variant_query_map[[source_contrast]])
  spec <- contrast_specs[[source_contrast]]

  add_signature(
    query_id = query_id,
    family = spec$family,
    source_contrast = source_contrast,
    maximum_requested_size = VARIANT_SIGNATURE_SIZE,
    signature = de_results[[source_contrast]]$signature_max150,
    description = spec$description,
    resolution = spec$resolution,
    patient6_included = spec$patient6_included
  )
}

if (length(signatures) != 11L) {
  stop("Expected 11 non-LOO query signatures.")
}

if (!identical(
  signatures$PRIMARY_150,
  primary_signature
)) {
  stop("Frozen PRIMARY_150 signature differs from the submitted signature.")
}

# -----------------------------------------------------------------------------
# Section 8 · Assemble summaries, gene tables, and DE exports
# -----------------------------------------------------------------------------

primary_150 <- signatures$PRIMARY_150
query_metadata_df <- dplyr::bind_rows(query_metadata)

summary_rows <- lapply(names(signatures), function(query_id) {
  metadata_row <- query_metadata_df[
    query_metadata_df$query_id == query_id,
    ,
    drop = FALSE
  ]
  source_contrast <- metadata_row$source_contrast
  source_result <- de_results[[source_contrast]]
  signature <- signatures[[query_id]]

  data.frame(
    query_id = query_id,
    family = metadata_row$family,
    source_contrast = source_contrast,
    description = metadata_row$description,
    resolution = metadata_row$resolution,
    patient6_included = metadata_row$patient6_included,
    patients_n = source_result$n_patients,
    high_cells_n = source_result$n_high_cells,
    low_cells_n = source_result$n_low_cells,
    fdr_significant_genes_n =
      source_result$n_fdr_significant,
    available_fdr_up = metadata_row$available_fdr_up,
    available_fdr_down = metadata_row$available_fdr_down,
    maximum_requested_size = metadata_row$maximum_requested_size,
    balanced_reduction_used = metadata_row$balanced_reduction_used,
    signature_up_n = length(signature$up),
    signature_down_n = length(signature$down),
    up_overlap_primary_n =
      length(intersect(signature$up, primary_150$up)),
    down_overlap_primary_n =
      length(intersect(signature$down, primary_150$down)),
    up_jaccard_primary =
      jaccard_similarity(signature$up, primary_150$up),
    down_jaccard_primary =
      jaccard_similarity(signature$down, primary_150$down),
    stringsAsFactors = FALSE
  )
})

signature_summary <- dplyr::bind_rows(summary_rows)

signature_genes <- dplyr::bind_rows(lapply(
  names(signatures),
  function(query_id) {
    signature <- signatures[[query_id]]
    dplyr::bind_rows(
      data.frame(
        query_id = query_id,
        direction = "UP",
        rank = seq_along(signature$up),
        gene = signature$up,
        stringsAsFactors = FALSE
      ),
      data.frame(
        query_id = query_id,
        direction = "DOWN",
        rank = seq_along(signature$down),
        gene = signature$down,
        stringsAsFactors = FALSE
      )
    )
  }
))

donor_counts <- dplyr::bind_rows(lapply(
  de_results,
  function(x) x$donor_counts
))

combined_de <- dplyr::bind_rows(lapply(
  names(de_results),
  function(contrast_id) {
    dplyr::mutate(
      de_results[[contrast_id]]$de_table,
      contrast_id = contrast_id,
      .before = 1
    )
  }
))

primary_de_common <- dplyr::inner_join(
  dplyr::select(
    as.data.frame(primary_de_saved),
    gene,
    submitted_lfc = log2FoldChange,
    submitted_padj = padj
  ),
  dplyr::select(
    as.data.frame(primary_recomputed$de_table),
    gene,
    recomputed_lfc = log2FoldChange,
    recomputed_padj = padj
  ),
  by = "gene"
)

primary_reproduction <- list(
  fdr_significant_genes =
    primary_recomputed$n_fdr_significant,
  submitted_150_up_identical = identical(
    primary_recomputed$signature_max150$up,
    primary_signature$up
  ),
  submitted_150_down_identical = identical(
    primary_recomputed$signature_max150$down,
    primary_signature$down
  ),
  lfc_pearson = stats::cor(
    primary_de_common$submitted_lfc,
    primary_de_common$recomputed_lfc,
    use = "complete.obs",
    method = "pearson"
  ),
  lfc_max_absolute_difference = max(
    abs(
      primary_de_common$submitted_lfc -
        primary_de_common$recomputed_lfc
    ),
    na.rm = TRUE
  ),
  padj_max_absolute_difference = max(
    abs(
      primary_de_common$submitted_padj -
        primary_de_common$recomputed_padj
    ),
    na.rm = TRUE
  )
)

# -----------------------------------------------------------------------------
# Section 9 · Save and validate frozen outputs
# -----------------------------------------------------------------------------

result <- list(
  analysis_id = "revision_02C_signature_variants",
  configuration = list(
    minimum_pseudobulk_cells =
      MIN_PSEUDOBULK_CELLS,
    fdr_threshold = FDR_THRESHOLD,
    primary_signature_sizes =
      PRIMARY_SIGNATURE_SIZES,
    variant_maximum_signature_size =
      VARIANT_SIGNATURE_SIZE,
    minimum_variant_signature_size =
      MIN_VARIANT_SIGNATURE_SIZE,
    design = "~ patient + state",
    lfc_shrinkage = "apeglm",
    extraction_rule = paste(
      "FDR < 0.05; ranked separately by shrunken log2FC;",
      "primary sizes fixed; other contrasts balanced to the smaller",
      "available direction, capped at 150 genes per direction"
    ),
    resolution_05_query_note = paste(
      "C5_ALLREST_MAX150 is also the resolution-0.5 all-rest query;",
      "it is not duplicated under a second query ID"
    )
  ),
  signatures = signatures,
  query_metadata = query_metadata_df,
  signature_summary = signature_summary,
  donor_counts = donor_counts,
  de_results = lapply(
    de_results,
    function(x) x$de_table
  ),
  primary_reproduction = primary_reproduction,
  input_provenance = list(
    seurat_md5 = unname(tools::md5sum(SEURAT_PATH)),
    primary_signature_md5 =
      unname(tools::md5sum(PRIMARY_SIGNATURE_PATH)),
    primary_de_md5 =
      unname(tools::md5sum(PRIMARY_DE_PATH)),
    resolution_result_md5 =
      unname(tools::md5sum(RESOLUTION_PATH)),
    marker_state_result_md5 =
      unname(tools::md5sum(MARKER_STATE_PATH)),
    uniform_monocytes_md5 =
      unname(tools::md5sum(UNIFORM_MONOCYTES_PATH))
  )
)

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(
  REVISION_RESULTS_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(result, OUTPUT_RDS)
readr::write_csv(signature_summary, OUTPUT_SUMMARY_CSV)
readr::write_csv(donor_counts, OUTPUT_DONOR_CSV)
readr::write_csv(signature_genes, OUTPUT_GENES_CSV)
readr::write_csv(combined_de, OUTPUT_DE_CSV_GZ)

if (!all(file.exists(OUTPUT_PATHS))) {
  stop("One or more signature-variant outputs were not created.")
}

if (any(is.na(file.info(OUTPUT_PATHS)$size)) ||
    any(file.info(OUTPUT_PATHS)$size <= 0)) {
  stop("One or more signature-variant outputs are empty.")
}

result_check <- readRDS(OUTPUT_RDS)

if (!identical(
  names(result_check$signatures),
  names(signatures)
)) {
  stop("Saved signature query order failed validation.")
}

if (!identical(
  result_check$signatures$PRIMARY_150,
  primary_signature
)) {
  stop("Saved PRIMARY_150 signature failed validation.")
}

if (nrow(readr::read_csv(
  OUTPUT_SUMMARY_CSV,
  show_col_types = FALSE
)) != 11L) {
  stop("Saved signature summary should contain 11 queries.")
}

# -----------------------------------------------------------------------------
# Section 10 · Closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 2C COMPLETE\n")
cat("============================================================\n")
cat("  Paired pseudobulk DE models : ", length(de_results), "\n", sep = "")
cat("  Frozen non-LOO queries      : ", length(signatures), "\n", sep = "")
cat("  Primary reproduction        : PASS\n")
cat("  Primary FDR-significant     : ",
    primary_recomputed$n_fdr_significant, "\n", sep = "")
cat("\nSignature summary:\n")
print(signature_summary, row.names = FALSE, digits = 4)
cat("\nPrimary reproduction diagnostics:\n")
print(primary_reproduction)
cat("\nOutputs saved:\n")
cat("  ", OUTPUT_RDS, "\n", sep = "")
cat("  ", OUTPUT_SUMMARY_CSV, "\n", sep = "")
cat("  ", OUTPUT_DONOR_CSV, "\n", sep = "")
cat("  ", OUTPUT_GENES_CSV, "\n", sep = "")
cat("  ", OUTPUT_DE_CSV_GZ, "\n", sep = "")
cat("============================================================\n")
