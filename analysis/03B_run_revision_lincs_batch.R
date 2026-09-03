# =============================================================================
# Revision Analysis 03B — Manifest-driven LINCS2 batch analysis
# Script: 03B_run_revision_lincs_batch.R
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
# Execute the locked reviewer-revision LINCS2 queries without changing any
# upstream cell-state or signature definition. The runner is checkpointed and
# resumable: each completed query is stored in its own validated RDS file and
# is skipped on later runs only when its manifest and execution hashes match.
#
# Default safety behavior
# -----------------------
# The default mode is "validate". It runs PRIMARY_150 once against the complete
# THP1 trt_cp subset and verifies that WTCS and NCS reproduce the corresponding
# THP1 rows from the submitted full-reference result. It does not launch the
# 19 long full-reference sensitivity searches or 100 null searches.
#
# Run modes
# ---------
# INPH_REVISION_LINCS_MODE=validate  subset equivalence test only (default)
# INPH_REVISION_LINCS_MODE=full      validation + 19 full-reference queries
# INPH_REVISION_LINCS_MODE=null      validation + 100 THP1-subset null queries
# INPH_REVISION_LINCS_MODE=all       validation + both query families
#
# Optional execution batching
# ---------------------------
# INPH_REVISION_LINCS_MAX_NEW_QUERIES can limit the number of new query
# checkpoints created in one invocation. Completed matching checkpoints are
# skipped and do not count toward the limit. Leave unset to run all pending
# queries selected by the mode.
#
# Overwrite behavior
# ------------------
# Matching completed checkpoints are always reused. A checkpoint whose stored
# hashes disagree with the locked manifest causes a hard stop. Intentional
# replacement of mismatched checkpoints requires:
# INPH_OVERWRITE_REVISION_LINCS=1
#
# Inputs
# ------
# data/processed/revision_signature_manifest.rds
# results/revision/revision_signature_lock.txt
# data/processed/lincs2_gess_result.rds
# data/processed/hallmark_ifn_sensitivity_result.rds
# ExperimentHub resource EH7297
#
# Outputs
# -------
# data/processed/revision_lincs_query_results/*.rds
# data/processed/revision_lincs_primary_thp1_subset_validation.rds
# data/processed/revision_lincs_batch_index.rds
# results/revision/revision_lincs_batch_status.csv
#
# Statistical interpretation
# --------------------------
# Full-reference biological queries retain WTCS, NCS, and WTCS FDR. Tau is not
# computed or interpreted for new queries. THP1-subset null searches are used
# only to build compound-specific empirical NCS distributions; their
# subset-adjusted WTCS FDR values must not be interpreted.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment, configuration, and paths
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "readr",
  "digest",
  "signatureSearch",
  "ExperimentHub",
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

PROCESSED_DIR <- here::here("data", "processed")
RAW_DIR <- here::here("data", "raw")
REVISION_RESULTS_DIR <- here::here("results", "revision")
QUERY_RESULT_DIR <- file.path(
  PROCESSED_DIR,
  "revision_lincs_query_results"
)

MANIFEST_PATH <- file.path(
  PROCESSED_DIR,
  "revision_signature_manifest.rds"
)
LOCK_PATH <- file.path(
  REVISION_RESULTS_DIR,
  "revision_signature_lock.txt"
)
PRIMARY_RESULT_PATH <- file.path(
  PROCESSED_DIR,
  "lincs2_gess_result.rds"
)
HALLMARK_RESULT_PATH <- file.path(
  PROCESSED_DIR,
  "hallmark_ifn_sensitivity_result.rds"
)
SUBSET_VALIDATION_PATH <- file.path(
  PROCESSED_DIR,
  "revision_lincs_primary_thp1_subset_validation.rds"
)
BATCH_INDEX_PATH <- file.path(
  PROCESSED_DIR,
  "revision_lincs_batch_index.rds"
)
BATCH_STATUS_PATH <- file.path(
  REVISION_RESULTS_DIR,
  "revision_lincs_batch_status.csv"
)

required_inputs <- c(
  MANIFEST_PATH,
  LOCK_PATH,
  PRIMARY_RESULT_PATH,
  HALLMARK_RESULT_PATH
)

input_info <- file.info(required_inputs)
bad_inputs <- required_inputs[
  !file.exists(required_inputs) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(bad_inputs) > 0L) {
  stop(
    "Required LINCS batch input(s) missing or empty: ",
    paste(bad_inputs, collapse = ", ")
  )
}

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REVISION_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QUERY_RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

run_mode <- tolower(trimws(Sys.getenv(
  "INPH_REVISION_LINCS_MODE",
  unset = "validate"
)))

allowed_modes <- c("validate", "full", "null", "all")

if (!run_mode %in% allowed_modes) {
  stop(
    "INPH_REVISION_LINCS_MODE must be one of: ",
    paste(allowed_modes, collapse = ", "),
    "."
  )
}

max_new_queries_text <- trimws(Sys.getenv(
  "INPH_REVISION_LINCS_MAX_NEW_QUERIES",
  unset = ""
))

if (nzchar(max_new_queries_text)) {
  max_new_queries <- suppressWarnings(as.integer(max_new_queries_text))

  if (
    is.na(max_new_queries) ||
      max_new_queries < 1L ||
      as.character(max_new_queries) != max_new_queries_text
  ) {
    stop(
      "INPH_REVISION_LINCS_MAX_NEW_QUERIES must be a positive integer."
    )
  }
} else {
  max_new_queries <- Inf
}

overwrite_mismatched <- identical(
  Sys.getenv("INPH_OVERWRITE_REVISION_LINCS"),
  "1"
)

cat("Revision LINCS run mode: ", run_mode, "\n", sep = "")
cat(
  "Maximum new query checkpoints this run: ",
  if (is.infinite(max_new_queries)) "all selected" else max_new_queries,
  "\n",
  sep = ""
)

BiocParallel::register(BiocParallel::SerialParam())


# -----------------------------------------------------------------------------
# Section 1 · Validate the signature lock and all upstream input hashes
# -----------------------------------------------------------------------------

sha256_text <- function(x) {
  digest::digest(
    paste(x, collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

manifest_result <- readRDS(MANIFEST_PATH)
lock_lines <- readLines(LOCK_PATH, warn = FALSE)

required_manifest_fields <- c(
  "manifest_version",
  "manifest_status",
  "manifest_content_sha256",
  "query_order",
  "queries",
  "manifest",
  "execution_plan",
  "input_provenance"
)

missing_manifest_fields <- setdiff(
  required_manifest_fields,
  names(manifest_result)
)

if (length(missing_manifest_fields) > 0L) {
  stop(
    "Signature manifest is missing required field(s): ",
    paste(missing_manifest_fields, collapse = ", ")
  )
}

if (
  !identical(manifest_result$manifest_version, 1L) ||
    !identical(manifest_result$manifest_status, "LOCKED") ||
    length(manifest_result$query_order) != 121L ||
    length(manifest_result$queries) != 121L ||
    nrow(manifest_result$manifest) != 121L ||
    !identical(
      manifest_result$query_order,
      manifest_result$manifest$query_id
    ) ||
    !identical(
      manifest_result$query_order,
      names(manifest_result$queries)
    )
) {
  stop("Signature manifest failed its locked 121-query structure check.")
}

provenance <- as.data.frame(manifest_result$input_provenance)

if (!all(c("input_name", "path", "md5") %in% colnames(provenance))) {
  stop("Manifest input provenance is missing required columns.")
}

missing_provenance_files <- provenance$path[
  !file.exists(provenance$path)
]

if (length(missing_provenance_files) > 0L) {
  stop(
    "Locked manifest input file(s) are missing: ",
    paste(missing_provenance_files, collapse = ", ")
  )
}

current_input_md5 <- unname(tools::md5sum(provenance$path))

if (!identical(current_input_md5, as.character(provenance$md5))) {
  changed_inputs <- provenance$input_name[
    current_input_md5 != as.character(provenance$md5)
  ]
  stop(
    "Locked signature input hash changed for: ",
    paste(changed_inputs, collapse = ", "),
    ". Rebuild the signature manifest before running LINCS."
  )
}

manifest_content_lines <- c(
  "revision_signature_manifest_version=1",
  "query_count=121",
  paste0(
    "input_md5:",
    provenance$input_name,
    "=",
    provenance$md5
  ),
  paste0(
    manifest_result$manifest$query_id,
    "=",
    manifest_result$manifest$execution_sha256
  )
)

recomputed_manifest_sha256 <- sha256_text(manifest_content_lines)

if (!identical(
  recomputed_manifest_sha256,
  manifest_result$manifest_content_sha256
)) {
  stop("Recomputed manifest content hash does not match the locked manifest.")
}

expected_lock_line <- paste0(
  "manifest_content_sha256=",
  manifest_result$manifest_content_sha256
)

if (!expected_lock_line %in% lock_lines) {
  stop("Signature lock text does not match the manifest content hash.")
}

manifest_rds_md5_line <- grep(
  "^manifest_rds_md5=",
  lock_lines,
  value = TRUE
)

if (length(manifest_rds_md5_line) != 1L) {
  stop("Signature lock text has no unique manifest RDS MD5 entry.")
}

locked_manifest_rds_md5 <- sub(
  "^manifest_rds_md5=",
  "",
  manifest_rds_md5_line
)

if (!identical(
  unname(tools::md5sum(MANIFEST_PATH)),
  locked_manifest_rds_md5
)) {
  stop("Current manifest RDS MD5 does not match the signature lock file.")
}

cat(
  "Signature lock validated: ",
  manifest_result$manifest_content_sha256,
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 2 · Resolve and validate the LINCS2 reference
# -----------------------------------------------------------------------------

ext_root <- Sys.getenv("INPH_EXTERNAL_DATA")

eh_cache <- if (nzchar(ext_root)) {
  file.path(ext_root, "experiment_hub_cache")
} else {
  file.path(RAW_DIR, "experiment_hub_cache")
}

dir.create(eh_cache, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(eh_cache)) {
  stop(
    "Unable to create or access the ExperimentHub cache directory: ",
    eh_cache
  )
}

ExperimentHub::setExperimentHubOption("CACHE", eh_cache)
eh <- ExperimentHub::ExperimentHub(ask = FALSE)
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

lincs_genes <- as.character(
  rhdf5::h5read(lincs2_path, "rownames")[, 1]
)
lincs_treatments <- as.character(
  rhdf5::h5read(lincs2_path, "colnames")[, 1]
)

if (
  length(lincs_genes) != 12328L ||
    length(lincs_treatments) != 136460L ||
    anyNA(lincs_genes) ||
    anyNA(lincs_treatments) ||
    anyDuplicated(lincs_genes) ||
    anyDuplicated(lincs_treatments)
) {
  stop("LINCS2 reference dimensions or identifiers failed validation.")
}

thp1_trt_cp_treatments <- grep(
  "__THP1__trt_cp$",
  lincs_treatments,
  value = TRUE,
  ignore.case = TRUE
)

if (
  length(thp1_trt_cp_treatments) != 1376L ||
    anyDuplicated(thp1_trt_cp_treatments)
) {
  stop(
    "Expected 1,376 unique THP1 trt_cp reference profiles, observed ",
    length(thp1_trt_cp_treatments),
    "."
  )
}

cat(
  "LINCS2 reference validated: 12,328 genes x 136,460 signatures; ",
  "1,376 THP1 trt_cp profiles\n"
)


# -----------------------------------------------------------------------------
# Section 3 · Validate existing primary and Hallmark results against manifest
# -----------------------------------------------------------------------------

primary_existing <- readRDS(PRIMARY_RESULT_PATH)
hallmark_existing <- readRDS(HALLMARK_RESULT_PATH)

if (
  !inherits(primary_existing, "gessResult") ||
    !inherits(hallmark_existing, "gessResult")
) {
  stop("Existing primary or Hallmark LINCS result is not a gessResult object.")
}

extract_executed_query <- function(gess_object) {
  query <- methods::slot(gess_object, "query")
  list(
    up = unname(as.character(query$upset)),
    down = if ("downset" %in% names(query)) {
      unname(as.character(query$downset))
    } else {
      character(0)
    }
  )
}

primary_executed_query <- extract_executed_query(primary_existing)
hallmark_executed_query <- extract_executed_query(hallmark_existing)

primary_manifest_query <- manifest_result$queries[["PRIMARY_150"]]
hallmark_manifest_query <-
  manifest_result$queries[["HALLMARK_IFNA_NONOVERLAP"]]

if (
  !identical(
    primary_executed_query$up,
    primary_manifest_query$up_entrez_in_lincs
  ) ||
    !identical(
      primary_executed_query$down,
      primary_manifest_query$down_entrez_in_lincs
    )
) {
  stop("Existing primary LINCS query does not match locked PRIMARY_150.")
}

if (
  !identical(
    hallmark_executed_query$up,
    hallmark_manifest_query$up_entrez_in_lincs
  ) ||
    length(hallmark_executed_query$down) != 0L
) {
  stop("Existing Hallmark LINCS query does not match the locked manifest.")
}

primary_full_table <- as.data.frame(
  signatureSearch::result(primary_existing)
)
hallmark_full_table <- as.data.frame(
  signatureSearch::result(hallmark_existing)
)

required_result_columns <- c(
  "pert",
  "type",
  "cell",
  "NCS",
  "WTCS",
  "WTCS_Pval",
  "WTCS_FDR"
)

if (
  !all(required_result_columns %in% colnames(primary_full_table)) ||
    !all(required_result_columns %in% colnames(hallmark_full_table)) ||
    nrow(primary_full_table) != 136460L ||
    nrow(hallmark_full_table) != 136460L
) {
  stop("Existing full-reference LINCS result tables failed validation.")
}

cat("Existing PRIMARY_150 and Hallmark full-reference results validated.\n")


# -----------------------------------------------------------------------------
# Section 4 · Query execution and checkpoint helpers
# -----------------------------------------------------------------------------

checkpoint_path <- function(query_id) {
  order_id <- manifest_result$manifest$manifest_order[
    manifest_result$manifest$query_id == query_id
  ]

  if (length(order_id) != 1L || is.na(order_id)) {
    stop("No unique manifest order found for query ", query_id, ".")
  }

  file.path(
    QUERY_RESULT_DIR,
    paste0(
      sprintf("%03d", as.integer(order_id)),
      "__",
      query_id,
      ".rds"
    )
  )
}

atomic_save_rds <- function(object, path) {
  temporary_path <- tempfile(
    pattern = paste0(basename(path), ".tmp_"),
    tmpdir = dirname(path),
    fileext = ".rds"
  )

  on.exit({
    if (file.exists(temporary_path)) {
      unlink(temporary_path)
    }
  }, add = TRUE)

  saveRDS(object, temporary_path, compress = TRUE)

  if (
    !file.exists(temporary_path) ||
      is.na(file.info(temporary_path)$size) ||
      file.info(temporary_path)$size <= 0
  ) {
    stop("Failed to create temporary checkpoint for ", path, ".")
  }

  if (file.exists(path)) {
    if (!overwrite_mismatched) {
      stop(
        "Refusing to replace existing checkpoint without ",
        "INPH_OVERWRITE_REVISION_LINCS=1: ",
        path
      )
    }
    unlink(path)
  }

  if (!file.rename(temporary_path, path)) {
    stop("Failed to atomically move checkpoint into place: ", path)
  }

  invisible(path)
}

validate_result_table <- function(
    result_table,
    reference_scope,
    query_id
) {
  missing_columns <- setdiff(
    required_result_columns,
    colnames(result_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "LINCS result for ",
      query_id,
      " is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  expected_rows <- if (
    identical(reference_scope, "full_lincs2_reference")
  ) {
    136460L
  } else if (identical(reference_scope, "thp1_trt_cp_subset")) {
    1376L
  } else {
    stop("Unsupported reference scope for ", query_id, ": ", reference_scope)
  }

  if (nrow(result_table) != expected_rows) {
    stop(
      "LINCS result row count for ",
      query_id,
      " differs from expected ",
      expected_rows,
      ": observed ",
      nrow(result_table),
      "."
    )
  }

  if (identical(reference_scope, "thp1_trt_cp_subset")) {
    if (
      any(toupper(result_table$cell) != "THP1") ||
        any(result_table$type != "trt_cp")
    ) {
      stop("Subset result for ", query_id, " contains non-THP1/trt_cp rows.")
    }
  }

  invisible(TRUE)
}

run_locked_query <- function(query_id, reference_scope) {
  query <- manifest_result$queries[[query_id]]

  if (is.null(query)) {
    stop("Locked query not found: ", query_id)
  }

  qsig <- signatureSearch::qSig(
    query = list(
      upset = query$up_entrez,
      downset = query$down_entrez
    ),
    gess_method = "LINCS",
    refdb = lincs2_path
  )

  executed_query <- methods::slot(qsig, "query")
  executed_up <- unname(as.character(executed_query$upset))
  executed_down <- if ("downset" %in% names(executed_query)) {
    unname(as.character(executed_query$downset))
  } else {
    character(0)
  }

  if (
    !identical(executed_up, query$up_entrez_in_lincs) ||
      !identical(executed_down, query$down_entrez_in_lincs)
  ) {
    stop(
      "qSig() executed gene sets do not match the locked LINCS-covered ",
      "sets for ",
      query_id,
      "."
    )
  }

  ref_subset <- if (
    identical(reference_scope, "thp1_trt_cp_subset")
  ) {
    thp1_trt_cp_treatments
  } else {
    NULL
  }

  start_time <- Sys.time()

  gess_result <- signatureSearch::gess_lincs(
    qsig,
    sortby = "NCS",
    tau = FALSE,
    ref_trts = ref_subset,
    workers = 1
  )

  elapsed_seconds <- as.numeric(difftime(
    Sys.time(),
    start_time,
    units = "secs"
  ))

  result_table <- as.data.frame(
    signatureSearch::result(gess_result)
  )

  validate_result_table(
    result_table,
    reference_scope = reference_scope,
    query_id = query_id
  )

  list(
    analysis_id = "revision_03B_lincs_query_checkpoint",
    query_id = query_id,
    manifest_content_sha256 = manifest_result$manifest_content_sha256,
    signature_sha256 = query$signature_sha256,
    execution_sha256 = query$execution_sha256,
    reference_scope = reference_scope,
    tau = FALSE,
    subset_fdr_interpretable =
      identical(reference_scope, "full_lincs2_reference"),
    elapsed_seconds = elapsed_seconds,
    result_table = result_table
  )
}

checkpoint_is_valid <- function(query_id, path) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  checkpoint <- tryCatch(
    readRDS(path),
    error = function(e) NULL
  )

  if (is.null(checkpoint)) {
    if (overwrite_mismatched) {
      return(FALSE)
    }
    stop(
      "Existing checkpoint cannot be read: ",
      path,
      ". Use INPH_OVERWRITE_REVISION_LINCS=1 only if replacement is intended."
    )
  }

  query <- manifest_result$queries[[query_id]]
  metadata_row <- manifest_result$manifest[
    manifest_result$manifest$query_id == query_id,
    ,
    drop = FALSE
  ]

  valid_identity <-
    identical(checkpoint$query_id, query_id) &&
    identical(
      checkpoint$manifest_content_sha256,
      manifest_result$manifest_content_sha256
    ) &&
    identical(checkpoint$signature_sha256, query$signature_sha256) &&
    identical(checkpoint$execution_sha256, query$execution_sha256) &&
    identical(
      checkpoint$reference_scope,
      metadata_row$reference_scope
    ) &&
    identical(checkpoint$tau, FALSE)

  valid_table <- FALSE

  if (valid_identity && is.data.frame(checkpoint$result_table)) {
    valid_table <- tryCatch(
      {
        validate_result_table(
          checkpoint$result_table,
          reference_scope = metadata_row$reference_scope,
          query_id = query_id
        )
        TRUE
      },
      error = function(e) FALSE
    )
  }

  if (!valid_identity || !valid_table) {
    if (overwrite_mismatched) {
      return(FALSE)
    }
    stop(
      "Existing checkpoint does not match the locked manifest: ",
      path,
      ". Use INPH_OVERWRITE_REVISION_LINCS=1 only if replacement is intended."
    )
  }

  TRUE
}


# -----------------------------------------------------------------------------
# Section 5 · Mandatory PRIMARY_150 THP1-subset equivalence validation
# -----------------------------------------------------------------------------

validate_numeric_equivalence <- function(
    comparison,
    column,
    tolerance = 1e-12
) {
  submitted_values <- comparison[[paste0(column, "_submitted")]]
  subset_values <- comparison[[paste0(column, "_subset")]]

  if (!identical(is.na(submitted_values), is.na(subset_values))) {
    stop(column, " missing-value pattern differs in subset validation.")
  }

  comparable <- !is.na(submitted_values)
  maximum_difference <- if (any(comparable)) {
    max(abs(
      submitted_values[comparable] - subset_values[comparable]
    ))
  } else {
    0
  }

  if (!is.finite(maximum_difference) || maximum_difference > tolerance) {
    stop(
      column,
      " failed THP1-subset equivalence: maximum absolute difference ",
      format(maximum_difference, digits = 16),
      " exceeds tolerance ",
      tolerance,
      "."
    )
  }

  maximum_difference
}

run_subset_validation <- function() {
  cat("\nRunning mandatory PRIMARY_150 THP1-subset equivalence test...\n")

  subset_checkpoint <- run_locked_query(
    "PRIMARY_150",
    reference_scope = "thp1_trt_cp_subset"
  )

  submitted_thp1 <- primary_full_table |>
    dplyr::filter(
      toupper(cell) == "THP1",
      type == "trt_cp"
    ) |>
    dplyr::select(
      pert,
      cell,
      type,
      WTCS_submitted = WTCS,
      NCS_submitted = NCS,
      WTCS_Pval_submitted = WTCS_Pval
    )

  subset_thp1 <- subset_checkpoint$result_table |>
    dplyr::select(
      pert,
      cell,
      type,
      WTCS_subset = WTCS,
      NCS_subset = NCS,
      WTCS_Pval_subset = WTCS_Pval
    )

  comparison <- dplyr::inner_join(
    submitted_thp1,
    subset_thp1,
    by = c("pert", "cell", "type")
  ) |>
    dplyr::arrange(pert)

  if (
    nrow(submitted_thp1) != 1376L ||
      nrow(subset_thp1) != 1376L ||
      nrow(comparison) != 1376L ||
      anyDuplicated(submitted_thp1$pert) ||
      anyDuplicated(subset_thp1$pert)
  ) {
    stop("PRIMARY_150 subset comparison did not resolve 1,376 unique profiles.")
  }

  max_wtcs_difference <- validate_numeric_equivalence(
    comparison,
    "WTCS"
  )
  max_ncs_difference <- validate_numeric_equivalence(
    comparison,
    "NCS"
  )
  max_pvalue_difference <- validate_numeric_equivalence(
    comparison,
    "WTCS_Pval"
  )

  validation_result <- list(
    analysis_id = "revision_03B_primary_thp1_subset_validation",
    status = "PASS",
    manifest_content_sha256 = manifest_result$manifest_content_sha256,
    primary_signature_sha256 =
      primary_manifest_query$signature_sha256,
    primary_full_result_md5 = unname(
      tools::md5sum(PRIMARY_RESULT_PATH)
    ),
    thp1_profiles_n = nrow(comparison),
    tolerance = 1e-12,
    maximum_absolute_difference = c(
      WTCS = max_wtcs_difference,
      NCS = max_ncs_difference,
      WTCS_Pval = max_pvalue_difference
    ),
    elapsed_seconds = subset_checkpoint$elapsed_seconds,
    note = paste(
      "Subset WTCS/NCS/nominal p-values reproduce the full-reference",
      "THP1 rows. Subset-adjusted WTCS FDR is intentionally not compared",
      "or interpreted."
    )
  )

  atomic_save_rds(validation_result, SUBSET_VALIDATION_PATH)
  validation_result
}

subset_validation_is_valid <- function() {
  if (!file.exists(SUBSET_VALIDATION_PATH)) {
    return(FALSE)
  }

  validation <- tryCatch(
    readRDS(SUBSET_VALIDATION_PATH),
    error = function(e) NULL
  )

  valid <-
    is.list(validation) &&
    identical(validation$status, "PASS") &&
    identical(
      validation$manifest_content_sha256,
      manifest_result$manifest_content_sha256
    ) &&
    identical(
      validation$primary_signature_sha256,
      primary_manifest_query$signature_sha256
    ) &&
    identical(
      validation$primary_full_result_md5,
      unname(tools::md5sum(PRIMARY_RESULT_PATH))
    ) &&
    identical(validation$thp1_profiles_n, 1376L) &&
    all(validation$maximum_absolute_difference <= validation$tolerance)

  if (!valid && !overwrite_mismatched) {
    stop(
      "Existing THP1-subset validation does not match current locked inputs. ",
      "Use INPH_OVERWRITE_REVISION_LINCS=1 only if replacement is intended."
    )
  }

  valid
}

if (subset_validation_is_valid()) {
  subset_validation <- readRDS(SUBSET_VALIDATION_PATH)
  cat("PRIMARY_150 THP1-subset equivalence checkpoint: PASS (reused)\n")
} else {
  subset_validation <- run_subset_validation()
  cat("PRIMARY_150 THP1-subset equivalence checkpoint: PASS (new)\n")
}


# -----------------------------------------------------------------------------
# Section 6 · Run selected pending queries with per-query checkpoints
# -----------------------------------------------------------------------------

full_query_ids <- manifest_result$execution_plan$new_full_reference_queries
null_query_ids <- manifest_result$execution_plan$thp1_subset_null_queries

if (
  length(full_query_ids) != 19L ||
    length(null_query_ids) != 100L ||
    anyDuplicated(c(full_query_ids, null_query_ids))
) {
  stop("Manifest execution plan does not contain the expected 19 + 100 queries.")
}

selected_query_ids <- switch(
  run_mode,
  validate = character(0),
  full = full_query_ids,
  null = null_query_ids,
  all = c(full_query_ids, null_query_ids)
)

new_queries_completed_this_run <- 0L

for (query_id in selected_query_ids) {
  path <- checkpoint_path(query_id)

  if (checkpoint_is_valid(query_id, path)) {
    cat("Checkpoint already valid; skipping ", query_id, ".\n", sep = "")
    next
  }

  if (new_queries_completed_this_run >= max_new_queries) {
    cat(
      "New-query limit reached; remaining selected queries stay pending.\n"
    )
    break
  }

  metadata_row <- manifest_result$manifest[
    manifest_result$manifest$query_id == query_id,
    ,
    drop = FALSE
  ]

  if (nrow(metadata_row) != 1L) {
    stop("No unique manifest metadata row for query ", query_id, ".")
  }

  cat(
    "\nRunning ",
    query_id,
    " [",
    metadata_row$reference_scope,
    "]...\n",
    sep = ""
  )

  checkpoint <- run_locked_query(
    query_id,
    reference_scope = metadata_row$reference_scope
  )

  atomic_save_rds(checkpoint, path)

  if (!checkpoint_is_valid(query_id, path)) {
    stop("New checkpoint failed read-back validation for ", query_id, ".")
  }

  new_queries_completed_this_run <-
    new_queries_completed_this_run + 1L

  cat(
    "Completed ",
    query_id,
    " in ",
    round(checkpoint$elapsed_seconds / 60, 2),
    " minutes.\n",
    sep = ""
  )

  rm(checkpoint)
  invisible(gc(full = TRUE))
}


# -----------------------------------------------------------------------------
# Section 7 · Build current batch status and index
# -----------------------------------------------------------------------------

status_rows <- lapply(
  manifest_result$query_order,
  function(query_id) {
    metadata_row <- manifest_result$manifest[
      manifest_result$manifest$query_id == query_id,
      ,
      drop = FALSE
    ]

    if (metadata_row$execution_action == "reuse_validated_existing_result") {
      result_status <- "validated_existing"
      result_path <- if (query_id == "PRIMARY_150") {
        PRIMARY_RESULT_PATH
      } else if (query_id == "HALLMARK_IFNA_NONOVERLAP") {
        HALLMARK_RESULT_PATH
      } else {
        NA_character_
      }
      elapsed_seconds <- NA_real_
    } else {
      result_path <- checkpoint_path(query_id)
      valid_checkpoint <- checkpoint_is_valid(query_id, result_path)
      result_status <- if (valid_checkpoint) "complete" else "pending"

      elapsed_seconds <- if (valid_checkpoint) {
        readRDS(result_path)$elapsed_seconds
      } else {
        NA_real_
      }
    }

    data.frame(
      manifest_order = metadata_row$manifest_order,
      query_id = query_id,
      query_role = metadata_row$query_role,
      execution_action = metadata_row$execution_action,
      reference_scope = metadata_row$reference_scope,
      status = result_status,
      elapsed_seconds = elapsed_seconds,
      result_path = result_path,
      signature_sha256 = metadata_row$signature_sha256,
      execution_sha256 = metadata_row$execution_sha256,
      stringsAsFactors = FALSE
    )
  }
)

batch_status <- dplyr::bind_rows(status_rows) |>
  dplyr::arrange(manifest_order)

full_complete_n <- sum(
  batch_status$query_id %in% full_query_ids &
    batch_status$status == "complete"
)
null_complete_n <- sum(
  batch_status$query_id %in% null_query_ids &
    batch_status$status == "complete"
)

selected_queries_complete <- if (length(selected_query_ids) == 0L) {
  TRUE
} else {
  all(
    batch_status$status[
      match(selected_query_ids, batch_status$query_id)
    ] == "complete"
  )
}

all_locked_queries_complete <-
  full_complete_n == length(full_query_ids) &&
  null_complete_n == length(null_query_ids)

batch_index <- list(
  analysis_id = "revision_03B_lincs_batch_index",
  manifest_content_sha256 = manifest_result$manifest_content_sha256,
  run_mode = run_mode,
  subset_validation = subset_validation,
  new_queries_completed_this_run = new_queries_completed_this_run,
  full_reference_complete_n = full_complete_n,
  full_reference_expected_n = length(full_query_ids),
  null_complete_n = null_complete_n,
  null_expected_n = length(null_query_ids),
  selected_queries_complete = selected_queries_complete,
  all_locked_queries_complete = all_locked_queries_complete,
  status = batch_status,
  reference = list(
    experimenthub_id = "EH7297",
    genes_n = length(lincs_genes),
    signatures_n = length(lincs_treatments),
    thp1_trt_cp_profiles_n = length(thp1_trt_cp_treatments)
  )
)

saveRDS(batch_index, BATCH_INDEX_PATH)
readr::write_csv(batch_status, BATCH_STATUS_PATH)

if (
  !file.exists(BATCH_INDEX_PATH) ||
    !file.exists(BATCH_STATUS_PATH) ||
    is.na(file.info(BATCH_INDEX_PATH)$size) ||
    is.na(file.info(BATCH_STATUS_PATH)$size) ||
    file.info(BATCH_INDEX_PATH)$size <= 0 ||
    file.info(BATCH_STATUS_PATH)$size <= 0
) {
  stop("Batch index or status output was not created correctly.")
}

index_check <- readRDS(BATCH_INDEX_PATH)
status_check <- readr::read_csv(
  BATCH_STATUS_PATH,
  show_col_types = FALSE
)

if (
  !identical(
    index_check$manifest_content_sha256,
    manifest_result$manifest_content_sha256
  ) ||
    nrow(status_check) != 121L
) {
  stop("Saved batch index or status failed read-back validation.")
}


# -----------------------------------------------------------------------------
# Section 8 · Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION LINCS BATCH CHECKPOINT COMPLETE\n")
cat("============================================================\n")
cat("  Run mode                       : ", run_mode, "\n", sep = "")
cat("  THP1 subset equivalence        : PASS\n")
cat(
  "  Maximum absolute NCS difference: ",
  format(
    subset_validation$maximum_absolute_difference[["NCS"]],
    digits = 6
  ),
  "\n",
  sep = ""
)
cat(
  "  New queries completed this run: ",
  new_queries_completed_this_run,
  "\n",
  sep = ""
)
cat(
  "  Full-reference progress       : ",
  full_complete_n,
  " / ",
  length(full_query_ids),
  "\n",
  sep = ""
)
cat(
  "  THP1-null progress             : ",
  null_complete_n,
  " / ",
  length(null_query_ids),
  "\n",
  sep = ""
)
cat("============================================================\n")

cat("\nOutputs saved and verified:\n")
cat("  ", SUBSET_VALIDATION_PATH, "\n", sep = "")
cat("  ", BATCH_INDEX_PATH, "\n", sep = "")
cat("  ", BATCH_STATUS_PATH, "\n", sep = "")

if (identical(run_mode, "validate")) {
  cat(
    "\nValidation-only mode completed. No biological sensitivity or null ",
    "query was launched.\n",
    sep = ""
  )
} else if (!batch_index$selected_queries_complete) {
  cat(
    "\nThe batch is incomplete but safely checkpointed. Rerun the same mode ",
    "to resume from the next pending query.\n",
    sep = ""
  )
} else if (!batch_index$all_locked_queries_complete) {
  cat(
    "\nAll queries selected by this mode are complete. Run the remaining ",
    "query family when ready.\n",
    sep = ""
  )
} else {
  cat(
    "\nAll locked revision LINCS queries are complete. Next: build the ",
    "candidate robustness and empirical-null summaries.\n",
    sep = ""
  )
}
