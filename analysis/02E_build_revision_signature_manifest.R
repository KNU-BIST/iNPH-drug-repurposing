# =============================================================================
# Revision Analysis 02E — Build and lock the complete signature manifest
# Script: 02E_build_revision_signature_manifest.R
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
# Assemble every frozen gene-expression query required for the reviewer
# revision into one versioned manifest before any new LINCS analysis is run.
# Each query receives a deterministic SHA-256 signature hash, a defined
# execution scope, and complete symbol-to-Entrez/LINCS coverage records.
#
# Query inventory
# ---------------
# 11 non-LOO biological signatures from Analysis 02C
#  9 leave-one-donor-out signatures from revised Script 04
#  1 existing non-overlapping Hallmark IFN-alpha query
# 100 expression-matched null signatures from Analysis 02D
# ---
# 121 total locked queries
#
# Connectivity execution plan
# ---------------------------
# - PRIMARY_150: reuse the submitted full-reference result after exact query
#   and provenance validation.
# - Other biological and LOO queries: run against the complete LINCS2
#   reference with tau=FALSE; retain NCS, WTCS, and WTCS FDR.
# - Hallmark: reuse the existing full-reference result after exact query and
#   provenance validation.
# - Matched nulls: run against all THP1 trt_cp reference profiles only with
#   tau=FALSE. These results are used only for compound-specific empirical NCS
#   null distributions, not subset-adjusted FDR inference.
#
# Inputs
# ------
# data/processed/revision_signature_variants.rds
# data/processed/revision_loo_signatures.rds
# data/processed/revision_expression_matched_null_signatures.rds
# data/processed/hallmark_ifn_sensitivity_result.rds
# ExperimentHub resource EH7297
#
# Outputs
# -------
# data/processed/revision_signature_manifest.rds
# results/revision/revision_signature_manifest.csv
# results/revision/revision_signature_manifest_genes.csv
# results/revision/revision_signature_lock.txt
#
# Lock rule
# ---------
# Once this script completes and the lock file is committed, no query gene set
# or execution scope may change without rebuilding the manifest and treating
# all dependent revision LINCS results as invalid.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment, paths, and overwrite protection
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "readr",
  "digest",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "ExperimentHub",
  "rhdf5",
  "signatureSearch"
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

VARIANT_PATH <- file.path(
  PROCESSED_DIR,
  "revision_signature_variants.rds"
)
LOO_PATH <- file.path(
  PROCESSED_DIR,
  "revision_loo_signatures.rds"
)
NULL_PATH <- file.path(
  PROCESSED_DIR,
  "revision_expression_matched_null_signatures.rds"
)
HALLMARK_RESULT_PATH <- file.path(
  PROCESSED_DIR,
  "hallmark_ifn_sensitivity_result.rds"
)

OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_signature_manifest.rds"
)
OUTPUT_MANIFEST_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "revision_signature_manifest.csv"
)
OUTPUT_GENES_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "revision_signature_manifest_genes.csv"
)
OUTPUT_LOCK_TXT <- file.path(
  REVISION_RESULTS_DIR,
  "revision_signature_lock.txt"
)

INPUT_PATHS <- c(
  signature_variants = VARIANT_PATH,
  loo_signatures = LOO_PATH,
  matched_null_signatures = NULL_PATH,
  hallmark_result = HALLMARK_RESULT_PATH
)

OUTPUT_PATHS <- c(
  manifest_rds = OUTPUT_RDS,
  manifest_csv = OUTPUT_MANIFEST_CSV,
  manifest_genes_csv = OUTPUT_GENES_CSV,
  lock_txt = OUTPUT_LOCK_TXT
)

input_info <- file.info(unname(INPUT_PATHS))
bad_inputs <- names(INPUT_PATHS)[
  !file.exists(unname(INPUT_PATHS)) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(bad_inputs) > 0L) {
  stop(
    "Required manifest input(s) missing or empty: ",
    paste(bad_inputs, collapse = ", ")
  )
}

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_SIGNATURE_MANIFEST"),
  "1"
)

existing_outputs <- unname(OUTPUT_PATHS)[
  file.exists(unname(OUTPUT_PATHS))
]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "One or more signature-lock outputs already exist: ",
    paste(existing_outputs, collapse = ", "),
    ". Refusing to overwrite them automatically. If replacement is ",
    "intentional, set INPH_OVERWRITE_SIGNATURE_MANIFEST=1 for that run."
  )
}

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REVISION_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Section 1 · Load and validate upstream frozen objects
# -----------------------------------------------------------------------------

variant_result <- readRDS(VARIANT_PATH)
loo_result <- readRDS(LOO_PATH)
null_result <- readRDS(NULL_PATH)
hallmark_result <- readRDS(HALLMARK_RESULT_PATH)

validate_signature_list <- function(
    x,
    expected_n,
    object_label,
    allow_empty_down = FALSE
) {
  if (!is.list(x) || length(x) != expected_n || is.null(names(x))) {
    stop(
      object_label,
      " should contain exactly ",
      expected_n,
      " named signatures."
    )
  }

  if (anyDuplicated(names(x))) {
    stop(object_label, " contains duplicated query IDs.")
  }

  valid <- vapply(
    x,
    function(signature) {
      is.list(signature) &&
        all(c("up", "down") %in% names(signature)) &&
        length(signature$up) > 0L &&
        (allow_empty_down || length(signature$down) > 0L) &&
        !anyNA(signature$up) &&
        !anyNA(signature$down) &&
        !anyDuplicated(signature$up) &&
        !anyDuplicated(signature$down) &&
        length(intersect(signature$up, signature$down)) == 0L
    },
    logical(1)
  )

  if (!all(valid)) {
    stop(
      object_label,
      " contains one or more malformed or overlapping signatures: ",
      paste(names(x)[!valid], collapse = ", ")
    )
  }

  invisible(TRUE)
}

validate_signature_list(
  variant_result$signatures,
  expected_n = 11L,
  object_label = "Analysis 02C"
)
validate_signature_list(
  loo_result$signatures,
  expected_n = 9L,
  object_label = "Revised Script 04"
)
validate_signature_list(
  null_result$signatures,
  expected_n = 100L,
  object_label = "Analysis 02D"
)

if (!inherits(hallmark_result, "gessResult")) {
  stop(
    "hallmark_ifn_sensitivity_result.rds should contain a gessResult object."
  )
}

hallmark_query <- methods::slot(hallmark_result, "query")

if (
  !is.list(hallmark_query) ||
    !all(c("upset", "downset") %in% names(hallmark_query)) ||
    length(hallmark_query$upset) != 49L ||
    length(hallmark_query$downset) != 0L ||
    anyNA(hallmark_query$upset) ||
    anyDuplicated(hallmark_query$upset)
) {
  stop(
    "Existing Hallmark result does not contain the expected unique ",
    "49-UP/0-DOWN executed Entrez query."
  )
}

all_upstream_query_ids <- c(
  names(variant_result$signatures),
  names(loo_result$signatures),
  names(null_result$signatures),
  "HALLMARK_IFNA_NONOVERLAP"
)

if (anyDuplicated(all_upstream_query_ids)) {
  stop("Query IDs overlap across upstream signature families.")
}

cat("Frozen upstream query families validated:\n")
cat("  Analysis 02C signatures : 11\n")
cat("  Revised Script 04 LOO   : 9\n")
cat("  Analysis 02D nulls      : 100\n")
cat("  Existing Hallmark query : 1\n")


# -----------------------------------------------------------------------------
# Section 2 · Load the exact LINCS2 gene universe
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

if (
  length(lincs_genes) != 12328L ||
    anyNA(lincs_genes) ||
    anyDuplicated(lincs_genes)
) {
  stop(
    "LINCS2 gene universe differs from the expected 12,328 unique identifiers."
  )
}


# -----------------------------------------------------------------------------
# Section 3 · Define frozen query order and metadata
# -----------------------------------------------------------------------------

variant_order <- c(
  "PRIMARY_150",
  "PRIMARY_50",
  "PRIMARY_100",
  "PRIMARY_200",
  setdiff(
    names(variant_result$signatures),
    c("PRIMARY_150", "PRIMARY_50", "PRIMARY_100", "PRIMARY_200")
  )
)

if (!setequal(variant_order, names(variant_result$signatures))) {
  stop("Could not resolve the frozen Analysis 02C query order.")
}

loo_order <- names(loo_result$signatures)
null_order <- names(null_result$signatures)

query_order <- c(
  variant_order,
  loo_order,
  "HALLMARK_IFNA_NONOVERLAP",
  null_order
)

if (
  length(query_order) != 121L ||
    anyDuplicated(query_order)
) {
  stop("Complete manifest should contain 121 unique ordered query IDs.")
}

variant_metadata <- as.data.frame(variant_result$query_metadata)
loo_metadata <- as.data.frame(loo_result$query_summary)

required_variant_metadata <- c(
  "query_id",
  "family",
  "source_contrast",
  "description",
  "resolution",
  "patient6_included"
)

if (!all(required_variant_metadata %in% colnames(variant_metadata))) {
  stop("Analysis 02C query metadata is missing required fields.")
}

if (!all(c("query_id", "excluded_patient") %in% colnames(loo_metadata))) {
  stop("LOO query metadata is missing required fields.")
}


# -----------------------------------------------------------------------------
# Section 4 · Build one symbol-to-Entrez mapping dictionary
# -----------------------------------------------------------------------------

symbol_signature_sources <- c(
  variant_result$signatures[variant_order],
  loo_result$signatures[loo_order],
  null_result$signatures[null_order]
)

all_symbols <- unique(unlist(
  lapply(symbol_signature_sources, function(x) c(x$up, x$down)),
  use.names = FALSE
))

if (anyNA(all_symbols) || any(!nzchar(all_symbols))) {
  stop("Symbol-based signature collection contains missing or empty genes.")
}

mapped_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = all_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

symbol_dictionary <- data.frame(
  symbol = all_symbols,
  entrez = as.character(mapped_entrez[all_symbols]),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    in_lincs = !is.na(entrez) & entrez %in% lincs_genes
  )

if (
  nrow(symbol_dictionary) != length(all_symbols) ||
    anyDuplicated(symbol_dictionary$symbol)
) {
  stop("Symbol-to-Entrez mapping dictionary failed validation.")
}

hallmark_up_entrez <- as.character(hallmark_query$upset)
hallmark_down_entrez <- as.character(hallmark_query$downset)

hallmark_symbols_mapped <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = hallmark_up_entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

hallmark_up_symbols <- as.character(
  hallmark_symbols_mapped[hallmark_up_entrez]
)


# -----------------------------------------------------------------------------
# Section 5 · Standardize queries and create gene-level records
# -----------------------------------------------------------------------------

sha256_text <- function(x) {
  digest::digest(
    paste(x, collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

standard_queries <- list()
manifest_rows <- list()
gene_rows <- list()

add_symbol_query <- function(
    query_id,
    signature,
    family,
    source_contrast,
    description,
    resolution = NA_real_,
    excluded_patient = NA_character_,
    patient6_included = FALSE,
    query_role,
    execution_action,
    reference_scope,
    tau
) {
  up_map <- data.frame(
    direction = "UP",
    input_rank = seq_along(signature$up),
    symbol = signature$up,
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(symbol_dictionary, by = "symbol")

  down_map <- data.frame(
    direction = "DOWN",
    input_rank = seq_along(signature$down),
    symbol = signature$down,
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(symbol_dictionary, by = "symbol")

  up_entrez <- unname(up_map$entrez[!is.na(up_map$entrez)])
  down_entrez <- unname(down_map$entrez[!is.na(down_map$entrez)])

  if (
    anyDuplicated(up_entrez) ||
      anyDuplicated(down_entrez) ||
      length(intersect(up_entrez, down_entrez)) != 0L
  ) {
    stop(
      "Query ",
      query_id,
      " contains duplicated or direction-overlapping Entrez identifiers."
    )
  }

  up_lincs <- up_entrez[up_entrez %in% lincs_genes]
  down_lincs <- down_entrez[down_entrez %in% lincs_genes]

  signature_hash <- sha256_text(c(
    "input_id_type=SYMBOL",
    paste0("UP=", paste(up_entrez, collapse = "|")),
    paste0("DOWN=", paste(down_entrez, collapse = "|"))
  ))

  execution_hash <- sha256_text(c(
    signature_hash,
    paste0("execution_action=", execution_action),
    paste0("reference_scope=", reference_scope),
    paste0("tau=", tau)
  ))

  standard_queries[[query_id]] <<- list(
    input_id_type = "SYMBOL",
    up_symbols = unname(as.character(signature$up)),
    down_symbols = unname(as.character(signature$down)),
    up_entrez = up_entrez,
    down_entrez = down_entrez,
    up_entrez_in_lincs = up_lincs,
    down_entrez_in_lincs = down_lincs,
    signature_sha256 = signature_hash,
    execution_sha256 = execution_hash
  )

  manifest_rows[[query_id]] <<- data.frame(
    query_id = query_id,
    manifest_order = match(query_id, query_order),
    family = family,
    query_role = query_role,
    source_contrast = source_contrast,
    description = description,
    resolution = resolution,
    excluded_patient = excluded_patient,
    patient6_included = patient6_included,
    input_id_type = "SYMBOL",
    input_up_n = length(signature$up),
    input_down_n = length(signature$down),
    entrez_up_n = length(up_entrez),
    entrez_down_n = length(down_entrez),
    lincs_up_n = length(up_lincs),
    lincs_down_n = length(down_lincs),
    execution_action = execution_action,
    reference_scope = reference_scope,
    tau = tau,
    interpret_subset_fdr = reference_scope == "full_lincs2_reference",
    signature_sha256 = signature_hash,
    execution_sha256 = execution_hash,
    stringsAsFactors = FALSE
  )

  gene_rows[[query_id]] <<- dplyr::bind_rows(up_map, down_map) |>
    dplyr::mutate(
      query_id = query_id,
      manifest_order = match(query_id, query_order),
      included_in_entrez_query = !is.na(entrez),
      .before = 1
    )

  invisible(NULL)
}

for (query_id in variant_order) {
  metadata_row <- variant_metadata[
    variant_metadata$query_id == query_id,
    ,
    drop = FALSE
  ]

  if (nrow(metadata_row) != 1L) {
    stop("Expected one Analysis 02C metadata row for ", query_id, ".")
  }

  is_primary <- identical(query_id, "PRIMARY_150")

  add_symbol_query(
    query_id = query_id,
    signature = variant_result$signatures[[query_id]],
    family = metadata_row$family,
    source_contrast = metadata_row$source_contrast,
    description = metadata_row$description,
    resolution = metadata_row$resolution,
    patient6_included = metadata_row$patient6_included,
    query_role = if (is_primary) "primary" else "biological_sensitivity",
    execution_action = if (is_primary) {
      "reuse_validated_existing_result"
    } else {
      "run_revision_batch"
    },
    reference_scope = "full_lincs2_reference",
    tau = is_primary
  )
}

for (query_id in loo_order) {
  metadata_row <- loo_metadata[
    loo_metadata$query_id == query_id,
    ,
    drop = FALSE
  ]

  if (nrow(metadata_row) != 1L) {
    stop("Expected one LOO metadata row for ", query_id, ".")
  }

  add_symbol_query(
    query_id = query_id,
    signature = loo_result$signatures[[query_id]],
    family = "leave_one_donor_out",
    source_contrast = "PRIMARY_C5_VS_C0",
    description = paste0(
      "cluster 5 versus cluster 0 after excluding patient ",
      metadata_row$excluded_patient
    ),
    resolution = 0.5,
    excluded_patient = as.character(metadata_row$excluded_patient),
    patient6_included = FALSE,
    query_role = "donor_sensitivity",
    execution_action = "run_revision_batch",
    reference_scope = "full_lincs2_reference",
    tau = FALSE
  )
}

for (query_id in null_order) {
  add_symbol_query(
    query_id = query_id,
    signature = null_result$signatures[[query_id]],
    family = "expression_matched_null",
    source_contrast = "MATCHED_RANDOM_PRIMARY_150",
    description = paste(
      "fixed-seed expression- and mapping-coverage-matched random",
      "150-UP/150-DOWN signature"
    ),
    query_role = "specificity_null",
    execution_action = "run_revision_batch",
    reference_scope = "thp1_trt_cp_subset",
    tau = FALSE
  )
}


# -----------------------------------------------------------------------------
# Section 6 · Add the existing Entrez-native Hallmark query
# -----------------------------------------------------------------------------

hallmark_signature_hash <- sha256_text(c(
  "input_id_type=ENTREZID",
  paste0("UP=", paste(hallmark_up_entrez, collapse = "|")),
  "DOWN="
))

hallmark_execution_hash <- sha256_text(c(
  hallmark_signature_hash,
  "execution_action=reuse_validated_existing_result",
  "reference_scope=full_lincs2_reference",
  "tau=TRUE"
))

standard_queries[["HALLMARK_IFNA_NONOVERLAP"]] <- list(
  input_id_type = "ENTREZID",
  up_symbols = unname(hallmark_up_symbols),
  down_symbols = character(0),
  up_entrez = unname(hallmark_up_entrez),
  down_entrez = character(0),
  up_entrez_in_lincs = unname(
    hallmark_up_entrez[hallmark_up_entrez %in% lincs_genes]
  ),
  down_entrez_in_lincs = character(0),
  signature_sha256 = hallmark_signature_hash,
  execution_sha256 = hallmark_execution_hash
)

manifest_rows[["HALLMARK_IFNA_NONOVERLAP"]] <- data.frame(
  query_id = "HALLMARK_IFNA_NONOVERLAP",
  manifest_order = match("HALLMARK_IFNA_NONOVERLAP", query_order),
  family = "hallmark_pathway",
  query_role = "pathway_sensitivity",
  source_contrast = "MSIGDB_HALLMARK_INTERFERON_ALPHA_RESPONSE",
  description = paste(
    "49 LINCS2-represented genes retained from the 54 Hallmark IFN-alpha",
    "genes not overlapping the complete PRIMARY_150 signature; UP-only query"
  ),
  resolution = NA_real_,
  excluded_patient = NA_character_,
  patient6_included = FALSE,
  input_id_type = "ENTREZID",
  input_up_n = length(hallmark_up_entrez),
  input_down_n = 0L,
  entrez_up_n = length(hallmark_up_entrez),
  entrez_down_n = 0L,
  lincs_up_n = sum(hallmark_up_entrez %in% lincs_genes),
  lincs_down_n = 0L,
  execution_action = "reuse_validated_existing_result",
  reference_scope = "full_lincs2_reference",
  tau = TRUE,
  interpret_subset_fdr = TRUE,
  signature_sha256 = hallmark_signature_hash,
  execution_sha256 = hallmark_execution_hash,
  stringsAsFactors = FALSE
)

gene_rows[["HALLMARK_IFNA_NONOVERLAP"]] <- data.frame(
  query_id = "HALLMARK_IFNA_NONOVERLAP",
  manifest_order = match("HALLMARK_IFNA_NONOVERLAP", query_order),
  direction = "UP",
  input_rank = seq_along(hallmark_up_entrez),
  symbol = unname(hallmark_up_symbols),
  entrez = unname(hallmark_up_entrez),
  in_lincs = hallmark_up_entrez %in% lincs_genes,
  included_in_entrez_query = TRUE,
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# Section 7 · Assemble and validate the complete manifest
# -----------------------------------------------------------------------------

standard_queries <- standard_queries[query_order]
manifest <- dplyr::bind_rows(manifest_rows) |>
  dplyr::arrange(manifest_order)
manifest_genes <- dplyr::bind_rows(gene_rows) |>
  dplyr::arrange(manifest_order, direction, input_rank)

if (
  length(standard_queries) != 121L ||
    nrow(manifest) != 121L ||
    !identical(names(standard_queries), query_order) ||
    !identical(manifest$query_id, query_order) ||
    anyDuplicated(manifest$query_id) ||
    anyNA(manifest$signature_sha256) ||
    anyNA(manifest$execution_sha256)
) {
  stop("Complete revision manifest failed its 121-query identity contract.")
}

family_counts <- table(manifest$query_role)

expected_role_counts <- c(
  biological_sensitivity = 10L,
  donor_sensitivity = 9L,
  pathway_sensitivity = 1L,
  primary = 1L,
  specificity_null = 100L
)

if (!identical(
  as.integer(family_counts[names(expected_role_counts)]),
  as.integer(expected_role_counts)
)) {
  stop("Manifest query-role counts differ from the frozen 121-query plan.")
}

primary_row <- manifest[manifest$query_id == "PRIMARY_150", , drop = FALSE]

if (
  nrow(primary_row) != 1L ||
    primary_row$input_up_n != 150L ||
    primary_row$input_down_n != 150L ||
    primary_row$entrez_up_n != 143L ||
    primary_row$entrez_down_n != 137L ||
    primary_row$lincs_up_n != 114L ||
    primary_row$lincs_down_n != 101L
) {
  stop("PRIMARY_150 coverage differs from the submitted LINCS query.")
}

hallmark_row <- manifest[
  manifest$query_id == "HALLMARK_IFNA_NONOVERLAP",
  ,
  drop = FALSE
]

if (
  nrow(hallmark_row) != 1L ||
    hallmark_row$entrez_up_n != 49L ||
    hallmark_row$entrez_down_n != 0L ||
    hallmark_row$lincs_up_n != 49L
) {
  stop("Hallmark query coverage differs from the submitted sensitivity query.")
}

null_rows <- manifest[manifest$query_role == "specificity_null", , drop = FALSE]

if (
  nrow(null_rows) != 100L ||
    any(null_rows$input_up_n != 150L) ||
    any(null_rows$input_down_n != 150L) ||
    any(null_rows$entrez_up_n != 143L) ||
    any(null_rows$entrez_down_n != 137L) ||
    any(null_rows$lincs_up_n != 114L) ||
    any(null_rows$lincs_down_n != 101L)
) {
  stop("Matched null query coverage differs from PRIMARY_150.")
}

new_full_reference_rows <- manifest |>
  dplyr::filter(
    execution_action == "run_revision_batch",
    reference_scope == "full_lincs2_reference"
  )

if (nrow(new_full_reference_rows) != 19L) {
  stop("Expected exactly 19 new full-reference biological LINCS queries.")
}

if (any(
  new_full_reference_rows$lincs_up_n < 20L |
    new_full_reference_rows$lincs_down_n < 20L
)) {
  insufficient <- new_full_reference_rows$query_id[
    new_full_reference_rows$lincs_up_n < 20L |
      new_full_reference_rows$lincs_down_n < 20L
  ]
  stop(
    "New biological query has fewer than 20 LINCS-represented genes in a ",
    "direction: ",
    paste(insufficient, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# Section 8 · Create manifest-level content hash and save outputs
# -----------------------------------------------------------------------------

input_md5 <- unname(tools::md5sum(unname(INPUT_PATHS)))
names(input_md5) <- names(INPUT_PATHS)

manifest_content_lines <- c(
  "revision_signature_manifest_version=1",
  "query_count=121",
  paste0(
    "input_md5:",
    names(input_md5),
    "=",
    input_md5
  ),
  paste0(
    manifest$query_id,
    "=",
    manifest$execution_sha256
  )
)

manifest_content_sha256 <- sha256_text(manifest_content_lines)

manifest_result <- list(
  analysis_id = "revision_02E_signature_manifest",
  manifest_version = 1L,
  manifest_status = "LOCKED",
  manifest_content_sha256 = manifest_content_sha256,
  query_order = query_order,
  queries = standard_queries,
  manifest = manifest,
  manifest_genes = manifest_genes,
  execution_plan = list(
    existing_results_to_validate_and_reuse = c(
      "PRIMARY_150",
      "HALLMARK_IFNA_NONOVERLAP"
    ),
    new_full_reference_queries =
      new_full_reference_rows$query_id,
    thp1_subset_null_queries = null_order,
    full_reference_tau = FALSE,
    null_subset_tau = FALSE,
    null_inference = paste(
      "compound-specific empirical lower-tail probability from THP1 NCS;",
      "subset-adjusted WTCS FDR is not interpreted"
    ),
    mandatory_subset_validation = paste(
      "PRIMARY_150 must be rerun once against the identical THP1 trt_cp",
      "subset and reproduce the THP1 NCS values from the existing",
      "full-reference result before null batch results are accepted"
    )
  ),
  input_provenance = data.frame(
    input_name = names(INPUT_PATHS),
    path = unname(INPUT_PATHS),
    md5 = input_md5,
    stringsAsFactors = FALSE
  ),
  lincs_reference = list(
    experimenthub_id = "EH7297",
    gene_universe_n = length(lincs_genes)
  )
)

saveRDS(manifest_result, OUTPUT_RDS)
readr::write_csv(manifest, OUTPUT_MANIFEST_CSV)
readr::write_csv(manifest_genes, OUTPUT_GENES_CSV)

initial_outputs <- c(
  OUTPUT_RDS,
  OUTPUT_MANIFEST_CSV,
  OUTPUT_GENES_CSV
)

if (!all(file.exists(initial_outputs))) {
  stop("One or more signature-manifest outputs were not created.")
}

initial_output_info <- file.info(initial_outputs)

if (
  any(is.na(initial_output_info$size)) ||
    any(initial_output_info$size <= 0)
) {
  stop("One or more signature-manifest outputs are empty.")
}

output_md5 <- unname(tools::md5sum(initial_outputs))

lock_lines <- c(
  "REVISION SIGNATURE LOCK",
  "=======================",
  "status=LOCKED",
  "manifest_version=1",
  "query_count=121",
  paste0("manifest_content_sha256=", manifest_content_sha256),
  paste0("manifest_rds_md5=", output_md5[1]),
  paste0("manifest_csv_md5=", output_md5[2]),
  paste0("manifest_genes_csv_md5=", output_md5[3]),
  "",
  "QUERY COUNTS",
  "primary=1",
  "biological_sensitivity=10",
  "donor_sensitivity=9",
  "pathway_sensitivity=1",
  "specificity_null=100",
  "",
  "EXECUTION PLAN",
  "reuse_existing_full_reference_results=2",
  "new_full_reference_queries=19",
  "thp1_subset_null_queries=100",
  "",
  "LOCK RULE",
  paste(
    "Changing any query gene set, direction, mapping, reference scope, or",
    "execution setting requires a new manifest version and invalidates all",
    "dependent revision LINCS results."
  )
)

writeLines(lock_lines, OUTPUT_LOCK_TXT, useBytes = TRUE)


# -----------------------------------------------------------------------------
# Section 9 · Read-back and file-integrity validation
# -----------------------------------------------------------------------------

if (!all(file.exists(unname(OUTPUT_PATHS)))) {
  stop("One or more signature-lock outputs were not created.")
}

output_info <- file.info(unname(OUTPUT_PATHS))

if (any(is.na(output_info$size)) || any(output_info$size <= 0)) {
  stop("One or more signature-lock outputs are empty.")
}

manifest_check <- readRDS(OUTPUT_RDS)
manifest_csv_check <- readr::read_csv(
  OUTPUT_MANIFEST_CSV,
  show_col_types = FALSE
)
manifest_genes_check <- readr::read_csv(
  OUTPUT_GENES_CSV,
  show_col_types = FALSE
)
lock_check <- readLines(OUTPUT_LOCK_TXT, warn = FALSE)

if (!identical(manifest_check, manifest_result)) {
  stop("Saved manifest RDS failed exact read-back validation.")
}

if (
  nrow(manifest_csv_check) != 121L ||
    !identical(manifest_csv_check$query_id, query_order)
) {
  stop("Saved manifest CSV failed its 121-query order check.")
}

if (
  nrow(manifest_genes_check) != nrow(manifest_genes) ||
    !identical(
      sort(unique(manifest_genes_check$query_id)),
      sort(query_order)
    )
) {
  stop("Saved manifest gene table failed validation.")
}

if (!any(grepl(
  paste0("^manifest_content_sha256=", manifest_content_sha256, "$"),
  lock_check
))) {
  stop("Signature lock text does not contain the expected manifest hash.")
}


# -----------------------------------------------------------------------------
# Section 10 · Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION SIGNATURE MANIFEST LOCKED\n")
cat("============================================================\n")
cat("  Manifest version                 : 1\n")
cat("  Total locked queries             : 121\n")
cat("  Primary query                    : 1\n")
cat("  Biological sensitivity queries  : 10\n")
cat("  Leave-one-donor-out queries      : 9\n")
cat("  Hallmark pathway query           : 1\n")
cat("  Expression-matched null queries  : 100\n")
cat("  Existing results to reuse        : 2\n")
cat("  New full-reference LINCS queries : 19\n")
cat("  THP1-subset null queries          : 100\n")
cat("  Manifest SHA-256                 : ", manifest_content_sha256, "\n", sep = "")
cat("============================================================\n")

cat("\nNew biological-query LINCS coverage:\n")
print(
  new_full_reference_rows |>
    dplyr::select(
      query_id,
      entrez_up_n,
      entrez_down_n,
      lincs_up_n,
      lincs_down_n
    ),
  row.names = FALSE
)

cat("\nOutputs saved and verified:\n")
for (path in unname(OUTPUT_PATHS)) {
  cat("  ", path, "\n", sep = "")
}

cat(
  "\nSIGNATURE LOCK ACTIVE. Do not modify upstream query definitions. ",
  "Next: implement the manifest-driven LINCS batch runner.\n",
  sep = ""
)
