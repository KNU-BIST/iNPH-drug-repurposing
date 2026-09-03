# =============================================================================
# Revision Analysis 02D — Expression-matched null signatures
# Script: 02D_generate_expression_matched_null_signatures.R
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
# Generate a frozen set of random 150-UP/150-DOWN signatures for the reviewer-
# requested LINCS specificity analysis. Each null signature:
#
# 1. excludes every gene in the submitted PRIMARY_150 signature;
# 2. contains 150 unique UP and 150 unique DOWN genes with no overlap;
# 3. matches the primary signature separately by direction for baseline
#    expression (DESeq2 baseMean); and
# 4. exactly matches the primary signature's three gene-mapping classes:
#    represented in LINCS2, Entrez-mapped but absent from LINCS2, and unmapped.
#
# The script generates signatures only. It does not run LINCS. The frozen null
# signatures will be added to the revision signature manifest before any new
# connectivity searches are performed.
#
# Inputs
# ------
# data/processed/inph_ifn_signature.rds
# data/processed/de_ifn_subtype_results.rds
# ExperimentHub resource EH7297
#
# Outputs
# -------
# data/processed/revision_expression_matched_null_signatures.rds
# results/revision/expression_matched_null_signature_summary.csv
# results/revision/expression_matched_null_signature_genes.csv
#
# Frozen revision-stage choices
# -----------------------------
# Null signatures             : 100
# Genes per direction         : 150
# Random seed                 : 20260903
# Nearest-expression pool     : 25 candidates per target slot
# Expression distance         : absolute difference in log1p(baseMean)
# Candidate reuse             : prohibited within a query; allowed across nulls
# Primary-signature genes     : excluded from every null query
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment, configuration, and overwrite protection
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "readr",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "ExperimentHub",
  "rhdf5"
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

PRIMARY_SIGNATURE_PATH <- file.path(
  PROCESSED_DIR,
  "inph_ifn_signature.rds"
)
PRIMARY_DE_PATH <- file.path(
  PROCESSED_DIR,
  "de_ifn_subtype_results.rds"
)

OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_expression_matched_null_signatures.rds"
)
OUTPUT_SUMMARY_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "expression_matched_null_signature_summary.csv"
)
OUTPUT_GENES_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "expression_matched_null_signature_genes.csv"
)

OUTPUT_PATHS <- c(
  OUTPUT_RDS,
  OUTPUT_SUMMARY_CSV,
  OUTPUT_GENES_CSV
)

N_NULL_SIGNATURES <- 100L
SIGNATURE_SIZE <- 150L
RANDOM_SEED <- 20260903L
NEAREST_POOL_SIZE <- 25L

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_NULL_SIGNATURES"),
  "1"
)

existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "One or more null-signature outputs already exist: ",
    paste(existing_outputs, collapse = ", "),
    ". Refusing to overwrite them automatically. If replacement is ",
    "intentional, set INPH_OVERWRITE_NULL_SIGNATURES=1 for that run."
  )
}

required_inputs <- c(
  PRIMARY_SIGNATURE_PATH,
  PRIMARY_DE_PATH
)

input_info <- file.info(required_inputs)
bad_inputs <- required_inputs[
  !file.exists(required_inputs) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(bad_inputs) > 0L) {
  stop(
    "Required input file(s) missing or empty: ",
    paste(bad_inputs, collapse = ", ")
  )
}

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REVISION_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Section 1 · Load and validate the submitted primary signature and DE table
# -----------------------------------------------------------------------------

primary_signature <- readRDS(PRIMARY_SIGNATURE_PATH)
primary_de <- as.data.frame(readRDS(PRIMARY_DE_PATH))

if (
  !is.list(primary_signature) ||
    !all(c("up", "down") %in% names(primary_signature))
) {
  stop(
    "inph_ifn_signature.rds must contain a list with 'up' and 'down' elements."
  )
}

if (
  length(primary_signature$up) != SIGNATURE_SIZE ||
    length(primary_signature$down) != SIGNATURE_SIZE ||
    anyDuplicated(primary_signature$up) ||
    anyDuplicated(primary_signature$down) ||
    length(intersect(primary_signature$up, primary_signature$down)) != 0L
) {
  stop(
    "Submitted primary signature failed the expected disjoint ",
    "150-UP/150-DOWN validation."
  )
}

required_de_columns <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "padj"
)

missing_de_columns <- setdiff(
  required_de_columns,
  colnames(primary_de)
)

if (length(missing_de_columns) > 0L) {
  stop(
    "Primary DE table is missing required column(s): ",
    paste(missing_de_columns, collapse = ", ")
  )
}

primary_de <- primary_de |>
  dplyr::filter(
    !is.na(gene),
    nzchar(gene),
    !is.na(baseMean),
    is.finite(baseMean),
    baseMean >= 0
  ) |>
  dplyr::distinct(gene, .keep_all = TRUE)

missing_primary_genes <- setdiff(
  c(primary_signature$up, primary_signature$down),
  primary_de$gene
)

if (length(missing_primary_genes) > 0L) {
  stop(
    "Primary signature gene(s) missing from the DE table: ",
    paste(missing_primary_genes, collapse = ", ")
  )
}

cat(
  "Primary signature validated: ",
  length(primary_signature$up),
  " UP / ",
  length(primary_signature$down),
  " DOWN genes\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 2 · Resolve LINCS2 gene coverage
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
  length(lincs_genes) == 0L ||
    anyNA(lincs_genes) ||
    anyDuplicated(lincs_genes)
) {
  stop("LINCS2 row-name identifiers failed validation.")
}

cat(
  "LINCS2 gene universe loaded: ",
  length(lincs_genes),
  " identifiers\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 3 · Annotate mapping class and baseline expression
# -----------------------------------------------------------------------------

mapped_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = primary_de$gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_universe <- primary_de |>
  dplyr::mutate(
    entrez = as.character(mapped_entrez[gene]),
    mapping_class = dplyr::case_when(
      is.na(entrez) ~ "unmapped",
      entrez %in% lincs_genes ~ "lincs_represented",
      TRUE ~ "mapped_not_lincs"
    ),
    log1p_base_mean = log1p(baseMean)
  ) |>
  dplyr::arrange(gene)

primary_gene_set <- unique(c(
  primary_signature$up,
  primary_signature$down
))

candidate_universe <- gene_universe |>
  dplyr::filter(!gene %in% primary_gene_set)

build_target_profile <- function(genes, direction) {
  target <- data.frame(
    gene = genes,
    target_rank = seq_along(genes),
    direction = direction,
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      gene_universe |>
        dplyr::select(
          gene,
          baseMean,
          log1p_base_mean,
          entrez,
          mapping_class
        ),
      by = "gene"
    )

  if (
    nrow(target) != SIGNATURE_SIZE ||
      anyNA(target$baseMean) ||
      anyNA(target$mapping_class)
  ) {
    stop("Failed to construct the primary ", direction, " target profile.")
  }

  target
}

target_up <- build_target_profile(primary_signature$up, "UP")
target_down <- build_target_profile(primary_signature$down, "DOWN")

mapping_profile <- dplyr::bind_rows(target_up, target_down) |>
  dplyr::count(direction, mapping_class, name = "n") |>
  dplyr::arrange(direction, mapping_class) |>
  as.data.frame()

expected_mapping_profile <- tibble::tribble(
  ~direction, ~mapping_class,       ~n,
  "DOWN",     "lincs_represented", 101L,
  "DOWN",     "mapped_not_lincs",   36L,
  "DOWN",     "unmapped",           13L,
  "UP",       "lincs_represented", 114L,
  "UP",       "mapped_not_lincs",   29L,
  "UP",       "unmapped",            7L
) |>
  dplyr::arrange(direction, mapping_class) |>
  as.data.frame()

if (!identical(mapping_profile, expected_mapping_profile)) {
  stop(
    "Primary signature mapping profile differs from the submitted LINCS ",
    "coverage: expected UP 114/29/7 and DOWN 101/36/13 across ",
    "LINCS-represented/mapped-not-LINCS/unmapped classes."
  )
}

candidate_class_counts <- candidate_universe |>
  dplyr::count(mapping_class, name = "candidate_n")

cat("Primary direction-specific mapping profile:\n")
print(mapping_profile)
cat("\nEligible background genes by mapping class:\n")
print(candidate_class_counts)


# -----------------------------------------------------------------------------
# Section 4 · Deterministic nearest-expression matching
# -----------------------------------------------------------------------------

match_direction <- function(
    target,
    candidates,
    excluded_genes = character(0),
    nearest_pool_size = NEAREST_POOL_SIZE
) {
  selected_genes <- character(0)
  target_order <- sample(seq_len(nrow(target)))
  match_rows <- vector("list", nrow(target))

  for (target_index in target_order) {
    target_row <- target[target_index, , drop = FALSE]

    available <- candidates |>
      dplyr::filter(
        mapping_class == target_row$mapping_class,
        !gene %in% excluded_genes,
        !gene %in% selected_genes
      ) |>
      dplyr::mutate(
        expression_distance = abs(
          log1p_base_mean - target_row$log1p_base_mean
        )
      ) |>
      dplyr::arrange(expression_distance, gene)

    if (nrow(available) == 0L) {
      stop(
        "No candidate remained for ",
        target_row$direction,
        " target gene ",
        target_row$gene,
        " in mapping class ",
        target_row$mapping_class,
        "."
      )
    }

    nearest <- available |>
      dplyr::slice_head(
        n = min(nearest_pool_size, nrow(available))
      )

    selected_index <- sample.int(nrow(nearest), size = 1L)
    selected <- nearest[selected_index, , drop = FALSE]
    selected_genes <- c(selected_genes, selected$gene)

    match_rows[[target_index]] <- data.frame(
      direction = target_row$direction,
      target_rank = target_row$target_rank,
      target_gene = target_row$gene,
      target_base_mean = target_row$baseMean,
      target_log1p_base_mean = target_row$log1p_base_mean,
      mapping_class = target_row$mapping_class,
      matched_gene = selected$gene,
      matched_base_mean = selected$baseMean,
      matched_log1p_base_mean = selected$log1p_base_mean,
      matched_entrez = selected$entrez,
      absolute_log1p_expression_difference =
        selected$expression_distance,
      stringsAsFactors = FALSE
    )
  }

  matches <- dplyr::bind_rows(match_rows) |>
    dplyr::arrange(target_rank)

  if (
    nrow(matches) != nrow(target) ||
      anyDuplicated(matches$matched_gene) ||
      any(matches$matched_gene %in% primary_gene_set) ||
      !identical(matches$mapping_class, target$mapping_class)
  ) {
    stop("Nearest-expression matching failed its direction-level contract.")
  }

  matches
}

set.seed(RANDOM_SEED)

null_signatures <- vector("list", N_NULL_SIGNATURES)
names(null_signatures) <- sprintf(
  "NULL_MATCHED_%03d",
  seq_len(N_NULL_SIGNATURES)
)

null_gene_rows <- vector("list", N_NULL_SIGNATURES)

for (null_index in seq_len(N_NULL_SIGNATURES)) {
  query_id <- names(null_signatures)[null_index]

  up_matches <- match_direction(
    target = target_up,
    candidates = candidate_universe
  )

  down_matches <- match_direction(
    target = target_down,
    candidates = candidate_universe,
    excluded_genes = up_matches$matched_gene
  )

  null_signatures[[query_id]] <- list(
    up = up_matches$matched_gene,
    down = down_matches$matched_gene
  )

  null_gene_rows[[null_index]] <- dplyr::bind_rows(
    up_matches,
    down_matches
  ) |>
    dplyr::mutate(
      query_id = query_id,
      .before = 1
    )
}

null_genes <- dplyr::bind_rows(null_gene_rows)


# -----------------------------------------------------------------------------
# Section 5 · Validate every null query and summarize matching quality
# -----------------------------------------------------------------------------

validate_null_signature <- function(query_id) {
  signature <- null_signatures[[query_id]]
  query_genes <- null_genes |>
    dplyr::filter(query_id == .env$query_id)

  if (
    length(signature$up) != SIGNATURE_SIZE ||
      length(signature$down) != SIGNATURE_SIZE ||
      anyDuplicated(signature$up) ||
      anyDuplicated(signature$down) ||
      length(intersect(signature$up, signature$down)) != 0L ||
      any(c(signature$up, signature$down) %in% primary_gene_set)
  ) {
    stop("Null signature ", query_id, " failed its 150/150 disjointness contract.")
  }

  observed_mapping_profile <- query_genes |>
    dplyr::count(direction, mapping_class, name = "n") |>
    dplyr::arrange(direction, mapping_class) |>
    as.data.frame()

  if (!identical(observed_mapping_profile, expected_mapping_profile)) {
    stop("Null signature ", query_id, " failed exact mapping-profile matching.")
  }

  up_rows <- query_genes |>
    dplyr::filter(direction == "UP")
  down_rows <- query_genes |>
    dplyr::filter(direction == "DOWN")

  data.frame(
    query_id = query_id,
    family = "expression_matched_null",
    signature_up_n = length(signature$up),
    signature_down_n = length(signature$down),
    up_entrez_mapped_n = sum(!is.na(up_rows$matched_entrez)),
    down_entrez_mapped_n = sum(!is.na(down_rows$matched_entrez)),
    up_lincs_represented_n =
      sum(up_rows$mapping_class == "lincs_represented"),
    down_lincs_represented_n =
      sum(down_rows$mapping_class == "lincs_represented"),
    up_median_abs_log1p_expression_difference = median(
      up_rows$absolute_log1p_expression_difference
    ),
    down_median_abs_log1p_expression_difference = median(
      down_rows$absolute_log1p_expression_difference
    ),
    up_p95_abs_log1p_expression_difference = unname(stats::quantile(
      up_rows$absolute_log1p_expression_difference,
      probs = 0.95,
      names = FALSE,
      type = 7
    )),
    down_p95_abs_log1p_expression_difference = unname(stats::quantile(
      down_rows$absolute_log1p_expression_difference,
      probs = 0.95,
      names = FALSE,
      type = 7
    )),
    overlap_primary_n = length(intersect(
      c(signature$up, signature$down),
      primary_gene_set
    )),
    stringsAsFactors = FALSE
  )
}

null_summary <- dplyr::bind_rows(lapply(
  names(null_signatures),
  validate_null_signature
))

if (
  nrow(null_summary) != N_NULL_SIGNATURES ||
    any(null_summary$signature_up_n != SIGNATURE_SIZE) ||
    any(null_summary$signature_down_n != SIGNATURE_SIZE) ||
    any(null_summary$up_entrez_mapped_n != 143L) ||
    any(null_summary$down_entrez_mapped_n != 137L) ||
    any(null_summary$up_lincs_represented_n != 114L) ||
    any(null_summary$down_lincs_represented_n != 101L) ||
    any(null_summary$overlap_primary_n != 0L)
) {
  stop("One or more null signatures failed the frozen design contract.")
}

all_null_up <- unlist(
  lapply(null_signatures, `[[`, "up"),
  use.names = FALSE
)
all_null_down <- unlist(
  lapply(null_signatures, `[[`, "down"),
  use.names = FALSE
)

null_diversity <- list(
  unique_up_genes = length(unique(all_null_up)),
  unique_down_genes = length(unique(all_null_down)),
  unique_genes_overall = length(unique(c(all_null_up, all_null_down))),
  maximum_up_gene_frequency = max(table(all_null_up)),
  maximum_down_gene_frequency = max(table(all_null_down))
)


# -----------------------------------------------------------------------------
# Section 6 · Save and read back the frozen null-signature set
# -----------------------------------------------------------------------------

result <- list(
  analysis_id = "revision_02D_expression_matched_null_signatures",
  configuration = list(
    n_null_signatures = N_NULL_SIGNATURES,
    signature_size_per_direction = SIGNATURE_SIZE,
    random_seed = RANDOM_SEED,
    nearest_expression_pool_size = NEAREST_POOL_SIZE,
    expression_measure = "DESeq2 baseMean",
    matching_distance = "absolute difference in log1p(baseMean)",
    mapping_classes = c(
      "lincs_represented",
      "mapped_not_lincs",
      "unmapped"
    ),
    candidate_reuse_within_query = FALSE,
    candidate_reuse_across_queries = TRUE,
    primary_signature_genes_excluded = TRUE,
    intended_lincs_scope = paste(
      "All THP1 trt_cp reference profiles only; tau=FALSE;",
      "compound-specific empirical null inference uses NCS"
    )
  ),
  signatures = null_signatures,
  query_summary = null_summary,
  gene_matches = null_genes,
  primary_mapping_profile = mapping_profile,
  candidate_mapping_class_counts = candidate_class_counts,
  diversity = null_diversity,
  input_provenance = list(
    primary_signature_md5 = unname(
      tools::md5sum(PRIMARY_SIGNATURE_PATH)
    ),
    primary_de_md5 = unname(
      tools::md5sum(PRIMARY_DE_PATH)
    ),
    lincs_experimenthub_id = "EH7297",
    lincs_gene_universe_n = length(lincs_genes)
  )
)

saveRDS(result, OUTPUT_RDS)
readr::write_csv(null_summary, OUTPUT_SUMMARY_CSV)
readr::write_csv(null_genes, OUTPUT_GENES_CSV)

if (!all(file.exists(OUTPUT_PATHS))) {
  stop("One or more null-signature outputs were not created.")
}

output_info <- file.info(OUTPUT_PATHS)

if (any(is.na(output_info$size)) || any(output_info$size <= 0)) {
  stop("One or more null-signature outputs are empty.")
}

result_check <- readRDS(OUTPUT_RDS)

if (!identical(result_check, result)) {
  stop("Saved null-signature RDS failed exact read-back validation.")
}

summary_check <- readr::read_csv(
  OUTPUT_SUMMARY_CSV,
  show_col_types = FALSE
)
genes_check <- readr::read_csv(
  OUTPUT_GENES_CSV,
  show_col_types = FALSE
)

if (nrow(summary_check) != N_NULL_SIGNATURES) {
  stop("Saved null-signature summary should contain 100 rows.")
}

if (
  nrow(genes_check) !=
    N_NULL_SIGNATURES * SIGNATURE_SIZE * 2L
) {
  stop("Saved null-signature gene table should contain 30,000 rows.")
}


# -----------------------------------------------------------------------------
# Section 7 · Console report
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 02D COMPLETE\n")
cat("============================================================\n")
cat("  Null signatures                 : ", N_NULL_SIGNATURES, "\n", sep = "")
cat("  Genes per direction             : ", SIGNATURE_SIZE, "\n", sep = "")
cat("  Random seed                     : ", RANDOM_SEED, "\n", sep = "")
cat("  Nearest-expression pool         : ", NEAREST_POOL_SIZE, "\n", sep = "")
cat("  UP Entrez/LINCS coverage        : 143 / 114\n")
cat("  DOWN Entrez/LINCS coverage      : 137 / 101\n")
cat(
  "  Median UP matching distance     : ",
  format(
    median(null_summary$up_median_abs_log1p_expression_difference),
    digits = 4
  ),
  "\n",
  sep = ""
)
cat(
  "  Median DOWN matching distance   : ",
  format(
    median(null_summary$down_median_abs_log1p_expression_difference),
    digits = 4
  ),
  "\n",
  sep = ""
)
cat(
  "  Unique null genes overall       : ",
  null_diversity$unique_genes_overall,
  "\n",
  sep = ""
)
cat("  Primary genes in null queries   : 0\n")
cat("============================================================\n")

cat("\nOutputs saved and verified:\n")
for (path in OUTPUT_PATHS) {
  cat("  ", path, "\n", sep = "")
}

cat(
  "\nNext: build and hash the complete revision signature manifest. ",
  "Do not run LINCS before the manifest is locked.\n",
  sep = ""
)
