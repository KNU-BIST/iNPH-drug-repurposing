# =============================================================================
# Candidate robustness and empirical-null summary
# Script: 04B_candidate_robustness_and_null_summary.R
# =============================================================================

# Purpose
# -------
# Consolidate the locked revision LINCS checkpoints into reviewer-facing THP1
# robustness tables. This script performs no clustering, differential
# expression, signature generation, or LINCS query.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(signatureSearch)
})

ANALYSIS_ID <- "revision_candidate_robustness_v1"
EXPECTED_LOCK <- paste0(
  "1e3a2b827a1429e406303871d547a9dec",
  "07ee4068078c666df69ca82e5428380"
)

MANIFEST_RDS <- "data/processed/revision_signature_manifest.rds"
BATCH_INDEX_RDS <- "data/processed/revision_lincs_batch_index.rds"
PRIMARY_RESULT_RDS <- "data/processed/lincs2_gess_result.rds"
HALLMARK_RESULT_RDS <-
  "data/processed/hallmark_ifn_sensitivity_result.rds"
COMPOUND_SUMMARY_RDS <- "data/processed/lincs2_compound_summary.rds"

OUTPUT_DIR <- "results/revision"
OUTPUT_RDS <-
  "data/processed/revision_candidate_robustness_and_null_summary.rds"

required_files <- c(
  MANIFEST_RDS,
  BATCH_INDEX_RDS,
  PRIMARY_RESULT_RDS,
  HALLMARK_RESULT_RDS,
  COMPOUND_SUMMARY_RDS
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Section 1. Validate the signature lock and complete LINCS batch
# -----------------------------------------------------------------------------

manifest_object <- readRDS(MANIFEST_RDS)
batch_index <- readRDS(BATCH_INDEX_RDS)
compound_summary <- readRDS(COMPOUND_SUMMARY_RDS)

if (
  !identical(manifest_object$manifest_content_sha256, EXPECTED_LOCK) ||
    !identical(batch_index$manifest_content_sha256, EXPECTED_LOCK)
) {
  stop("The active signature manifest does not match the frozen revision lock.")
}

if (
  !isTRUE(batch_index$all_locked_queries_complete) ||
    batch_index$full_reference_complete_n != 19L ||
    batch_index$null_complete_n != 100L
) {
  stop("The complete locked LINCS batch has not been successfully finished.")
}

manifest <- manifest_object$manifest |>
  arrange(manifest_order)
status <- batch_index$status |>
  arrange(manifest_order)

if (
  nrow(manifest) != 121L ||
    anyDuplicated(manifest$query_id) ||
    anyDuplicated(status$query_id) ||
    !setequal(manifest$query_id, status$query_id)
) {
  stop("Manifest/status query accounting failed.")
}

required_annotation_columns <- c("pert", "pert_iname")
if (!all(required_annotation_columns %in% colnames(compound_summary))) {
  stop("Compound summary lacks pert and/or pert_iname.")
}

compound_annotation <- compound_summary |>
  transmute(
    pert = as.character(pert),
    pert_iname = as.character(pert_iname),
    canonical_smiles = if ("canonical_smiles" %in% colnames(compound_summary)) {
      as.character(canonical_smiles)
    } else {
      NA_character_
    },
    inchi_key = if ("inchi_key" %in% colnames(compound_summary)) {
      as.character(inchi_key)
    } else {
      NA_character_
    }
  ) |>
  distinct(pert, .keep_all = TRUE)


# -----------------------------------------------------------------------------
# Section 2. Helpers for extracting and ranking THP1 results
# -----------------------------------------------------------------------------

extract_result_table <- function(result_object) {
  if (isS4(result_object)) {
    return(as.data.frame(signatureSearch::result(result_object)))
  }

  if (is.list(result_object) && "result_table" %in% names(result_object)) {
    return(as.data.frame(result_object$result_table))
  }

  if (is.data.frame(result_object)) {
    return(as.data.frame(result_object))
  }

  stop("Unsupported LINCS result object structure.")
}

validate_result_columns <- function(result_table, query_id) {
  required <- c(
    "pert", "cell", "type", "trend", "WTCS", "WTCS_Pval",
    "WTCS_FDR", "NCS", "N_upset", "N_downset"
  )
  missing <- setdiff(required, colnames(result_table))
  if (length(missing) > 0L) {
    stop(
      "Query ", query_id, " result lacks column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

extract_thp1 <- function(result_table, query_row) {
  query_id <- query_row$query_id[[1]]
  validate_result_columns(result_table, query_id)

  thp1 <- result_table |>
    filter(
      toupper(as.character(cell)) == "THP1",
      as.character(type) == "trt_cp",
      grepl("^BRD-", as.character(pert))
    ) |>
    transmute(
      pert = as.character(pert),
      trend = as.character(trend),
      WTCS = as.numeric(WTCS),
      WTCS_Pval = as.numeric(WTCS_Pval),
      WTCS_FDR = as.numeric(WTCS_FDR),
      NCS = as.numeric(NCS),
      NCSct = if ("NCSct" %in% colnames(result_table)) {
        as.numeric(NCSct)
      } else {
        NA_real_
      },
      N_upset = as.integer(N_upset),
      N_downset = as.integer(N_downset)
    )

  if (nrow(thp1) != 1376L || anyDuplicated(thp1$pert)) {
    stop(
      query_id,
      " does not contain exactly 1,376 unique BRD THP1 trt_cp entries."
    )
  }

  n_rankable <- sum(!is.na(thp1$NCS))
  n_nonzero <- sum(!is.na(thp1$NCS) & thp1$NCS != 0)

  thp1 |>
    mutate(
      rank_all = if_else(
        !is.na(NCS),
        as.numeric(rank(NCS, ties.method = "min", na.last = "keep")),
        NA_real_
      ),
      percentile_all = if (n_rankable > 1L) {
        (rank_all - 1) / (n_rankable - 1)
      } else {
        NA_real_
      },
      rank_nonzero = if_else(
        !is.na(NCS) & NCS != 0,
        as.numeric(rank(if_else(NCS != 0, NCS, NA_real_),
                        ties.method = "min", na.last = "keep")),
        NA_real_
      ),
      percentile_nonzero = if (n_nonzero > 1L) {
        (rank_nonzero - 1) / (n_nonzero - 1)
      } else {
        NA_real_
      },
      n_rankable = n_rankable,
      n_nonzero = n_nonzero,
      query_id = query_id,
      manifest_order = query_row$manifest_order[[1]],
      family = query_row$family[[1]],
      query_role = query_row$query_role[[1]],
      source_contrast = query_row$source_contrast[[1]],
      resolution = query_row$resolution[[1]],
      excluded_patient = as.character(query_row$excluded_patient[[1]]),
      patient6_included = query_row$patient6_included[[1]],
      reference_scope = query_row$reference_scope[[1]],
      subset_fdr_interpretable = isTRUE(query_row$interpret_subset_fdr[[1]])
    ) |>
    left_join(compound_annotation, by = "pert") |>
    mutate(
      pert_iname = if_else(
        is.na(pert_iname) | pert_iname == "",
        "UNIDENTIFIED",
        pert_iname
      )
    )
}

read_query_result <- function(query_row) {
  query_id <- query_row$query_id[[1]]

  result_path <- status$result_path[match(query_id, status$query_id)]
  if (length(result_path) != 1L || is.na(result_path) || !file.exists(result_path)) {
    stop("Missing result path for ", query_id, ".")
  }

  extract_result_table(readRDS(result_path))
}


# -----------------------------------------------------------------------------
# Section 3. Consolidate primary, biological, LOO, and Hallmark THP1 results
# -----------------------------------------------------------------------------

observed_manifest <- manifest |>
  filter(query_role != "specificity_null") |>
  arrange(manifest_order)

observed_tables <- lapply(seq_len(nrow(observed_manifest)), function(i) {
  query_row <- observed_manifest[i, , drop = FALSE]
  extract_thp1(read_query_result(query_row), query_row)
})

observed_thp1 <- bind_rows(observed_tables) |>
  arrange(manifest_order, rank_all, pert)

if (
  nrow(observed_thp1) != nrow(observed_manifest) * 1376L ||
    anyDuplicated(observed_thp1[c("query_id", "pert")])
) {
  stop("Observed-query THP1 consolidation failed.")
}

# Preserve the exact submitted rank convention: ascending non-zero THP1 NCS.
primary_thp1 <- observed_thp1 |>
  filter(query_id == "PRIMARY_150")

if (
  nrow(primary_thp1) != 1376L ||
    sum(!is.na(primary_thp1$NCS) & primary_thp1$NCS != 0) != 597L
) {
  stop("PRIMARY_150 does not reproduce the submitted 597 non-zero compounds.")
}

priority_check <- primary_thp1 |>
  filter(tolower(pert_iname) %in% c("rosuvastatin", "roflumilast")) |>
  dplyr::select(pert, pert_iname, rank_nonzero, NCS, WTCS, WTCS_FDR)

expected_priority_ranks <- c(rosuvastatin = 1, roflumilast = 14)
observed_priority_ranks <- setNames(
  priority_check$rank_nonzero,
  tolower(priority_check$pert_iname)
)

if (
  nrow(priority_check) != 2L ||
    !all(observed_priority_ranks[names(expected_priority_ranks)] ==
         expected_priority_ranks)
) {
  stop("Submitted rosuvastatin/roflumilast rank reproduction failed.")
}


# -----------------------------------------------------------------------------
# Section 4. Core biological robustness and objective ranking
# -----------------------------------------------------------------------------

core_query_ids <- manifest |>
  filter(query_role %in% c("primary", "biological_sensitivity")) |>
  arrange(manifest_order) |>
  pull(query_id)

if (length(core_query_ids) != 11L || core_query_ids[[1]] != "PRIMARY_150") {
  stop("Unexpected core biological query set.")
}

core_long <- observed_thp1 |>
  filter(query_id %in% core_query_ids)

core_robustness <- core_long |>
  group_by(pert, pert_iname, canonical_smiles, inchi_key) |>
  summarise(
    core_query_n = n(),
    core_nonmissing_n = sum(!is.na(NCS)),
    core_nonzero_n = sum(!is.na(NCS) & NCS != 0),
    median_rank_all = median(rank_all, na.rm = TRUE),
    min_rank_all = min(rank_all, na.rm = TRUE),
    max_rank_all = max(rank_all, na.rm = TRUE),
    median_percentile_all = median(percentile_all, na.rm = TRUE),
    min_percentile_all = min(percentile_all, na.rm = TRUE),
    max_percentile_all = max(percentile_all, na.rm = TRUE),
    median_NCS = median(NCS, na.rm = TRUE),
    min_NCS = min(NCS, na.rm = TRUE),
    max_NCS = max(NCS, na.rm = TRUE),
    inverse_n = sum(NCS < 0, na.rm = TRUE),
    inverse_proportion = mean(NCS < 0, na.rm = TRUE),
    strong_inverse_n = sum(NCS < -0.5, na.rm = TRUE),
    strong_inverse_proportion = mean(NCS < -0.5, na.rm = TRUE),
    fdr_significant_n = sum(WTCS_FDR < 0.05, na.rm = TRUE),
    fdr_significant_proportion = mean(WTCS_FDR < 0.05, na.rm = TRUE),
    operational_hit_n = sum(NCS < -0.5 & WTCS_FDR < 0.05, na.rm = TRUE),
    operational_hit_proportion =
      mean(NCS < -0.5 & WTCS_FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    primary_thp1 |>
      transmute(
        pert,
        primary_rank_nonzero = rank_nonzero,
        primary_rank_all = rank_all,
        primary_percentile_all = percentile_all,
        primary_NCS = NCS,
        primary_WTCS = WTCS,
        primary_FDR = WTCS_FDR
      ),
    by = "pert"
  ) |>
  arrange(
    median_percentile_all,
    desc(inverse_proportion),
    primary_percentile_all,
    pert
  ) |>
  mutate(objective_robustness_rank = row_number(), .before = 1)

if (
  nrow(core_robustness) != 1376L ||
    any(core_robustness$core_query_n != 11L) ||
    anyDuplicated(core_robustness$pert)
) {
  stop("Core robustness table validation failed.")
}


# -----------------------------------------------------------------------------
# Section 5. Reviewer-specific LOO and sensitivity summaries
# -----------------------------------------------------------------------------

loo_long <- observed_thp1 |>
  filter(query_role == "donor_sensitivity")

loo_summary <- loo_long |>
  group_by(pert, pert_iname) |>
  summarise(
    loo_query_n = n(),
    loo_median_rank_all = median(rank_all, na.rm = TRUE),
    loo_min_rank_all = min(rank_all, na.rm = TRUE),
    loo_max_rank_all = max(rank_all, na.rm = TRUE),
    loo_median_NCS = median(NCS, na.rm = TRUE),
    loo_min_NCS = min(NCS, na.rm = TRUE),
    loo_max_NCS = max(NCS, na.rm = TRUE),
    loo_inverse_n = sum(NCS < 0, na.rm = TRUE),
    loo_strong_inverse_n = sum(NCS < -0.5, na.rm = TRUE),
    loo_fdr_significant_n = sum(WTCS_FDR < 0.05, na.rm = TRUE),
    loo_operational_hit_n =
      sum(NCS < -0.5 & WTCS_FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

if (any(loo_summary$loo_query_n != 9L)) {
  stop("LOO summary does not contain nine iterations per compound.")
}

signature_size_ids <- c(
  "PRIMARY_50", "PRIMARY_100", "PRIMARY_150", "PRIMARY_200"
)
signature_size_long <- observed_thp1 |>
  filter(query_id %in% signature_size_ids) |>
  mutate(
    signature_size = as.integer(sub("PRIMARY_", "", query_id)),
    query_id = factor(query_id, levels = signature_size_ids)
  ) |>
  arrange(pert, query_id)

state_definition_ids <- c(
  "PRIMARY_150",
  "C5_ALLREST_MAX150",
  "C5_MARKERLOW_MAX150",
  "MARKER_MAX150",
  "RES03_ALLREST_MAX150",
  "RES07_ALLREST_MAX150",
  "RES10_ALLREST_MAX150",
  "P6_INCLUDED_MARKER_MAX150"
)
state_definition_long <- observed_thp1 |>
  filter(query_id %in% state_definition_ids) |>
  mutate(query_id = factor(query_id, levels = state_definition_ids)) |>
  arrange(pert, query_id)


# -----------------------------------------------------------------------------
# Section 6. Compound-specific empirical null analysis
# -----------------------------------------------------------------------------

null_manifest <- manifest |>
  filter(query_role == "specificity_null") |>
  arrange(manifest_order)

null_tables <- lapply(seq_len(nrow(null_manifest)), function(i) {
  query_row <- null_manifest[i, , drop = FALSE]
  extract_thp1(read_query_result(query_row), query_row)
})

null_thp1 <- bind_rows(null_tables) |>
  arrange(manifest_order, pert)

if (
  nrow(null_thp1) != 100L * 1376L ||
    anyDuplicated(null_thp1[c("query_id", "pert")])
) {
  stop("Null-query THP1 consolidation failed.")
}

null_summary <- null_thp1 |>
  group_by(pert, pert_iname) |>
  summarise(
    null_query_n = n(),
    null_nonmissing_n = sum(!is.na(NCS)),
    null_median_NCS = median(NCS, na.rm = TRUE),
    null_mean_NCS = mean(NCS, na.rm = TRUE),
    null_sd_NCS = sd(NCS, na.rm = TRUE),
    null_q025_NCS = quantile(NCS, 0.025, na.rm = TRUE, names = FALSE),
    null_q975_NCS = quantile(NCS, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) |>
  left_join(
    primary_thp1 |>
      dplyr::select(
        pert,
        primary_NCS = NCS,
        primary_rank_nonzero = rank_nonzero,
        primary_FDR = WTCS_FDR
      ),
    by = "pert"
  )

lower_tail_counts <- null_thp1 |>
  dplyr::select(pert, null_NCS = NCS) |>
  inner_join(
    primary_thp1 |> dplyr::select(pert, primary_NCS = NCS),
    by = "pert"
  ) |>
  group_by(pert) |>
  summarise(
    null_le_observed_n = sum(null_NCS <= primary_NCS, na.rm = TRUE),
    comparable_null_n = sum(!is.na(null_NCS) & !is.na(primary_NCS)),
    .groups = "drop"
  )

null_summary <- null_summary |>
  left_join(lower_tail_counts, by = "pert") |>
  mutate(
    observed_minus_null_median = primary_NCS - null_median_NCS,
    empirical_lower_tail_p =
      (1 + null_le_observed_n) / (1 + comparable_null_n),
    empirical_percentile = null_le_observed_n / comparable_null_n
  ) |>
  arrange(empirical_lower_tail_p, primary_NCS, pert)

if (
  nrow(null_summary) != 1376L ||
    any(null_summary$null_query_n != 100L) ||
    any(null_summary$comparable_null_n != 100L) ||
    any(null_summary$empirical_lower_tail_p < 1 / 101 - 1e-12)
) {
  stop("Empirical-null summary validation failed.")
}


# -----------------------------------------------------------------------------
# Section 7. Integrated candidate table and focused inspection set
# -----------------------------------------------------------------------------

candidate_summary <- core_robustness |>
  left_join(loo_summary |> dplyr::select(-pert_iname), by = "pert") |>
  left_join(
    null_summary |>
      dplyr::select(
        pert,
        null_median_NCS,
        null_q025_NCS,
        null_q975_NCS,
        observed_minus_null_median,
        empirical_lower_tail_p,
        empirical_percentile
      ),
    by = "pert"
  ) |>
  arrange(objective_robustness_rank)

inspection_perts <- unique(c(
  candidate_summary$pert[
    seq_len(min(30L, nrow(candidate_summary)))
  ],
  priority_check$pert
))

priority_query_results <- observed_thp1 |>
  filter(pert %in% inspection_perts) |>
  dplyr::select(
    pert,
    pert_iname,
    query_id,
    family,
    query_role,
    excluded_patient,
    rank_nonzero,
    rank_all,
    percentile_all,
    NCS,
    WTCS,
    WTCS_FDR,
    N_upset,
    N_downset
  )

priority_candidate_summary <- candidate_summary |>
  filter(pert %in% inspection_perts)


# -----------------------------------------------------------------------------
# Section 8. Save machine-readable outputs and verify them
# -----------------------------------------------------------------------------

output_paths <- c(
  observed_long = file.path(
    OUTPUT_DIR, "revision_lincs_thp1_observed_query_results.csv"
  ),
  core_robustness = file.path(
    OUTPUT_DIR, "revision_candidate_core_robustness.csv"
  ),
  loo_long = file.path(
    OUTPUT_DIR, "revision_loo_lincs_results.csv"
  ),
  loo_summary = file.path(
    OUTPUT_DIR, "revision_loo_candidate_summary.csv"
  ),
  signature_size = file.path(
    OUTPUT_DIR, "revision_signature_size_lincs_results.csv"
  ),
  state_definition = file.path(
    OUTPUT_DIR, "revision_state_definition_lincs_results.csv"
  ),
  null_summary = file.path(
    OUTPUT_DIR, "revision_empirical_null_summary.csv"
  ),
  candidate_summary = file.path(
    OUTPUT_DIR, "revision_integrated_candidate_summary.csv"
  ),
  priority_queries = file.path(
    OUTPUT_DIR, "revision_priority_candidate_query_results.csv"
  ),
  priority_summary = file.path(
    OUTPUT_DIR, "revision_priority_candidate_summary.csv"
  )
)

write_csv(observed_thp1, output_paths[["observed_long"]], na = "")
write_csv(core_robustness, output_paths[["core_robustness"]], na = "")
write_csv(loo_long, output_paths[["loo_long"]], na = "")
write_csv(loo_summary, output_paths[["loo_summary"]], na = "")
write_csv(signature_size_long, output_paths[["signature_size"]], na = "")
write_csv(state_definition_long, output_paths[["state_definition"]], na = "")
write_csv(null_summary, output_paths[["null_summary"]], na = "")
write_csv(candidate_summary, output_paths[["candidate_summary"]], na = "")
write_csv(priority_query_results, output_paths[["priority_queries"]], na = "")
write_csv(priority_candidate_summary, output_paths[["priority_summary"]], na = "")

result <- list(
  analysis_id = ANALYSIS_ID,
  manifest_content_sha256 = EXPECTED_LOCK,
  ranking_definitions = list(
    rank_nonzero =
      "Ascending THP1 NCS among non-zero compounds; reproduces submitted rank",
    percentile_all = paste(
      "Ascending THP1 NCS rank among all 1,376 THP1 trt_cp compounds,",
      "scaled from 0 (most negative) to 1 (most positive)"
    ),
    objective_ranking = paste(
      "Median percentile_all across the 11 primary/biological queries;",
      "inverse-direction consistency second; primary percentile third"
    ),
    operational_hit = "THP1 NCS < -0.5 and WTCS FDR < 0.05",
    empirical_p = paste(
      "(1 + number of matched-null NCS <= PRIMARY_150 NCS) / 101;",
      "compound-specific lower-tail test"
    )
  ),
  core_query_ids = core_query_ids,
  primary_priority_check = priority_check,
  observed_thp1 = observed_thp1,
  core_robustness = core_robustness,
  loo_summary = loo_summary,
  null_summary = null_summary,
  candidate_summary = candidate_summary,
  priority_candidate_summary = priority_candidate_summary,
  output_paths = output_paths,
  input_provenance = data.frame(
    input = required_files,
    md5 = unname(tools::md5sum(required_files)),
    stringsAsFactors = FALSE
  )
)

saveRDS(result, OUTPUT_RDS, compress = TRUE)

if (
  !file.exists(OUTPUT_RDS) ||
    any(!file.exists(output_paths)) ||
    any(file.info(output_paths)$size <= 0)
) {
  stop("One or more candidate-summary outputs were not created correctly.")
}

result_check <- readRDS(OUTPUT_RDS)
if (
  !identical(result_check$manifest_content_sha256, EXPECTED_LOCK) ||
    nrow(result_check$candidate_summary) != 1376L ||
    nrow(result_check$null_summary) != 1376L
) {
  stop("Post-write candidate-summary validation failed.")
}

cat("\nSubmitted priority-compound reproduction:\n")
print(priority_check, row.names = FALSE)

cat("\nTop 20 candidates by predefined robustness ranking:\n")
print(
  candidate_summary |>
    dplyr::select(
      objective_robustness_rank,
      pert,
      pert_iname,
      primary_rank_nonzero,
      primary_NCS,
      median_percentile_all,
      inverse_proportion,
      operational_hit_n,
      loo_operational_hit_n,
      empirical_lower_tail_p
    ) |>
    slice_head(n = 20L),
  row.names = FALSE
)

cat("\nRosuvastatin and roflumilast integrated results:\n")
print(
  candidate_summary |>
    filter(tolower(pert_iname) %in% c("rosuvastatin", "roflumilast")) |>
    dplyr::select(
      objective_robustness_rank,
      pert,
      pert_iname,
      primary_rank_nonzero,
      primary_NCS,
      primary_FDR,
      median_percentile_all,
      inverse_n,
      strong_inverse_n,
      operational_hit_n,
      loo_median_rank_all,
      loo_operational_hit_n,
      null_median_NCS,
      empirical_lower_tail_p
    ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 04B COMPLETE\n")
cat("============================================================\n")
cat("  Observed biological/pathway queries :", nrow(observed_manifest), "\n")
cat("  Core robustness queries             :", length(core_query_ids), "\n")
cat("  LOO queries                         :", 9L, "\n")
cat("  Matched-null queries                :", 100L, "\n")
cat("  THP1 compounds per query            :", 1376L, "\n")
cat("  Submitted non-zero primary compounds:", 597L, "\n")
cat("  Output RDS                          :", OUTPUT_RDS, "\n")
cat("============================================================\n")
cat("Next: inspect the exported tables and generate diagnostic figures.\n")
