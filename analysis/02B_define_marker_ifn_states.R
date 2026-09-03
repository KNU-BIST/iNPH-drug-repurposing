# =============================================================================
# 02B_define_marker_ifn_states.R
# =============================================================================
#
# Purpose
# -------
# Freeze a clustering-independent IFN-high/IFN-low definition for the reviewer
# revision. The canonical nine-patient arm reuses the exact cell-level AUCell
# scores retained by 02_cross_disease_contextualization.R. The uniform
# ten-patient arm applies the same 15-gene selector, aucMaxRank, and externally
# derived healthy-control threshold to the population reconstructed by
# 01C_patient6_sensitivity.R.
#
# This script defines cell states only. It performs no differential-expression
# analysis and no LINCS query.
#
# Inputs
# ------
# data/processed/phase2_contextualization.rds
# data/processed/revision_patient6_uniform_monocytes.rds
# data/processed/revision_patient6_sensitivity.rds
#
# Outputs
# -------
# data/processed/revision_marker_ifn_states.rds
# results/revision/marker_ifn_state_population_summary.csv
# results/revision/marker_ifn_state_donor_counts.csv
# results/revision/marker_ifn_state_cluster_concordance.csv
# results/revision/marker_ifn_state_scoring_concordance.csv
# =============================================================================

# -----------------------------------------------------------------------------
# Section 1 · Packages, paths, and overwrite protection
# -----------------------------------------------------------------------------

required_packages <- c(
  "AUCell",
  "Matrix",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "dplyr",
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

PHASE2_PATH <- file.path(
  PROCESSED_DIR,
  "phase2_contextualization.rds"
)
UNIFORM_MONOCYTES_PATH <- file.path(
  PROCESSED_DIR,
  "revision_patient6_uniform_monocytes.rds"
)
PATIENT6_RESULT_PATH <- file.path(
  PROCESSED_DIR,
  "revision_patient6_sensitivity.rds"
)

OUTPUT_RDS <- file.path(
  PROCESSED_DIR,
  "revision_marker_ifn_states.rds"
)
OUTPUT_POPULATION_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "marker_ifn_state_population_summary.csv"
)
OUTPUT_DONOR_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "marker_ifn_state_donor_counts.csv"
)
OUTPUT_CLUSTER_CONCORDANCE_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "marker_ifn_state_cluster_concordance.csv"
)
OUTPUT_SCORING_CONCORDANCE_CSV <- file.path(
  REVISION_RESULTS_DIR,
  "marker_ifn_state_scoring_concordance.csv"
)

OUTPUT_PATHS <- c(
  OUTPUT_RDS,
  OUTPUT_POPULATION_CSV,
  OUTPUT_DONOR_CSV,
  OUTPUT_CLUSTER_CONCORDANCE_CSV,
  OUTPUT_SCORING_CONCORDANCE_CSV
)
required_inputs <- c(
  PHASE2_PATH,
  UNIFORM_MONOCYTES_PATH,
  PATIENT6_RESULT_PATH
)

if (!all(file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

overwrite_allowed <- identical(
  Sys.getenv("INPH_OVERWRITE_MARKER_IFN_STATES"),
  "1"
)
existing_outputs <- OUTPUT_PATHS[file.exists(OUTPUT_PATHS)]

if (length(existing_outputs) > 0L && !overwrite_allowed) {
  stop(
    "Marker-state output(s) already exist: ",
    paste(existing_outputs, collapse = ", "),
    ". Refusing to overwrite automatically. If replacement is intentional, ",
    "set INPH_OVERWRITE_MARKER_IFN_STATES=1 for that run."
  )
}

# -----------------------------------------------------------------------------
# Section 2 · Frozen revision-stage configuration
# -----------------------------------------------------------------------------

RANDOM_SEED <- 1L
MIN_CELLS_PER_STATE_PER_DONOR <- 30L
EXPECTED_CANONICAL_CELLS <- 22197L
EXPECTED_UNIFORM_CELLS <- 24703L
EXPECTED_CANONICAL_DONORS <- 9L
EXPECTED_UNIFORM_DONORS <- 10L
EXPECTED_CANONICAL_IFN_HIGH <- 3195L
EXPECTED_PRIMARY_CLUSTER5 <- 940L
EXPECTED_HIGH_AND_CLUSTER5 <- 898L

EXPECTED_ISG_SELECTOR <- c(
  "ISG15", "IFI6", "IFI44", "IFI44L", "MX1",
  "MX2", "OAS1", "OAS2", "OAS3", "IFIT1",
  "IFIT3", "RSAD2", "STAT1", "USP18", "LY6E"
)

# -----------------------------------------------------------------------------
# Section 3 · Load the frozen Phase 2 state definition
# -----------------------------------------------------------------------------

phase2 <- readRDS(PHASE2_PATH)

required_phase2_fields <- c(
  "engine_A_primary",
  "engine_B_sensitivity",
  "inph_cell_scores",
  "params",
  "excluded_genes",
  "input_provenance"
)
missing_phase2_fields <- setdiff(required_phase2_fields, names(phase2))

if (length(missing_phase2_fields) > 0L) {
  stop(
    "Phase 2 result is missing required field(s): ",
    paste(missing_phase2_fields, collapse = ", "),
    ". Run the revised 02_cross_disease_contextualization.R first."
  )
}

if (!identical(
  as.character(phase2$excluded_genes),
  EXPECTED_ISG_SELECTOR
)) {
  stop("The saved Phase 2 ISG selector does not match the frozen panel.")
}

AUC_MAXRANK <- as.integer(phase2$params$auc_maxrank)
IFN_THRESHOLD <- as.numeric(phase2$engine_A_primary$thr_shared)

if (length(AUC_MAXRANK) != 1L ||
    is.na(AUC_MAXRANK) ||
    AUC_MAXRANK != 1500L) {
  stop("Expected the saved Phase 2 aucMaxRank to equal 1500.")
}

if (length(IFN_THRESHOLD) != 1L || !is.finite(IFN_THRESHOLD)) {
  stop("The saved Phase 2 AUCell threshold is missing or non-finite.")
}

canonical_scores <- as.data.frame(
  phase2$inph_cell_scores,
  stringsAsFactors = FALSE
)
required_score_columns <- c(
  "cell_id", "donor", "is_c5", "A_score", "A_isg",
  "B_score", "B_isg"
)
missing_score_columns <- setdiff(
  required_score_columns,
  colnames(canonical_scores)
)

if (length(missing_score_columns) > 0L) {
  stop(
    "Saved canonical cell scores are missing column(s): ",
    paste(missing_score_columns, collapse = ", ")
  )
}

if (nrow(canonical_scores) != EXPECTED_CANONICAL_CELLS) {
  stop(
    "Expected ", EXPECTED_CANONICAL_CELLS,
    " canonical scores; observed ", nrow(canonical_scores), "."
  )
}

if (anyDuplicated(canonical_scores$cell_id)) {
  stop("Canonical cell-score identifiers are not unique.")
}

if (length(unique(as.character(canonical_scores$donor))) !=
    EXPECTED_CANONICAL_DONORS) {
  stop("Canonical cell scores should represent nine donors.")
}

canonical_assignments <- dplyr::mutate(
  canonical_scores,
  patient = as.character(donor),
  marker_ifn_score = as.numeric(A_isg),
  marker_ifn_state = ifelse(
    marker_ifn_score >= IFN_THRESHOLD,
    "marker_ifn_high",
    "marker_ifn_low"
  ),
  score_source = "exact_phase2_canonical_score"
)

canonical_high <- canonical_assignments$marker_ifn_state ==
  "marker_ifn_high"
canonical_c5 <- as.logical(canonical_assignments$is_c5)

canonical_overlap_counts <- c(
  sum(canonical_high),
  sum(canonical_c5),
  sum(canonical_high & canonical_c5)
)

if (!identical(
  as.integer(canonical_overlap_counts),
  c(
    EXPECTED_CANONICAL_IFN_HIGH,
    EXPECTED_PRIMARY_CLUSTER5,
    EXPECTED_HIGH_AND_CLUSTER5
  )
)) {
  stop(
    "Canonical marker-state assignments do not reproduce the saved ",
    "Phase 2 overlap counts."
  )
}

# -----------------------------------------------------------------------------
# Section 4 · Score the uniformly reconstructed ten-patient population
# -----------------------------------------------------------------------------

patient6_result <- readRDS(PATIENT6_RESULT_PATH)

if (!isTRUE(
  patient6_result$overall_validation$eligible_for_downstream
)) {
  stop(
    "Analysis 01C did not mark the uniform population as eligible ",
    "for downstream analysis."
  )
}

uniform_monocytes <- readRDS(UNIFORM_MONOCYTES_PATH)

if (!inherits(uniform_monocytes, "SingleCellExperiment")) {
  stop(
    "Uniform monocyte input must be a SingleCellExperiment; observed: ",
    paste(class(uniform_monocytes), collapse = ", ")
  )
}

if (ncol(uniform_monocytes) != EXPECTED_UNIFORM_CELLS) {
  stop(
    "Expected ", EXPECTED_UNIFORM_CELLS,
    " uniform candidate monocytes; observed ",
    ncol(uniform_monocytes), "."
  )
}

uniform_patient <- as.character(uniform_monocytes$patient)

if (length(unique(uniform_patient)) != EXPECTED_UNIFORM_DONORS) {
  stop("Uniform candidate-monocyte input should represent ten patients.")
}

if (!"6" %in% uniform_patient) {
  stop("Patient 6 is absent from the uniform candidate-monocyte input.")
}

if (!"counts" %in%
    SummarizedExperiment::assayNames(uniform_monocytes)) {
  stop("Uniform candidate-monocyte input lacks a counts assay.")
}

uniform_counts <- SummarizedExperiment::assay(
  uniform_monocytes,
  "counts"
)

missing_selector_genes <- setdiff(
  EXPECTED_ISG_SELECTOR,
  rownames(uniform_counts)
)

if (length(missing_selector_genes) > 0L) {
  stop(
    "Uniform counts are missing selector gene(s): ",
    paste(missing_selector_genes, collapse = ", ")
  )
}

cat(
  "Scoring ", ncol(uniform_counts),
  " uniform candidate monocytes with AUCell ...\n",
  sep = ""
)

set.seed(RANDOM_SEED)

uniform_rankings <- AUCell::AUCell_buildRankings(
  uniform_counts,
  plotStats = FALSE,
  verbose = FALSE
)
uniform_auc <- AUCell::AUCell_calcAUC(
  list(isg_selector = EXPECTED_ISG_SELECTOR),
  uniform_rankings,
  aucMaxRank = AUC_MAXRANK,
  verbose = FALSE
)
uniform_auc_matrix <- AUCell::getAUC(uniform_auc)

if (!"isg_selector" %in% rownames(uniform_auc_matrix)) {
  stop("AUCell output does not contain the ISG selector score.")
}

uniform_cell_ids <- colnames(uniform_counts)
uniform_scores <- as.numeric(
  uniform_auc_matrix["isg_selector", uniform_cell_ids]
)

if (length(uniform_scores) != ncol(uniform_monocytes) ||
    any(!is.finite(uniform_scores))) {
  stop("Uniform AUCell scoring returned missing or non-finite values.")
}

uniform_moca <- if ("moca_group" %in% colnames(
  SummarizedExperiment::colData(uniform_monocytes)
)) {
  as.character(uniform_monocytes$moca_group)
} else {
  rep(NA_character_, ncol(uniform_monocytes))
}

uniform_assignments <- data.frame(
  cell_id = uniform_cell_ids,
  patient = uniform_patient,
  moca_group = uniform_moca,
  marker_ifn_score = uniform_scores,
  marker_ifn_state = ifelse(
    uniform_scores >= IFN_THRESHOLD,
    "marker_ifn_high",
    "marker_ifn_low"
  ),
  score_source = paste0(
    "uniform_ten_patient_AUCell_seed_",
    RANDOM_SEED
  ),
  is_canonical_cell =
    uniform_cell_ids %in% canonical_assignments$cell_id,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Section 5 · Donor summaries and concordance diagnostics
# -----------------------------------------------------------------------------

summarize_donors <- function(assignments, cohort_label) {
  patient_ids <- unique(as.character(assignments$patient))
  patient_ids <- patient_ids[order(as.numeric(patient_ids))]

  dplyr::bind_rows(lapply(patient_ids, function(patient_id) {
    donor_data <- assignments[
      as.character(assignments$patient) == patient_id,
      ,
      drop = FALSE
    ]

    n_high <- sum(
      donor_data$marker_ifn_state == "marker_ifn_high"
    )
    n_low <- sum(
      donor_data$marker_ifn_state == "marker_ifn_low"
    )

    data.frame(
      cohort = cohort_label,
      patient = patient_id,
      total_cells = nrow(donor_data),
      marker_ifn_high_n = n_high,
      marker_ifn_low_n = n_low,
      marker_ifn_high_pct = 100 * n_high / nrow(donor_data),
      mean_marker_ifn_score =
        mean(donor_data$marker_ifn_score),
      median_marker_ifn_score =
        stats::median(donor_data$marker_ifn_score),
      high_meets_30_cell_minimum =
        n_high >= MIN_CELLS_PER_STATE_PER_DONOR,
      low_meets_30_cell_minimum =
        n_low >= MIN_CELLS_PER_STATE_PER_DONOR,
      both_states_meet_30_cell_minimum =
        n_high >= MIN_CELLS_PER_STATE_PER_DONOR &&
        n_low >= MIN_CELLS_PER_STATE_PER_DONOR,
      stringsAsFactors = FALSE
    )
  }))
}

canonical_donor_summary <- summarize_donors(
  canonical_assignments,
  "canonical_nine_patient"
)
uniform_donor_summary <- summarize_donors(
  uniform_assignments,
  "uniform_ten_patient"
)
donor_summary <- dplyr::bind_rows(
  canonical_donor_summary,
  uniform_donor_summary
)

canonical_cluster_concordance <- data.frame(
  cohort = "canonical_nine_patient",
  marker_ifn_high_n = sum(canonical_high),
  cluster5_n = sum(canonical_c5),
  marker_high_and_cluster5_n =
    sum(canonical_high & canonical_c5),
  marker_high_outside_cluster5_n =
    sum(canonical_high & !canonical_c5),
  cluster5_below_marker_threshold_n =
    sum(!canonical_high & canonical_c5),
  pct_marker_high_in_cluster5 =
    100 * sum(canonical_high & canonical_c5) /
      sum(canonical_high),
  pct_cluster5_marker_high =
    100 * sum(canonical_high & canonical_c5) /
      sum(canonical_c5),
  jaccard =
    sum(canonical_high & canonical_c5) /
      sum(canonical_high | canonical_c5),
  stringsAsFactors = FALSE
)

common_ids <- intersect(
  canonical_assignments$cell_id,
  uniform_assignments$cell_id
)
canonical_common <- canonical_assignments[
  match(common_ids, canonical_assignments$cell_id),
  ,
  drop = FALSE
]
uniform_common <- uniform_assignments[
  match(common_ids, uniform_assignments$cell_id),
  ,
  drop = FALSE
]

if (!identical(canonical_common$cell_id, uniform_common$cell_id)) {
  stop("Failed to align shared canonical and uniform cell identifiers.")
}

state_agreement <- canonical_common$marker_ifn_state ==
  uniform_common$marker_ifn_state

scoring_concordance <- data.frame(
  shared_cells_n = length(common_ids),
  exact_state_agreement_n = sum(state_agreement),
  exact_state_agreement_pct = 100 * mean(state_agreement),
  canonical_high_n = sum(
    canonical_common$marker_ifn_state == "marker_ifn_high"
  ),
  uniform_high_n = sum(
    uniform_common$marker_ifn_state == "marker_ifn_high"
  ),
  high_in_both_n = sum(
    canonical_common$marker_ifn_state == "marker_ifn_high" &
      uniform_common$marker_ifn_state == "marker_ifn_high"
  ),
  score_pearson = stats::cor(
    canonical_common$marker_ifn_score,
    uniform_common$marker_ifn_score,
    method = "pearson"
  ),
  score_spearman = stats::cor(
    canonical_common$marker_ifn_score,
    uniform_common$marker_ifn_score,
    method = "spearman"
  ),
  mean_absolute_score_difference = mean(abs(
    canonical_common$marker_ifn_score -
      uniform_common$marker_ifn_score
  )),
  maximum_absolute_score_difference = max(abs(
    canonical_common$marker_ifn_score -
      uniform_common$marker_ifn_score
  )),
  stringsAsFactors = FALSE
)

patient6_mask <- uniform_assignments$patient == "6"

population_summary <- dplyr::bind_rows(
  data.frame(
    cohort = "canonical_nine_patient",
    patients_n = length(unique(canonical_assignments$patient)),
    cells_n = nrow(canonical_assignments),
    marker_ifn_high_n = sum(canonical_high),
    marker_ifn_low_n = sum(!canonical_high),
    marker_ifn_high_pct = 100 * mean(canonical_high),
    all_donors_both_states_ge_30 = all(
      canonical_donor_summary$both_states_meet_30_cell_minimum
    ),
    patient6_cells_n = 0L,
    patient6_marker_ifn_high_n = 0L,
    patient6_marker_ifn_high_pct = NA_real_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "uniform_ten_patient",
    patients_n = length(unique(uniform_assignments$patient)),
    cells_n = nrow(uniform_assignments),
    marker_ifn_high_n = sum(
      uniform_assignments$marker_ifn_state ==
        "marker_ifn_high"
    ),
    marker_ifn_low_n = sum(
      uniform_assignments$marker_ifn_state ==
        "marker_ifn_low"
    ),
    marker_ifn_high_pct = 100 * mean(
      uniform_assignments$marker_ifn_state ==
        "marker_ifn_high"
    ),
    all_donors_both_states_ge_30 = all(
      uniform_donor_summary$both_states_meet_30_cell_minimum
    ),
    patient6_cells_n = sum(patient6_mask),
    patient6_marker_ifn_high_n = sum(
      patient6_mask &
        uniform_assignments$marker_ifn_state ==
          "marker_ifn_high"
    ),
    patient6_marker_ifn_high_pct = 100 * mean(
      uniform_assignments$marker_ifn_state[patient6_mask] ==
        "marker_ifn_high"
    ),
    stringsAsFactors = FALSE
  )
)

# -----------------------------------------------------------------------------
# Section 6 · Save the locked state definitions
# -----------------------------------------------------------------------------

result <- list(
  analysis_id = "revision_02B_marker_ifn_states",
  configuration = list(
    primary_scoring_engine = "AUCell on raw counts",
    isg_selector = EXPECTED_ISG_SELECTOR,
    auc_maxrank = AUC_MAXRANK,
    random_seed_uniform_arm = RANDOM_SEED,
    threshold = IFN_THRESHOLD,
    threshold_basis = phase2$params$threshold_basis,
    classification_rule = paste(
      "marker_ifn_high if AUCell ISG-selector score >= threshold;",
      "otherwise marker_ifn_low"
    ),
    minimum_cells_per_state_per_donor =
      MIN_CELLS_PER_STATE_PER_DONOR,
    canonical_score_source = paste(
      "Exact cell-level scores retained from",
      "02_cross_disease_contextualization.R"
    ),
    uniform_score_source = paste(
      "Recomputed once across the complete uniformly reconstructed",
      "ten-patient candidate-monocyte population"
    )
  ),
  canonical_assignments = canonical_assignments,
  uniform_ten_patient_assignments = uniform_assignments,
  population_summary = population_summary,
  donor_summary = donor_summary,
  canonical_cluster_concordance =
    canonical_cluster_concordance,
  canonical_uniform_scoring_concordance =
    scoring_concordance,
  downstream_eligibility = list(
    canonical_marker_de_eligible = all(
      canonical_donor_summary$both_states_meet_30_cell_minimum
    ),
    uniform_ten_patient_marker_de_eligible = all(
      uniform_donor_summary$both_states_meet_30_cell_minimum
    ),
    patient6_present =
      "6" %in% uniform_assignments$patient,
    patient6_both_states_meet_30_cell_minimum =
      uniform_donor_summary$both_states_meet_30_cell_minimum[
        uniform_donor_summary$patient == "6"
      ]
  ),
  input_provenance = list(
    phase2_path =
      "data/processed/phase2_contextualization.rds",
    phase2_md5 = unname(tools::md5sum(PHASE2_PATH)),
    uniform_monocytes_path = paste0(
      "data/processed/",
      "revision_patient6_uniform_monocytes.rds"
    ),
    uniform_monocytes_md5 = unname(
      tools::md5sum(UNIFORM_MONOCYTES_PATH)
    ),
    patient6_result_path = paste0(
      "data/processed/",
      "revision_patient6_sensitivity.rds"
    ),
    patient6_result_md5 = unname(
      tools::md5sum(PATIENT6_RESULT_PATH)
    )
  )
)

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(
  REVISION_RESULTS_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(result, OUTPUT_RDS)
utils::write.csv(
  population_summary,
  OUTPUT_POPULATION_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  donor_summary,
  OUTPUT_DONOR_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  canonical_cluster_concordance,
  OUTPUT_CLUSTER_CONCORDANCE_CSV,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  scoring_concordance,
  OUTPUT_SCORING_CONCORDANCE_CSV,
  row.names = FALSE,
  na = ""
)

if (!all(file.exists(OUTPUT_PATHS))) {
  stop("One or more marker-state outputs were not created.")
}

if (any(is.na(file.info(OUTPUT_PATHS)$size)) ||
    any(file.info(OUTPUT_PATHS)$size <= 0)) {
  stop("One or more marker-state outputs are empty.")
}

result_check <- readRDS(OUTPUT_RDS)

if (!identical(
  result_check$canonical_assignments$cell_id,
  canonical_assignments$cell_id
)) {
  stop("Saved canonical marker-state order failed validation.")
}

if (!identical(
  result_check$uniform_ten_patient_assignments$cell_id,
  uniform_assignments$cell_id
)) {
  stop("Saved uniform marker-state order failed validation.")
}

# -----------------------------------------------------------------------------
# Section 7 · Closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  REVISION ANALYSIS 2B COMPLETE\n")
cat("============================================================\n")
cat(
  "  Frozen AUCell threshold           : ",
  format(IFN_THRESHOLD, digits = 7), "\n",
  sep = ""
)
cat(
  "  ISG selector genes                : ",
  length(EXPECTED_ISG_SELECTOR), "\n",
  sep = ""
)
cat(
  "  Canonical marker IFN-high cells   : ",
  sum(canonical_high), " / ", nrow(canonical_assignments), "\n",
  sep = ""
)
cat(
  "  Cluster 5 marker IFN-high         : ",
  sum(canonical_high & canonical_c5),
  " / ", sum(canonical_c5), "\n",
  sep = ""
)
cat(
  "  Uniform ten-patient IFN-high cells: ",
  sum(
    uniform_assignments$marker_ifn_state ==
      "marker_ifn_high"
  ),
  " / ", nrow(uniform_assignments), "\n",
  sep = ""
)
cat(
  "  Patient 6 IFN-high cells          : ",
  population_summary$patient6_marker_ifn_high_n[
    population_summary$cohort == "uniform_ten_patient"
  ],
  " / ",
  population_summary$patient6_cells_n[
    population_summary$cohort == "uniform_ten_patient"
  ],
  "\n",
  sep = ""
)
cat(
  "  Shared-cell state agreement       : ",
  round(scoring_concordance$exact_state_agreement_pct, 2),
  "%\n",
  sep = ""
)
cat(
  "  Canonical marker-DE eligible      : ",
  result$downstream_eligibility$canonical_marker_de_eligible,
  "\n",
  sep = ""
)
cat(
  "  Ten-patient marker-DE eligible    : ",
  result$downstream_eligibility[[
    "uniform_ten_patient_marker_de_eligible"
  ]],
  "\n",
  sep = ""
)
cat("\nDonor-level state counts:\n")
print(donor_summary, row.names = FALSE)
cat("\nCanonical cluster/marker concordance:\n")
print(canonical_cluster_concordance, row.names = FALSE)
cat("\nCanonical-versus-uniform scoring concordance:\n")
print(scoring_concordance, row.names = FALSE)
cat("\nOutputs saved:\n")
cat("  ", OUTPUT_RDS, "\n", sep = "")
cat("  ", OUTPUT_POPULATION_CSV, "\n", sep = "")
cat("  ", OUTPUT_DONOR_CSV, "\n", sep = "")
cat("  ", OUTPUT_CLUSTER_CONCORDANCE_CSV, "\n", sep = "")
cat("  ", OUTPUT_SCORING_CONCORDANCE_CSV, "\n", sep = "")
cat("============================================================\n")
