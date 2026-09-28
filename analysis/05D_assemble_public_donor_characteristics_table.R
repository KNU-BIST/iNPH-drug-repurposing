# =============================================================================
# Revision analysis 05D — Public donor-characteristics availability table
# Install as analysis/05D_assemble_public_donor_characteristics_table.R
# Build: 05D_PUBLIC_DONOR_CHARACTERISTICS_V1_20260907
# =============================================================================
#
# Purpose
# -------
# Assemble the supplementary donor table requested during peer review using
# only characteristics that can be linked reproducibly to the public
# GSE292141 sample identifiers. Missing clinical descriptors are reported as
# not publicly available; they are never inferred or imputed.
#
# This script does not rerun single-cell analysis, test MoCA associations, or
# change any existing result. It combines saved revision outputs and records
# the limits of the public donor metadata.
# =============================================================================

local({
  BUILD <- "05D_PUBLIC_DONOR_CHARACTERISTICS_V1_20260907"

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  required_inputs <- c(
    patient6 = "data/processed/revision_patient6_sensitivity.rds",
    marker_states = "data/processed/revision_marker_ifn_states.rds",
    phase2 = "data/processed/phase2_contextualization.rds",
    csf = "data/processed/revision_csf_ifn_contextualization.rds"
  )

  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  check(
    length(missing_inputs) == 0L,
    paste(
      "Missing required saved result(s):",
      paste(unname(missing_inputs), collapse = ", ")
    )
  )

  input_md5_before <- unname(tools::md5sum(required_inputs))
  check(!anyNA(input_md5_before), "An input checksum could not be calculated.")

  patient6 <- readRDS(required_inputs[["patient6"]])
  marker <- readRDS(required_inputs[["marker_states"]])
  phase2 <- readRDS(required_inputs[["phase2"]])
  csf <- readRDS(required_inputs[["csf"]])

  check(
    identical(patient6$analysis_id, "revision_01C_patient6_sensitivity"),
    "Unexpected Patient 6 result."
  )
  check(
    identical(marker$analysis_id, "revision_02B_marker_ifn_states"),
    "Unexpected marker-state result."
  )
  check(
    identical(csf$analysis_id, "revision_05B_csf_monocyte_ifn_contextualization"),
    "Unexpected CSF result."
  )
  check(
    "inph_cell_scores" %in% names(phase2),
    "Phase 2 result does not contain the preserved cell-level iNPH table."
  )

  # Public sample identifiers and group assignments are the verified metadata
  # encoded in R/data_acquisition_and_loading.R and GEO GSE292141.
  sample_map <- data.frame(
    patient = as.character(1:10),
    csf_GSM = c(
      "GSM8851048", "GSM8851050", "GSM8851052", "GSM8851054",
      "GSM8851056", "GSM8851058", "GSM8851060", "GSM8851062",
      "GSM8851064", "GSM8851066"
    ),
    pbmc_GSM = c(
      "GSM8851049", "GSM8851051", "GSM8851053", "GSM8851055",
      "GSM8851057", "GSM8851059", "GSM8851061", "GSM8851063",
      "GSM8851065", "GSM8851067"
    ),
    moca_group = c(
      rep("High", 5), "Unknown", rep("Low", 4)
    ),
    stringsAsFactors = FALSE
  )

  # The paper reports the exact MoCA values only as two unordered group-level
  # sets. It does not publicly link an exact value to each patient/GSM.
  sample_map$exact_moca_score <- NA_integer_
  sample_map$public_group_moca_values <- ifelse(
    sample_map$moca_group == "High",
    "23; 23; 25; 27; 29 (group-level values; patient mapping unavailable)",
    ifelse(
      sample_map$moca_group == "Low",
      "19; 20; 20; 20 (group-level values; patient mapping unavailable)",
      "Not reported"
    )
  )
  sample_map$age_years <- NA_real_
  sample_map$sex <- NA_character_
  sample_map$relevant_medications <- NA_character_

  p6_summary <- as.data.frame(patient6$patient_summary)
  required_p6_columns <- c(
    "patient", "moca_group", "raw_pbmc_n", "post_qc_pbmc_n",
    "candidate_monocyte_n"
  )
  check(
    all(required_p6_columns %in% names(p6_summary)),
    "Unexpected Patient 6 donor-summary schema."
  )
  p6_summary$patient <- as.character(p6_summary$patient)
  check(
    setequal(p6_summary$patient, sample_map$patient),
    "Patient 6 result does not contain patients 1 through 10 exactly once."
  )
  p6_summary <- p6_summary[
    match(sample_map$patient, p6_summary$patient),
    required_p6_columns,
    drop = FALSE
  ]
  check(
    identical(as.character(p6_summary$moca_group), sample_map$moca_group),
    "Public MoCA-group mapping disagrees with the saved ten-patient result."
  )

  marker_donors <- as.data.frame(marker$donor_summary)
  required_marker_columns <- c(
    "cohort", "patient", "total_cells", "marker_ifn_high_n",
    "marker_ifn_low_n", "marker_ifn_high_pct"
  )
  check(
    all(required_marker_columns %in% names(marker_donors)),
    "Unexpected marker-state donor-summary schema."
  )
  marker_donors$patient <- as.character(marker_donors$patient)

  canonical <- marker_donors[
    marker_donors$cohort == "canonical_nine_patient",
    required_marker_columns,
    drop = FALSE
  ]
  uniform <- marker_donors[
    marker_donors$cohort == "uniform_ten_patient",
    required_marker_columns,
    drop = FALSE
  ]
  check(
    setequal(canonical$patient, setdiff(sample_map$patient, "6")) &&
      setequal(uniform$patient, sample_map$patient),
    "Marker-state donor universes are inconsistent with the locked cohorts."
  )

  inph_scores <- as.data.frame(phase2$inph_cell_scores)
  check(
    "is_c5" %in% names(inph_scores) &&
      any(c("patient", "donor") %in% names(inph_scores)),
    "Unexpected Phase 2 iNPH cell-score schema."
  )
  if (!"patient" %in% names(inph_scores)) {
    inph_scores$patient <- inph_scores$donor
  }
  inph_scores$patient <- as.character(inph_scores$patient)
  c5_counts <- aggregate(
    as.integer(inph_scores$is_c5),
    by = list(patient = inph_scores$patient),
    FUN = sum
  )
  names(c5_counts)[2] <- "canonical_cluster5_n"

  csf_donors <- as.data.frame(csf$csf_donor_summary)
  required_csf_columns <- c(
    "patient", "n_cells", "marker_ifn_high_n", "marker_ifn_high_pct",
    "median_transfer_score", "median_marker_ifn_score"
  )
  check(
    all(required_csf_columns %in% names(csf_donors)),
    "Unexpected CSF donor-summary schema."
  )
  csf_donors$patient <- as.character(csf_donors$patient)
  check(
    setequal(csf_donors$patient, sample_map$patient),
    "CSF result does not contain patients 1 through 10 exactly once."
  )

  pick_by_patient <- function(data, column) {
    data[[column]][match(sample_map$patient, data$patient)]
  }

  donor_table <- sample_map
  donor_table$age_public_availability <- "Not reported at donor level"
  donor_table$sex_public_availability <- "Not reported at donor level"
  donor_table$medication_public_availability <- "Not reported at donor level"
  donor_table$exact_moca_public_availability <- ifelse(
    donor_table$moca_group == "Unknown",
    "Neither exact score nor group publicly reported",
    "Exact values reported only as an unordered group-level set"
  )
  donor_table$raw_PBMC_cells_loaded <- p6_summary$raw_pbmc_n
  donor_table$post_QC_PBMC_cells <- p6_summary$post_qc_pbmc_n
  donor_table$uniform_PBMC_candidate_monocytes <-
    p6_summary$candidate_monocyte_n
  donor_table$canonical_PBMC_monocytes <- pick_by_patient(
    canonical, "total_cells"
  )
  donor_table$canonical_marker_IFN_high_cells <- pick_by_patient(
    canonical, "marker_ifn_high_n"
  )
  donor_table$canonical_cluster5_cells <- pick_by_patient(
    c5_counts, "canonical_cluster5_n"
  )
  donor_table$CSF_candidate_monocytes <- pick_by_patient(
    csf_donors, "n_cells"
  )
  donor_table$CSF_marker_IFN_high_cells <- pick_by_patient(
    csf_donors, "marker_ifn_high_n"
  )
  donor_table$CSF_marker_IFN_high_pct <- pick_by_patient(
    csf_donors, "marker_ifn_high_pct"
  )
  donor_table$CSF_median_transfer_score <- pick_by_patient(
    csf_donors, "median_transfer_score"
  )
  donor_table$CSF_candidates_ge30 <-
    donor_table$CSF_candidate_monocytes >= 30L
  donor_table$primary_nine_patient_analysis <-
    donor_table$patient != "6"

  check(
    sum(donor_table$canonical_PBMC_monocytes, na.rm = TRUE) == 22197L &&
      sum(donor_table$canonical_cluster5_cells, na.rm = TRUE) == 940L &&
      sum(donor_table$canonical_marker_IFN_high_cells, na.rm = TRUE) == 3195L &&
      sum(donor_table$uniform_PBMC_candidate_monocytes) == 24703L &&
      sum(donor_table$CSF_candidate_monocytes) == 2683L &&
      sum(donor_table$CSF_marker_IFN_high_cells) == 408L,
    "One or more locked cell-count totals failed reconciliation."
  )

  availability <- data.frame(
    characteristic = c(
      "iNPH diagnosis", "paired tissue", "patient-linked MoCA group",
      "exact MoCA score", "age", "sex", "relevant medications",
      "PBMC cell numbers", "ventricular-CSF cell numbers"
    ),
    public_availability = c(
      "Available for all 10 donors",
      "Available for all 10 donors (PBMC and ventricular CSF)",
      "Available for 9 donors; Patient 6 is Unknown",
      paste(
        "Nine values are reported as unordered group-level sets;",
        "not linkable to individual patient/GSM IDs"
      ),
      "Not reported at donor level",
      "Not reported at donor level",
      "Not reported at donor level",
      "Available from reproducible reprocessing outputs",
      "Available from reproducible reprocessing outputs"
    ),
    reporting_action = c(
      "Report", "Report", "Report", "Do not assign or impute",
      "Mark NR", "Mark NR", "Mark NR", "Report exact analyzed counts",
      "Report exact analyzed counts with heterogeneity warning"
    ),
    stringsAsFactors = FALSE
  )

  manuscript_note <- paste(
    "Supplementary Table Sx summarizes all donor descriptors that could be",
    "linked reproducibly to the public GSE292141 sample identifiers. MoCA",
    "group was available for nine donors, whereas Patient 6 had unknown MoCA",
    "status. The source article reported the nine exact MoCA values only as",
    "unordered High- and Low-MoCA group-level sets, precluding assignment to",
    "individual patient identifiers. Donor-level age, sex, and medication",
    "information was not publicly reported and is marked NR rather than",
    "inferred. Exact analyzed PBMC and ventricular-CSF cell counts are shown."
  )

  reviewer_response <- paste(
    "We have added a supplementary donor-characteristics table containing",
    "the publicly linkable MoCA-group, sample-accession, tissue, and analyzed",
    "cell-count information for all ten donors. We agree that age, sex, exact",
    "MoCA score, and relevant medication data would aid interpretation.",
    "However, donor-level age, sex, and medication information was not",
    "available in the public source, and exact MoCA values were reported only",
    "as unordered values within the High- and Low-MoCA groups rather than",
    "linked to patient or GSM identifiers. We therefore mark these fields as",
    "not reported and do not infer or impute them."
  )

  parent <- "results/revision/public_donor_characteristics"
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  output_directory <- tempfile(
    paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = parent
  )
  check(dir.create(output_directory), "Cannot create the 05D output directory.")

  result <- list(
    analysis_id = "revision_05D_public_donor_characteristics",
    build_id = BUILD,
    status = "PASS_PUBLICLY_AVAILABLE_FIELDS_ONLY_NO_IMPUTATION",
    donor_characteristics = donor_table,
    availability_audit = availability,
    manuscript_note = manuscript_note,
    response_to_reviewer = reviewer_response,
    public_sources = data.frame(
      source = c("GSE292141", "Duy et al., PNAS 2025"),
      url = c(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292141",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC12087963/"
      ),
      use = c(
        "Patient, tissue, GSM, and MoCA-group linkage",
        "Cohort description and unordered group-level MoCA values"
      ),
      stringsAsFactors = FALSE
    ),
    input_provenance = data.frame(
      input = unname(required_inputs),
      md5 = input_md5_before,
      stringsAsFactors = FALSE
    ),
    output_directory = output_directory,
    session_info = utils::sessionInfo()
  )

  donor_csv <- file.path(
    output_directory,
    "supplementary_donor_characteristics_public_data.csv"
  )
  availability_csv <- file.path(
    output_directory,
    "donor_characteristic_availability_audit.csv"
  )
  wording_txt <- file.path(
    output_directory,
    "donor_characteristics_reportable_text.txt"
  )
  result_rds <- file.path(
    output_directory,
    "revision_public_donor_characteristics.rds"
  )

  utils::write.csv(donor_table, donor_csv, row.names = FALSE, na = "NR")
  utils::write.csv(availability, availability_csv, row.names = FALSE, na = "NR")
  writeLines(
    c(
      "MANUSCRIPT NOTE", manuscript_note, "",
      "RESPONSE TO REVIEWER", reviewer_response
    ),
    wording_txt
  )
  saveRDS(result, result_rds)
  check(identical(readRDS(result_rds), result), "05D RDS read-back failed.")

  check(
    identical(input_md5_before, unname(tools::md5sum(required_inputs))),
    "An input changed during Analysis 05D."
  )
  exports <- list.files(output_directory, full.names = TRUE)
  check(all(file.info(exports)$size > 0), "An Analysis 05D export is empty.")
  utils::write.csv(
    data.frame(
      file = basename(exports),
      md5 = unname(tools::md5sum(exports)),
      stringsAsFactors = FALSE
    ),
    file.path(output_directory, "output_checksums.csv"),
    row.names = FALSE
  )
  writeLines(
    c(result$status, BUILD),
    file.path(output_directory, "_SUCCESS.txt")
  )

  cat("\n============================================================\n")
  cat("  REVISION ANALYSIS 05D COMPLETE\n")
  cat("============================================================\n")
  cat("  Status                         :", result$status, "\n")
  cat("  Donors tabulated               :", nrow(donor_table), "\n")
  cat("  MoCA group available           : 9 / 10\n")
  cat("  Exact MoCA linked by patient   : 0 / 10\n")
  cat("  Age linked by patient          : 0 / 10\n")
  cat("  Sex linked by patient          : 0 / 10\n")
  cat("  Medications linked by patient  : 0 / 10\n")
  cat("  Cell-count totals reconciled   : PASS\n")
  cat("  Clinical fields inferred       : 0\n")
  cat("\nOutput directory:", output_directory, "\n")
  cat("============================================================\n")
})
