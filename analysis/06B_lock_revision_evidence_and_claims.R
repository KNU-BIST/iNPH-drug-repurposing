# =============================================================================
# Revision analysis 06B — Lock revision evidence and claims
# Install as analysis/06B_lock_revision_evidence_and_claims.R
# Build: 06B_REVISION_EVIDENCE_CLAIM_LOCK_V1_20260907
# =============================================================================
#
# This is a read-only integration step. It runs no biological scoring, no
# differential expression, and no LINCS query. It validates completed revision
# outputs, freezes reportable values and claim boundaries, and maps the evidence
# to reviewer-response actions.
# =============================================================================

local({
  BUILD <- "06B_REVISION_EVIDENCE_CLAIM_LOCK_V1_20260907"
  MANIFEST_SHA <-
    "1e3a2b827a1429e406303871d547a9dec07ee4068078c666df69ca82e5428380"

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }
  scalar <- function(x, message) {
    check(length(x) == 1L && !is.na(x), message)
    x[[1L]]
  }
  fmt <- function(x, digits = 3L) {
    if (is.numeric(x)) format(round(x, digits), trim = TRUE, scientific = FALSE) else as.character(x)
  }
  latest_successful <- function(parent, filename, analysis_id, status = NULL) {
    paths <- list.files(
      parent, pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
      recursive = TRUE, full.names = TRUE
    )
    paths <- paths[file.exists(file.path(dirname(paths), "_SUCCESS.txt"))]
    check(length(paths) > 0L, paste("No successful", filename, "found under", parent))
    paths <- paths[order(file.info(paths)$mtime, decreasing = TRUE)]
    for (path in paths) {
      object <- tryCatch(readRDS(path), error = function(e) NULL)
      if (is.list(object) && identical(object$analysis_id, analysis_id) &&
          (is.null(status) || identical(object$status, status))) {
        return(list(path = path, object = object))
      }
    }
    stop("No successful output matched the expected analysis ID/status: ", analysis_id,
         call. = FALSE)
  }
  row_value <- function(data, column, filter_column = NULL, filter_value = NULL,
                        message = "Expected one value") {
    check(is.data.frame(data) && column %in% names(data), paste(message, "(missing column)"))
    if (!is.null(filter_column)) {
      check(filter_column %in% names(data), paste(message, "(missing filter column)"))
      data <- data[as.character(data[[filter_column]]) == as.character(filter_value), , drop = FALSE]
    }
    scalar(data[[column]], message)
  }
  add_result <- function(id, section, metric, value, unit, scope, source_id,
                         source_path, claim_status, reportable_wording) {
    data.frame(
      result_id = id, section = section, metric = metric,
      value = as.character(value), unit = unit, scope = scope,
      source_analysis = source_id, source_path = source_path,
      claim_status = claim_status, reportable_wording = reportable_wording,
      stringsAsFactors = FALSE
    )
  }

  fixed_paths <- c(
    resolution = "data/processed/revision_state_resolution_sensitivity.rds",
    patient6 = "data/processed/revision_patient6_sensitivity.rds",
    marker = "data/processed/revision_marker_ifn_states.rds",
    variants = "data/processed/revision_signature_variants.rds",
    null_definitions = "data/processed/revision_expression_matched_null_signatures.rds",
    manifest = "data/processed/revision_signature_manifest.rds",
    batch = "data/processed/revision_lincs_batch_index.rds",
    robustness_original = "data/processed/revision_candidate_robustness_and_null_summary.rds",
    csf = "data/processed/revision_csf_ifn_contextualization.rds"
  )
  check(all(file.exists(fixed_paths)), paste(
    "Missing required fixed output(s):",
    paste(names(fixed_paths)[!file.exists(fixed_paths)], collapse = ", ")
  ))

  run_outputs <- list(
    rank_audit = latest_successful(
      "results/revision/zero_score_rank_tie_audit",
      "revision_zero_score_rank_tie_audit.rds",
      "revision_04D_zero_score_rank_tie_audit_v1",
      "PASS_AUDIT_COMPLETE_NOT_A_RANKING_REPLACEMENT"
    ),
    average_ties = latest_successful(
      "results/revision/average_tie_finalization",
      "revision_candidate_robustness_average_ties.rds",
      "revision_candidate_robustness_average_ties_v1", "PASS"
    ),
    curated = latest_successful(
      "results/revision/curated_compound_annotations",
      "revision_curated_compound_annotations.rds",
      "revision_curated_compound_annotations_v1",
      "PASS_TOP20_RESOLVED_QUALIFIED_REPORTING"
    ),
    exposure = latest_successful(
      "results/revision/lincs_exposure_metadata_audit",
      "revision_lincs_exposure_metadata_audit.rds",
      "revision_lincs_exposure_metadata_audit_v1",
      "PASS_EXPOSURE_METADATA_UNAVAILABLE_NO_UNSUPPORTED_LINKAGE"
    ),
    csf_audit = latest_successful(
      "results/revision/csf_cell_count_heterogeneity_audit",
      "revision_csf_cell_count_heterogeneity_audit.rds",
      "revision_05C_csf_cell_count_heterogeneity_audit"
    ),
    donor_table = latest_successful(
      "results/revision/public_donor_characteristics",
      "revision_public_donor_characteristics.rds",
      "revision_05D_public_donor_characteristics",
      "PASS_PUBLICLY_AVAILABLE_FIELDS_ONLY_NO_IMPUTATION"
    ),
    csf_figure = latest_successful(
      "results/revision/csf_contextualization_figures",
      "revision_csf_contextualization_figures.rds",
      "revision_05E_visualize_csf_contextualization",
      "PASS_SUPPORTIVE_SUPPLEMENTARY_FIGURE_CREATED"
    )
  )

  objects <- lapply(fixed_paths, readRDS)
  expected_ids <- c(
    resolution = "revision_01B_clustering_resolution_sensitivity",
    patient6 = "revision_01C_patient6_sensitivity",
    marker = "revision_02B_marker_ifn_states",
    variants = "revision_02C_signature_variants",
    null_definitions = "revision_02D_expression_matched_null_signatures",
    manifest = "revision_02E_signature_manifest",
    batch = "revision_03B_lincs_batch_index"
  )
  for (nm in names(expected_ids)) {
    check(identical(objects[[nm]]$analysis_id, expected_ids[[nm]]),
          paste("Unexpected analysis ID in", fixed_paths[[nm]]))
  }
  check(identical(objects$robustness_original$analysis_id, "revision_candidate_robustness_v1"),
        "Unexpected original robustness-summary analysis ID.")
  check(identical(objects$csf$analysis_id, "revision_05B_csf_monocyte_ifn_contextualization"),
        "Unexpected CSF contextualization analysis ID.")
  check(identical(objects$manifest$manifest_status, "LOCKED"), "Signature manifest is not locked.")
  check(identical(objects$manifest$manifest_content_sha256, MANIFEST_SHA),
        "Signature-manifest SHA-256 does not match the frozen lock.")
  check(identical(objects$batch$manifest_content_sha256, MANIFEST_SHA),
        "LINCS batch does not reference the frozen manifest.")
  check(isTRUE(objects$batch$all_locked_queries_complete),
        "Not all locked LINCS queries are complete.")
  check(objects$batch$full_reference_complete_n == 19L && objects$batch$null_complete_n == 100L,
        "LINCS completion counts are not 19 full-reference plus 100 null queries.")
  check(grepl("SUPPORTIVE_CSF_DETECTABILITY_ONLY", run_outputs$csf_audit$object$status,
              fixed = TRUE), "Unexpected CSF audit decision.")

  input_paths <- c(unname(fixed_paths), vapply(run_outputs, `[[`, character(1), "path"))
  own_script <- "analysis/06B_lock_revision_evidence_and_claims.R"
  if (file.exists(own_script)) input_paths <- c(input_paths, own_script)
  input_paths <- unique(input_paths)
  input_md5_before <- unname(tools::md5sum(input_paths))
  check(!anyNA(input_md5_before), "Could not calculate every input checksum.")

  resolution <- objects$resolution
  patient6 <- objects$patient6
  marker <- objects$marker
  rank_audit <- run_outputs$rank_audit$object
  average_ties <- run_outputs$average_ties$object
  curated <- run_outputs$curated$object
  exposure <- run_outputs$exposure$object
  csf <- objects$csf
  csf_audit <- run_outputs$csf_audit$object
  donor_table <- run_outputs$donor_table$object

  canonical <- marker$population_summary[
    marker$population_summary$cohort == "canonical_nine_patient", , drop = FALSE]
  uniform <- marker$population_summary[
    marker$population_summary$cohort == "uniform_ten_patient", , drop = FALSE]
  primary_zero <- rank_audit$score_counts[
    rank_audit$score_counts$query_id == "PRIMARY_150" &
      rank_audit$score_counts$population == "THP1_trt_cp", , drop = FALSE]
  priority <- average_ties$priority_candidate_summary
  rosu <- priority[tolower(priority$pert_iname) == "rosuvastatin", , drop = FALSE]
  rofl <- priority[tolower(priority$pert_iname) == "roflumilast", , drop = FALSE]
  check(nrow(rosu) == 1L && nrow(rofl) == 1L, "Priority compounds are not uniquely represented.")
  check(
    rosu$objective_robustness_rank == 1L && rofl$objective_robustness_rank == 15L &&
      abs(rosu$primary_NCS - (-1.845)) < 5e-4 &&
      abs(rofl$primary_NCS - (-1.375)) < 5e-4,
    "Priority compound scores or average-tie ranks differ from the reviewed results."
  )
  check(
    abs(rosu$empirical_lower_tail_p - 1 / 101) < 1e-12 &&
      abs(rofl$empirical_lower_tail_p - 1 / 101) < 1e-12,
    "Priority matched-null probabilities are not the reviewed plus-one values."
  )
  check(
    sum(rank_audit$all_compound_comparison$rank_change_average_minus_min != 0) == 885L &&
      all(rank_audit$top_k_overlap$overlap_n == rank_audit$top_k_overlap$top_k),
    "Average-tie all-compound or top-k results differ from the reviewed audit."
  )
  check(
    row_value(primary_zero, "NCS_negative_n") == 288L &&
      row_value(primary_zero, "NCS_zero_n") == 779L &&
      row_value(primary_zero, "NCS_positive_n") == 309L &&
      row_value(primary_zero, "NCS_missing_n") == 0L,
    "PRIMARY_150 zero-block accounting differs from the reviewed audit."
  )
  check(
    csf_audit$summary$primary_donors_ge30_n == 6L &&
      identical(as.character(csf_audit$summary$largest_donor), "10") &&
      abs(csf_audit$summary$largest_donor_share_pct - 77.12) < 0.01 &&
      abs(csf_audit$summary$effective_donor_number - 1.64) < 0.01 &&
      !isTRUE(csf_audit$summary$pooled_cell_weighted_percentage_for_cohort_reporting) &&
      !isTRUE(csf_audit$summary$definitive_CSF_enrichment_claim_supported),
    "CSF heterogeneity or inference guardrails differ from the reviewed audit."
  )
  check(nrow(curated$curated_top20) == 20L,
        "Curated top-20 annotation result differs from the reviewed output.")
  check(exposure$score_and_rank_changes_n == 0L && exposure$new_lincs_queries_n == 0L,
        "Exposure audit unexpectedly changed scores/ranks or ran a LINCS query.")

  reportable <- do.call(rbind, list(
    add_result("R01", "Population definition", "Canonical candidate monocytes", 22197, "cells",
      "Nine-patient canonical PBMC population", marker$analysis_id, fixed_paths[["marker"]], "report",
      "The canonical PBMC analysis contained 22,197 classical-monocyte candidates from nine donors."),
    add_result("R02", "IFN state", "Canonical marker-defined IFN-high", row_value(canonical, "marker_ifn_high_n"), "cells",
      "Frozen AUCell threshold in canonical PBMC candidates", marker$analysis_id, fixed_paths[["marker"]], "report",
      paste0("The frozen marker rule classified ", row_value(canonical, "marker_ifn_high_n"),
             " of ", row_value(canonical, "cells_n"), " canonical PBMC candidates as IFN-response-high.")),
    add_result("R03", "IFN state", "Cluster 5 marker-defined IFN-high", 898, "of 940 cells",
      "Concordance check; cluster 5 is not the primary state definition", marker$analysis_id,
      fixed_paths[["marker"]], "report_qualified",
      "The marker rule classified 898 of 940 cluster-5 cells as IFN-response-high (95.53%)."),
    add_result("R04", "Patient 6", "Uniform candidate population", patient6$overall_validation$candidate_n_9_patients +
      row_value(patient6$patient_summary, "candidate_monocyte_n", "patient", "6"), "cells",
      "All ten patients; sensitivity population", patient6$analysis_id, fixed_paths[["patient6"]], "report",
      "Uniform reconstruction yielded 24,703 candidate monocytes across all ten patients, including 2,004 from Patient 6."),
    add_result("R05", "Patient 6", "Canonical overlap Jaccard", fmt(patient6$overall_validation$jaccard, 3), "Jaccard",
      "Nine shared patients", patient6$analysis_id, fixed_paths[["patient6"]], "report",
      "The uniform reconstruction closely reproduced the canonical population (Jaccard 0.955)."),
    add_result("R06", "Patient 6", "Patient 6 marker-defined IFN-high", row_value(uniform, "patient6_marker_ifn_high_n"),
      "of 2004 cells", "Sensitivity analysis only", marker$analysis_id, fixed_paths[["marker"]], "report_qualified",
      "Patient 6 contributed 2,004 candidates, of which 567 met the frozen IFN-response-high rule."),
    add_result("R07", "LINCS", "Locked queries completed", 121, "queries",
      "1 primary, 10 biological sensitivity, 9 leave-one-out, 1 Hallmark, 100 matched null",
      objects$batch$analysis_id, fixed_paths[["batch"]], "report",
      "All 121 manifest-locked connectivity queries completed successfully."),
    add_result("R08", "Candidate", "Rosuvastatin primary THP1 NCS", fmt(rosu$primary_NCS, 3), "NCS",
      "PRIMARY_150; computational prediction", average_ties$analysis_id,
      run_outputs$average_ties$path, "report_qualified",
      paste0("Rosuvastatin had PRIMARY_150 THP1 NCS ", fmt(rosu$primary_NCS, 3),
             " and ranked first in the integrated robustness analysis.")),
    add_result("R09", "Candidate", "Roflumilast primary THP1 NCS", fmt(rofl$primary_NCS, 3), "NCS",
      "PRIMARY_150; computational prediction", average_ties$analysis_id,
      run_outputs$average_ties$path, "report_qualified",
      paste0("Roflumilast had PRIMARY_150 THP1 NCS ", fmt(rofl$primary_NCS, 3),
             " and ranked 15th in the integrated robustness analysis.")),
    add_result("R10", "Null", "Priority observed scores reached by matched nulls", 0, "of 100 per compound",
      "Nominal compound-specific empirical lower-tail comparison", average_ties$analysis_id,
      run_outputs$average_ties$path, "report_qualified",
      "Neither priority compound's observed PRIMARY_150 NCS was reached in 100 expression-matched null queries (nominal empirical P=1/101 for each)."),
    add_result("R11", "Rank ties", "PRIMARY_150 exact zero scores", row_value(primary_zero, "NCS_zero_n"), "of 1376 compounds",
      "Exact saved THP1 NCS; zeros are ties, not missingness", rank_audit$analysis_id,
      run_outputs$rank_audit$path, "report",
      "PRIMARY_150 contained 779 exact-zero THP1 NCS values; average ties were therefore adopted for percentile reporting."),
    add_result("R12", "Rank ties", "Top-100 membership retained", "100/100", "compounds",
      "Minimum-tie versus average-tie sensitivity", rank_audit$analysis_id,
      run_outputs$rank_audit$path, "report",
      "Average-tie ranking changed 885 ordinal positions but retained identical top-10, top-20, top-50, and top-100 membership; priority ranks were unchanged."),
    add_result("R13", "CSF", "Paired transfer-score donors", 9, "donors",
      "Exploratory paired PBMC/ventricular-CSF comparison", csf$analysis_id,
      fixed_paths[["csf"]], "supportive_only",
      "The blood-derived program was detectable in ventricular-CSF candidates, with positive paired transfer-score differences in all nine matched donors (nominal P=0.009)."),
    add_result("R14", "CSF", "At least 30 CSF candidates", csf_audit$summary$primary_donors_ge30_n,
      "of 9 paired donors", "Minimum-cell sensitivity", csf_audit$analysis_id,
      run_outputs$csf_audit$path, "supportive_only",
      "Six of nine paired donors had at least 30 CSF candidates; the restricted sensitivity remained directionally concordant (nominal P=0.036)."),
    add_result("R15", "CSF", "Largest donor contribution", fmt(csf_audit$summary$largest_donor_share_pct, 2), "%",
      paste("Selected CSF candidates; Patient", csf_audit$summary$largest_donor), csf_audit$analysis_id,
      run_outputs$csf_audit$path, "limitation",
      "CSF cell recovery was highly unequal: Patient 10 contributed 77.12% of selected candidates (effective donor number 1.64 of 10)."),
    add_result("R16", "Donor characteristics", "Publicly linked clinical covariates", 0,
      "exact patient-level age/sex/medication/MoCA values", "No imputation; MoCA group available for 9 of 10",
      donor_table$analysis_id, run_outputs$donor_table$path, "limitation",
      "Only publicly linkable fields are reported; unavailable patient-level age, sex, medication, and exact MoCA data were not inferred."),
    add_result("R17", "Exposure metadata", "Linked dose and duration fields", 0, "fields",
      "Saved executed THP1 results", exposure$analysis_id, run_outputs$exposure$path, "limitation",
      "The saved perturbagen-level results retained neither dose nor exposure duration; no unsupported compound-only linkage was performed."),
    add_result("R18", "Compound metadata", "Top-20 identities resolved", nrow(curated$curated_top20), "of 20",
      "Qualified BRD-level reporting; physical assay material not independently verified", curated$analysis_id,
      run_outputs$curated$path, "report_qualified",
      "All top-20 BRD identities were curated for qualified reporting while preserving BRD-level scores and ranks." )
  ))

  claim_guardrails <- data.frame(
    claim_id = sprintf("C%02d", 1:12),
    domain = c("State", "Specificity", "Comparator", "Candidates", "Candidates", "Null",
               "CSF", "CSF", "Cognition", "Clinical", "Exposure", "Donor data"),
    disposition = c("ALLOW", "PROHIBIT", "REQUIRE_QUALIFIER", "ALLOW_QUALIFIED",
                    "PROHIBIT", "REQUIRE_QUALIFIER", "ALLOW_SUPPORTIVE", "PROHIBIT",
                    "PROHIBIT", "PROHIBIT", "REQUIRE_LIMITATION", "REQUIRE_LIMITATION"),
    wording = c(
      "Use 'IFN-response-high monocyte state' or 'marker-defined IFN-response-high state'.",
      "Do not call the signature iNPH-specific or an iNPH-versus-control disease signature.",
      "State that contrasts are within iNPH monocytes and that SLE/healthy comparisons are contextualization.",
      "Describe rosuvastatin and roflumilast as computationally prioritized reversal candidates.",
      "Do not claim efficacy, clinical benefit, repurposing validation, or iNPH regulatory approval.",
      "Call matched-null P values nominal, compound-specific empirical probabilities; do not label them FDR.",
      "Report CSF only as supportive cross-compartment detectability/contextualization.",
      "Do not report pooled 15.21% as cohort burden or claim CSF enrichment, CNS origin, or causality.",
      "Do not infer association with cognition or MoCA from the available donor data.",
      "Do not infer missing age, sex, medication, exact MoCA, or other clinical characteristics.",
      "Explicitly state that executed-result dose and duration could not be linked.",
      "Present MoCA group for nine donors and mark Patient 6/other unavailable fields as unavailable."
    ),
    evidence_source = c("02B", "02B/Phase 2", "02B/Phase 2", "04E/04G", "04G",
                        "02D/04B", "05B/05C", "05C", "05D", "05D", "04H", "05D"),
    stringsAsFactors = FALSE
  )

  reviewer_map <- data.frame(
    map_id = sprintf("M%02d", 1:16),
    reviewer = c("R1/R2", "R1/R2", "R1/R2", "R1/R2", "R1/R2", "R1/R2",
                 "R1/R2", "R1/R2", "R1/R2", "R1/R2", "R1/R2", "R1/R2",
                 "R1/R2", "R1/R2", "R1/R2", "R1/R2"),
    topic = c("Disease specificity", "Cell-state definition", "Clustering resolution",
              "Patient 6 exclusion", "Leave-one-donor-out", "Signature-size sensitivity",
              "Null specificity", "Zero NCS and ties", "Candidate selection",
              "THP1 model limitation", "Exposure metadata", "Compound identity",
              "CSF relevance", "Cognition/MoCA", "Donor characteristics", "Claims and terminology"),
    evidence = c(
      "02B marker-defined state plus Phase 2 SLE/HC contextualization",
      "Frozen AUCell marker threshold; 3,195/22,197 canonical high cells; cluster concordance",
      "01B resolution 0.3/0.5/0.7/1.0 sensitivity",
      "01C uniform ten-patient reconstruction and 02B state transfer",
      "Nine donor-exclusion signatures and completed full-reference LINCS queries",
      "PRIMARY_50/100/150/200 completed under locked manifest",
      "100 expression-matched null signatures; compound-specific empirical probabilities",
      "04D source audit and all-compound minimum-versus-average-tie comparison",
      "04E integrated robustness plus explicit a priori/post-hoc distinction",
      "THP1 is a leukemia-derived monocytic model; computational perturbation context only",
      "04H found no retained/linkable dose or duration and performed no unsupported join",
      "04F/04G BRD-level curation; top 20 resolved; scores/ranks unchanged",
      "05B transfer plus 05C cell-yield audit and 05E transparent donor figure",
      "05D shows exact linked MoCA unavailable; no cognition analysis supported",
      "05D public-only table; no imputation",
      "06B claim guardrails and locked reportable wording"
    ),
    response_action = c(
      "Reframe throughout as within-iNPH IFN-response state; remove disease-specific language.",
      "Make marker-defined state primary and cluster 5 a concordance/sensitivity result.",
      "Report resolution sensitivity in Supplement and cite it in Methods/Results.",
      "Report Patient 6 inclusion as sensitivity, without retroactively making it primary.",
      "Report candidate NCS/rank stability across every donor exclusion.",
      "Report stability across signature sizes and retain PRIMARY_150 as submitted analysis.",
      "Report nominal matched-null evidence and avoid multiplicity-adjusted claims.",
      "Use average ties for reporting; explain exact zeros and unchanged top-k membership.",
      "Describe candidate hierarchy transparently; avoid implying roflumilast was pre-specified.",
      "Add model-system limitation and avoid direct therapeutic extrapolation.",
      "State metadata unavailability explicitly; do not invent dose/time specificity.",
      "Use curated names/qualifiers; retain BRD IDs in supplement.",
      "Use supportive Supplementary figure; prohibit prevalence/enrichment claims.",
      "Remove cognition association claims and identify unavailable exact scores.",
      "Add public donor table with unavailable fields shown explicitly.",
      "Audit abstract, results, discussion, legends, and response letter against guardrails."
    ),
    manuscript_destination = c("Title/Abstract/Introduction/Discussion", "Methods/Results",
      "Supplementary Methods/Figure/Table", "Supplementary Results/Table", "Results/Supplement",
      "Results/Supplement", "Results/Supplement", "Methods/Results/Supplement", "Results/Discussion",
      "Discussion limitations", "Methods/limitations", "Supplementary candidate table",
      "Supplementary Results/Figure", "Results/limitations", "Supplementary Table", "All sections"),
    status = "EVIDENCE_LOCKED__TEXT_EDIT_PENDING",
    stringsAsFactors = FALSE
  )

  placement <- data.frame(
    item_id = c("F_MAIN_ROBUSTNESS", "F_SUPP_HEATMAP", "F_SUPP_LOO", "F_SUPP_NULL",
                "F_SUPP_PERCENTILE", "F_SUPP_CSF", "T_SUPP_DONORS", "T_SUPP_CANDIDATES",
                "T_SUPP_SENSITIVITY", "T_SUPP_AUDIT"),
    item = c("Priority candidate robustness summary", "Top-20 connectivity heatmap",
             "Priority leave-one-donor-out stability", "Priority matched-null distributions",
             "Primary versus integrated average-tie percentiles", "CSF donor-transparency figure",
             "Public donor characteristics", "Curated candidate annotations and BRD IDs",
             "Resolution, state, signature-size, Patient 6 and LOO results", "Zero/tie and exposure audit"),
    recommended_placement = c("MAIN", rep("SUPPLEMENT", 9)),
    rationale = c(
      "Directly answers robustness of the principal candidate conclusions.",
      "Dense diagnostic overview; supports rather than replaces the candidate summary.",
      "Direct reviewer-requested sensitivity with detailed donor exclusions.",
      "Specificity diagnostic requiring nominal empirical interpretation.",
      "Technical tie-handling diagnostic; bands are mathematical consequences of tied scores.",
      "Supportive only; prominently displays unequal cell recovery and uncertainty.",
      "Publicly available fields only, with missingness explicit.",
      "Preserves BRD identity and qualified annotation details.",
      "Complete audit trail is too detailed for the main narrative.",
      "Methodological transparency and limitations."
    ),
    decision_status = "RECOMMENDED_FOR_DOCUMENT_INTEGRATION",
    stringsAsFactors = FALSE
  )

  checklist <- data.frame(
    order = 1:12,
    task = c("Freeze this 06B evidence lock in Git", "Revise title and abstract terminology",
      "Rewrite cell-state Methods", "Rewrite robustness and null Methods", "Rewrite Results",
      "Shorten and qualify Discussion", "Insert main robustness panel", "Assemble Supplement figures",
      "Assemble Supplement tables", "Draft point-by-point response letter",
      "Audit citations, ethics, data/code availability and references", "Final cross-document number audit"),
    required_input = c("06B outputs", "Claim guardrails", "02B/Phase 2", "02E/03B/04D/04E",
      "Reportable results", "Claim guardrails/limitations", "04E figures", "01B-05E figures",
      "01B-05D tables", "Reviewer evidence map", "Original manuscript and source paper",
      "Evidence lock plus final manuscript/response/supplement"),
    status = "PENDING",
    stringsAsFactors = FALSE
  )

  text_blocks <- c(
    "LOCKED INTERPRETIVE FRAME",
    "The primary finding is a marker-defined IFN-response-high monocyte state within iNPH PBMCs. It is not an iNPH-versus-control disease signature and is not claimed to be iNPH-specific.",
    "",
    "LOCKED CANDIDATE FRAME",
    paste0(rosu$pert_iname, " (PRIMARY_150 THP1 NCS ", fmt(rosu$primary_NCS, 3),
           "; integrated robustness rank 1) and ", rofl$pert_iname,
           " (NCS ", fmt(rofl$primary_NCS, 3), "; rank 15) are computationally prioritized transcriptional-reversal candidates. These results do not establish efficacy in iNPH."),
    "",
    "LOCKED NULL/TIE FRAME",
    "Neither priority observed score was reached in 100 expression-matched null queries; the plus-one lower-tail probability is 1/101 for each and is nominal, compound-specific, and not an FDR. Exact zero NCS values were treated as ties using average ranks; candidate scores and top-100 membership were unchanged.",
    "",
    "LOCKED CSF FRAME",
    "The blood-derived program was detectable in ventricular-CSF candidate monocytes, but markedly unequal cell recovery limits prevalence and enrichment inference. This is supportive cross-compartment contextualization only.",
    "",
    "LOCKED LIMITATIONS",
    "The LINCS readout is from THP1, a leukemia-derived monocytic model. Saved executed results lacked linkable dose and exposure-duration fields. Public data did not permit patient-level linkage of age, sex, medications, or exact MoCA scores; these values were not inferred."
  )

  parent <- "results/revision/evidence_claim_lock"
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  check(dir.create(out), "Cannot create the 06B output directory.")
  complete <- FALSE
  on.exit(if (!complete) message(
    "06B did not finish. Do not use a directory without _SUCCESS.txt: ", out), add = TRUE)

  write_csv <- function(x, filename) utils::write.csv(
    x, file.path(out, filename), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write_csv(reportable, "reportable_results.csv")
  write_csv(reviewer_map, "reviewer_comment_evidence_map.csv")
  write_csv(claim_guardrails, "claim_guardrails.csv")
  write_csv(placement, "figure_table_placement_plan.csv")
  write_csv(checklist, "manuscript_revision_checklist.csv")
  writeLines(text_blocks, file.path(out, "reportable_text_blocks.txt"), useBytes = TRUE)

  input_inventory <- data.frame(
    input = input_paths, md5 = input_md5_before,
    size_bytes = unname(file.info(input_paths)$size),
    modified = format(file.info(input_paths)$mtime, "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )
  write_csv(input_inventory, "input_inventory.csv")

  canonical_lock <- c(
    paste("BUILD", BUILD, sep = "\t"),
    paste("MANIFEST_SHA256", MANIFEST_SHA, sep = "\t"),
    apply(reportable, 1L, function(x) paste(c("RESULT", x), collapse = "\t")),
    apply(claim_guardrails, 1L, function(x) paste(c("CLAIM", x), collapse = "\t")),
    apply(reviewer_map, 1L, function(x) paste(c("MAP", x), collapse = "\t"))
  )
  lock_txt <- file.path(out, "revision_evidence_lock.txt")
  writeLines(canonical_lock, lock_txt, useBytes = TRUE)
  check(requireNamespace("digest", quietly = TRUE),
        "Package 'digest' is required only to calculate the evidence-lock SHA-256.")
  evidence_sha <- digest::digest(file = lock_txt, algo = "sha256")

  result <- list(
    analysis_id = "revision_06B_evidence_and_claim_lock",
    build_id = BUILD,
    status = "PASS_REVISION_EVIDENCE_AND_CLAIMS_LOCKED",
    manifest_content_sha256 = MANIFEST_SHA,
    evidence_lock_sha256 = evidence_sha,
    reportable_results = reportable,
    claim_guardrails = claim_guardrails,
    reviewer_comment_evidence_map = reviewer_map,
    figure_table_placement_plan = placement,
    manuscript_revision_checklist = checklist,
    reportable_text_blocks = text_blocks,
    selected_source_outputs = vapply(run_outputs, `[[`, character(1), "path"),
    input_provenance = input_inventory,
    new_biological_analyses_n = 0L,
    new_lincs_queries_n = 0L,
    output_directory = out,
    session_info = utils::sessionInfo()
  )
  rds_path <- file.path(out, "revision_evidence_claim_lock.rds")
  saveRDS(result, rds_path, compress = TRUE)
  check(identical(readRDS(rds_path), result), "06B RDS read-back failed.")
  check(identical(input_md5_before, unname(tools::md5sum(input_paths))),
        "An input changed while the evidence lock was being built.")

  exported <- list.files(out, full.names = TRUE)
  check(length(exported) == 9L && all(file.info(exported)$size > 0),
        "The pre-checksum 06B export set is incomplete or contains an empty file.")
  checksums <- data.frame(
    file = basename(exported), md5 = unname(tools::md5sum(exported)),
    stringsAsFactors = FALSE
  )
  write_csv(checksums, "output_checksums.csv")
  writeLines(c(
    "PASS_REVISION_EVIDENCE_AND_CLAIMS_LOCKED",
    paste("Build:", BUILD), paste("Manifest SHA-256:", MANIFEST_SHA),
    paste("Evidence lock SHA-256:", evidence_sha),
    "New biological analyses: 0", "New LINCS queries: 0"
  ), file.path(out, "_SUCCESS.txt"), useBytes = TRUE)
  complete <- TRUE

  cat("\n============================================================\n")
  cat("  REVISION ANALYSIS 06B COMPLETE\n")
  cat("============================================================\n")
  cat("  Status                    : PASS_REVISION_EVIDENCE_AND_CLAIMS_LOCKED\n")
  cat("  Reportable results locked :", nrow(reportable), "\n")
  cat("  Claim guardrails locked   :", nrow(claim_guardrails), "\n")
  cat("  Reviewer evidence mappings:", nrow(reviewer_map), "\n")
  cat("  Placement decisions       :", nrow(placement), "\n")
  cat("  Manifest SHA-256          :", MANIFEST_SHA, "\n")
  cat("  Evidence lock SHA-256     :", evidence_sha, "\n")
  cat("  New biological analyses   : 0\n")
  cat("  New LINCS queries         : 0\n")
  cat("  Output directory          :", out, "\n")
  cat("============================================================\n")
  cat("Next: edit the manuscript and point-by-point response from this locked evidence set.\n")
})
