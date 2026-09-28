# =============================================================================
# Revision analysis 05C — CSF cell-count heterogeneity and inference audit
# Install as analysis/05C_audit_csf_cell_count_heterogeneity.R
# Build: 05C_CSF_CELL_COUNT_HETEROGENEITY_AUDIT_V1_20260907
# =============================================================================
#
# Purpose
# -------
# Audit whether the ventricular-CSF contextualization from Analysis 05B is
# robust to highly unequal donor cell yields. This script uses only saved 05B
# results. It does not reload raw counts, recluster cells, rescore signatures,
# test MoCA groups, or run LINCS.
#
# Reporting objective
# -------------------
# Separate evidence that the transferred IFN-response-high program is present
# in CSF monocyte candidates from stronger claims about its prevalence or
# enrichment. Donor-level results are re-evaluated at prespecified minimum CSF
# candidate counts of 1, 10, 20, 30, 50, and 100 cells per donor.
# =============================================================================

local({
  BUILD <- "05C_CSF_CELL_COUNT_HETEROGENEITY_AUDIT_V1_20260907"

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  wilson_interval <- function(successes, trials, confidence = 0.95) {
    if (trials <= 0L) return(c(lower = NA_real_, upper = NA_real_))
    z <- stats::qnorm(1 - (1 - confidence) / 2)
    observed <- successes / trials
    denominator <- 1 + z^2 / trials
    center <- (observed + z^2 / (2 * trials)) / denominator
    half_width <- z * sqrt(
      observed * (1 - observed) / trials + z^2 / (4 * trials^2)
    ) / denominator
    c(
      lower = max(0, center - half_width),
      upper = min(1, center + half_width)
    )
  }

  paired_test_row <- function(data, csf_column, pbmc_column, metric) {
    csf <- as.numeric(data[[csf_column]])
    pbmc <- as.numeric(data[[pbmc_column]])
    difference <- csf - pbmc

    if (length(csf) >= 3L) {
      test <- suppressWarnings(stats::wilcox.test(
        csf,
        pbmc,
        paired = TRUE,
        exact = FALSE,
        alternative = "two.sided"
      ))
      statistic <- unname(test$statistic)
      p_value <- test$p.value
    } else {
      statistic <- NA_real_
      p_value <- NA_real_
    }

    data.frame(
      metric = metric,
      eligible_donors_n = length(csf),
      PBMC_median = if (length(pbmc)) stats::median(pbmc) else NA_real_,
      CSF_median = if (length(csf)) stats::median(csf) else NA_real_,
      median_paired_difference_CSF_minus_PBMC = if (length(difference)) {
        stats::median(difference)
      } else {
        NA_real_
      },
      donors_with_positive_difference_n = sum(difference > 0),
      donors_with_negative_difference_n = sum(difference < 0),
      donors_with_zero_difference_n = sum(difference == 0),
      wilcoxon_V = statistic,
      nominal_unadjusted_p = p_value,
      stringsAsFactors = FALSE
    )
  }

  run <- function() {
    started <- Sys.time()

    input_path <- "data/processed/revision_csf_ifn_contextualization.rds"
    check(file.exists(input_path), "Run Analysis 05B before Analysis 05C.")

    input_md5_before <- unname(tools::md5sum(input_path))
    check(!is.na(input_md5_before), "Cannot calculate the 05B input checksum.")

    result_05b <- readRDS(input_path)
    check(
      identical(
        result_05b$status,
        "PASS_ANALYSIS_COMPLETE_INTERPRET_FROM_EXPORTED_RESULTS"
      ),
      "The selected 05B result is not a completed contextualization result."
    )

    required_fields <- c(
      "configuration",
      "contextualization_summary",
      "csf_cell_scores",
      "csf_donor_summary",
      "canonical_pbmc_donor_summary",
      "paired_donor_summary",
      "paired_tests",
      "cluster_diagnostics",
      "isg_detection",
      "input_provenance"
    )
    check(
      all(required_fields %in% names(result_05b)),
      "The 05B result is missing one or more required fields."
    )

    cells <- as.data.frame(result_05b$csf_cell_scores)
    csf_donors <- as.data.frame(result_05b$csf_donor_summary)
    paired <- as.data.frame(result_05b$paired_donor_summary)
    primary_patients <- as.character(
      result_05b$configuration$primary_matched_patients
    )

    required_cell_columns <- c(
      "cell_id",
      "patient",
      "csf_cluster",
      "transfer_score_exact",
      "transfer_score_selection_independent",
      "marker_ifn_score",
      "marker_ifn_state"
    )
    required_donor_columns <- c(
      "patient",
      "n_cells",
      "marker_ifn_high_n",
      "marker_ifn_high_pct",
      "median_transfer_score",
      "median_marker_ifn_score"
    )
    required_paired_columns <- c(
      "patient",
      "n_cells_PBMC",
      "n_cells_CSF",
      "median_transfer_score_PBMC",
      "median_transfer_score_CSF",
      "median_marker_ifn_score_PBMC",
      "median_marker_ifn_score_CSF",
      "marker_ifn_high_pct_PBMC",
      "marker_ifn_high_pct_CSF"
    )

    check(
      all(required_cell_columns %in% names(cells)) &&
        all(required_donor_columns %in% names(csf_donors)) &&
        all(required_paired_columns %in% names(paired)),
      "Unexpected 05B cell, donor, or paired-summary schema."
    )

    check(
      nrow(cells) == sum(csf_donors$n_cells) &&
        !anyDuplicated(cells$cell_id) &&
        nrow(csf_donors) == 10L &&
        nrow(paired) == 9L &&
        setequal(as.character(paired$patient), primary_patients),
      "The 05B cell and donor universes do not reconcile."
    )

    # Add uncertainty and donor-contribution diagnostics.
    intervals <- t(vapply(
      seq_len(nrow(csf_donors)),
      function(i) wilson_interval(
        csf_donors$marker_ifn_high_n[i],
        csf_donors$n_cells[i]
      ),
      numeric(2)
    ))

    donor_audit <- csf_donors
    donor_audit$selected_cell_share_pct <-
      100 * donor_audit$n_cells / sum(donor_audit$n_cells)
    donor_audit$marker_ifn_high_wilson_lower_pct <- 100 * intervals[, "lower"]
    donor_audit$marker_ifn_high_wilson_upper_pct <- 100 * intervals[, "upper"]
    donor_audit$meets_10_cell_minimum <- donor_audit$n_cells >= 10L
    donor_audit$meets_20_cell_minimum <- donor_audit$n_cells >= 20L
    donor_audit$meets_30_cell_minimum <- donor_audit$n_cells >= 30L
    donor_audit$meets_50_cell_minimum <- donor_audit$n_cells >= 50L
    donor_audit$meets_100_cell_minimum <- donor_audit$n_cells >= 100L
    donor_audit$is_primary_nine_patient_donor <-
      as.character(donor_audit$patient) %in% primary_patients

    donor_audit <- donor_audit[
      order(as.numeric(as.character(donor_audit$patient))),
      ,
      drop = FALSE
    ]

    cell_shares <- donor_audit$n_cells / sum(donor_audit$n_cells)
    effective_donor_number <- 1 / sum(cell_shares^2)
    largest_donor_index <- which.max(donor_audit$n_cells)

    # Quantify whether either selected cluster is donor-specific.
    cluster_counts <- as.data.frame(
      table(
        patient = as.character(cells$patient),
        cluster = as.character(cells$csf_cluster)
      ),
      stringsAsFactors = FALSE
    )
    names(cluster_counts)[3L] <- "n_cells"
    cluster_counts <- cluster_counts[cluster_counts$n_cells > 0L, , drop = FALSE]

    cluster_totals <- aggregate(
      n_cells ~ cluster,
      data = cluster_counts,
      FUN = sum
    )
    names(cluster_totals)[2L] <- "cluster_total_n"
    cluster_counts <- merge(
      cluster_counts,
      cluster_totals,
      by = "cluster",
      all.x = TRUE,
      sort = FALSE
    )
    cluster_counts$within_cluster_donor_share_pct <-
      100 * cluster_counts$n_cells / cluster_counts$cluster_total_n

    cluster_balance <- do.call(rbind, lapply(
      unique(as.character(cluster_counts$cluster)),
      function(cluster_id) {
        d <- cluster_counts[
          as.character(cluster_counts$cluster) == cluster_id,
          ,
          drop = FALSE
        ]
        shares <- d$n_cells / sum(d$n_cells)
        maximum_index <- which.max(d$n_cells)
        data.frame(
          cluster = cluster_id,
          cluster_cells_n = sum(d$n_cells),
          contributing_donors_n = nrow(d),
          largest_contributing_donor = as.character(d$patient[maximum_index]),
          largest_donor_share_pct = 100 * max(shares),
          effective_donor_number = 1 / sum(shares^2),
          stringsAsFactors = FALSE
        )
      }
    ))

    donor_specific_clusters <- as.character(
      cluster_balance$cluster[
        cluster_balance$contributing_donors_n == 1L
      ]
    )
    shared_clusters <- setdiff(
      as.character(cluster_balance$cluster),
      donor_specific_clusters
    )

    check(
      length(shared_clusters) > 0L,
      "Every selected CSF cluster is donor-specific; supportive transfer analysis is not possible."
    )

    shared_cells <- cells[
      as.character(cells$csf_cluster) %in% shared_clusters,
      ,
      drop = FALSE
    ]
    shared_patients <- unique(as.character(shared_cells$patient))
    shared_patients <- shared_patients[order(as.numeric(shared_patients))]

    shared_cluster_donor_summary <- do.call(rbind, lapply(
      shared_patients,
      function(patient_id) {
        d <- shared_cells[
          as.character(shared_cells$patient) == patient_id,
          ,
          drop = FALSE
        ]
        high_state <- d$marker_ifn_state == "marker_ifn_high"
        data.frame(
          patient = patient_id,
          shared_cluster_CSF_cells_n = nrow(d),
          shared_cluster_marker_ifn_high_n = sum(high_state),
          shared_cluster_marker_ifn_high_pct = 100 * mean(high_state),
          shared_cluster_median_transfer_score = stats::median(
            d$transfer_score_exact
          ),
          shared_cluster_median_marker_ifn_score = stats::median(
            d$marker_ifn_score
          ),
          stringsAsFactors = FALSE
        )
      }
    ))

    shared_cluster_paired <- merge(
      paired[, c(
        "patient",
        "n_cells_PBMC",
        "median_transfer_score_PBMC",
        "median_marker_ifn_score_PBMC",
        "marker_ifn_high_pct_PBMC"
      )],
      shared_cluster_donor_summary,
      by = "patient",
      all = FALSE,
      sort = FALSE
    )
    shared_cluster_paired <- shared_cluster_paired[
      match(
        intersect(primary_patients, shared_cluster_paired$patient),
        shared_cluster_paired$patient
      ),
      ,
      drop = FALSE
    ]

    shared_cluster_paired_for_test <- data.frame(
      median_transfer_score_CSF =
        shared_cluster_paired$shared_cluster_median_transfer_score,
      median_transfer_score_PBMC =
        shared_cluster_paired$median_transfer_score_PBMC,
      median_marker_ifn_score_CSF =
        shared_cluster_paired$shared_cluster_median_marker_ifn_score,
      median_marker_ifn_score_PBMC =
        shared_cluster_paired$median_marker_ifn_score_PBMC,
      marker_ifn_high_pct_CSF =
        shared_cluster_paired$shared_cluster_marker_ifn_high_pct,
      marker_ifn_high_pct_PBMC =
        shared_cluster_paired$marker_ifn_high_pct_PBMC,
      stringsAsFactors = FALSE
    )

    shared_cluster_tests <- rbind(
      paired_test_row(
        shared_cluster_paired_for_test,
        "median_transfer_score_CSF",
        "median_transfer_score_PBMC",
        "median blood-derived transfer score"
      ),
      paired_test_row(
        shared_cluster_paired_for_test,
        "median_marker_ifn_score_CSF",
        "median_marker_ifn_score_PBMC",
        "median frozen IFN-selector score"
      ),
      paired_test_row(
        shared_cluster_paired_for_test,
        "marker_ifn_high_pct_CSF",
        "marker_ifn_high_pct_PBMC",
        "percentage above frozen IFN threshold"
      )
    )
    shared_cluster_tests$included_CSF_clusters <- paste(
      shared_clusters,
      collapse = " | "
    )
    shared_cluster_tests$excluded_donor_specific_CSF_clusters <- paste(
      donor_specific_clusters,
      collapse = " | "
    )
    shared_cluster_tests$interpretation_scope <- paste(
      "sensitivity excluding selected clusters contributed by only one donor;",
      "nominal P values remain descriptive"
    )

    # Repeat the paired analyses at prespecified minimum CSF cell counts.
    minimums <- c(1L, 10L, 20L, 30L, 50L, 100L)
    sensitivity_rows <- list()
    row_index <- 1L

    for (minimum in minimums) {
      eligible <- paired[
        paired$n_cells_CSF >= minimum,
        ,
        drop = FALSE
      ]

      metric_rows <- rbind(
        paired_test_row(
          eligible,
          "median_transfer_score_CSF",
          "median_transfer_score_PBMC",
          "median blood-derived transfer score"
        ),
        paired_test_row(
          eligible,
          "median_marker_ifn_score_CSF",
          "median_marker_ifn_score_PBMC",
          "median frozen IFN-selector score"
        ),
        paired_test_row(
          eligible,
          "marker_ifn_high_pct_CSF",
          "marker_ifn_high_pct_PBMC",
          "percentage above frozen IFN threshold"
        )
      )
      metric_rows$minimum_CSF_candidates_per_donor <- minimum
      metric_rows$eligible_patients <- paste(
        eligible$patient,
        collapse = " | "
      )
      metric_rows$nominal_p_is_descriptive_only <- TRUE
      metric_rows$interpretation_scope <- paste(
        "minimum-cell sensitivity for an exploratory compartment comparison;",
        "not a disease-specificity, cognition, or causal test"
      )
      sensitivity_rows[[row_index]] <- metric_rows
      row_index <- row_index + 1L
    }

    minimum_cell_sensitivity <- do.call(rbind, sensitivity_rows)
    minimum_cell_sensitivity <- minimum_cell_sensitivity[, c(
      "minimum_CSF_candidates_per_donor",
      "metric",
      "eligible_donors_n",
      "eligible_patients",
      "PBMC_median",
      "CSF_median",
      "median_paired_difference_CSF_minus_PBMC",
      "donors_with_positive_difference_n",
      "donors_with_negative_difference_n",
      "donors_with_zero_difference_n",
      "wilcoxon_V",
      "nominal_unadjusted_p",
      "nominal_p_is_descriptive_only",
      "interpretation_scope"
    )]

    # Donor-balanced summaries avoid the 2,069-cell donor dominating the cohort.
    primary_donor_audit <- donor_audit[
      donor_audit$is_primary_nine_patient_donor,
      ,
      drop = FALSE
    ]
    primary_ge30 <- primary_donor_audit[
      primary_donor_audit$n_cells >= 30L,
      ,
      drop = FALSE
    ]

    summarize_group <- function(d, label) {
      data.frame(
        donor_group = label,
        donors_n = nrow(d),
        total_selected_cells_n = sum(d$n_cells),
        minimum_cells_per_donor = min(d$n_cells),
        maximum_cells_per_donor = max(d$n_cells),
        donor_balanced_median_IFN_high_pct = stats::median(
          d$marker_ifn_high_pct
        ),
        donor_balanced_IQR_IFN_high_pct = stats::IQR(
          d$marker_ifn_high_pct
        ),
        donor_balanced_median_transfer_score = stats::median(
          d$median_transfer_score
        ),
        donors_with_at_least_1_IFN_high_cell_n = sum(
          d$marker_ifn_high_n >= 1L
        ),
        donors_with_at_least_10_IFN_high_cells_n = sum(
          d$marker_ifn_high_n >= 10L
        ),
        stringsAsFactors = FALSE
      )
    }

    donor_balanced_summary <- rbind(
      summarize_group(donor_audit, "all_10_CSF_donors"),
      summarize_group(primary_donor_audit, "primary_9_matched_donors"),
      summarize_group(
        primary_ge30,
        "primary_donors_with_at_least_30_CSF_candidates"
      )
    )

    primary_ge30_n <- nrow(primary_ge30)
    all_primary_ge30 <- primary_ge30_n == length(primary_patients)
    largest_donor_share_pct <- donor_audit$selected_cell_share_pct[
      largest_donor_index
    ]

    inference_audit <- data.frame(
      audit_item = c(
        "All primary donors have at least 30 CSF candidates",
        "At least two thirds of primary donors have at least 30 CSF candidates",
        "Largest donor contributes no more than 50% of selected cells",
        "All ten donors contain at least one frozen-threshold IFN-high cell",
        "Transfer score is stable after removing selection-panel overlaps",
        "Every selected CSF cluster is represented by at least two donors",
        "Cell-weighted CSF percentage suitable as cohort burden estimate",
        "CSF presence/detectability suitable for supportive reporting",
        "CSF enrichment suitable for definitive reporting"
      ),
      observed = c(
        paste0(primary_ge30_n, "/", length(primary_patients)),
        paste0(primary_ge30_n, "/", length(primary_patients)),
        sprintf("%.2f%% from Patient %s", largest_donor_share_pct,
                donor_audit$patient[largest_donor_index]),
        paste0(sum(donor_audit$marker_ifn_high_n > 0L), "/10"),
        sprintf(
          "Spearman rho %.4f",
          result_05b$contextualization_summary$
            transfer_score_independence_spearman
        ),
        if (length(donor_specific_clusters)) {
          paste(
            "donor-specific selected cluster(s):",
            paste(donor_specific_clusters, collapse = " | ")
          )
        } else {
          "all selected clusters represented by multiple donors"
        },
        sprintf(
          "effective donor number %.2f of 10",
          effective_donor_number
        ),
        paste0(
          sum(donor_audit$marker_ifn_high_n > 0L),
          "/10 donors with at least one high cell"
        ),
        paste0(
          "minimum-cell sensitivity exported across thresholds; ",
          "compartment and sampling limitations remain"
        )
      ),
      pass = c(
        all_primary_ge30,
        primary_ge30_n >= ceiling(2 * length(primary_patients) / 3),
        largest_donor_share_pct <= 50,
        all(donor_audit$marker_ifn_high_n > 0L),
        result_05b$contextualization_summary$
          transfer_score_independence_spearman >= 0.95,
        length(donor_specific_clusters) == 0L,
        FALSE,
        TRUE,
        FALSE
      ),
      consequence = c(
        "Do not treat all nine paired donor estimates as equally precise",
        "A restricted donor-level sensitivity analysis is feasible",
        "Do not report the pooled 15.21% as a donor-balanced cohort burden",
        "Supports presence, subject to marker-selection and low-count caveats",
        "One selection-panel overlap does not explain the transfer signal",
        "Exclude donor-specific selected clusters in a dedicated sensitivity analysis",
        "Use donor-balanced summaries and show individual donor counts",
        "May be reported as exploratory cross-compartment contextualization",
        "Do not claim CSF enrichment or CNS origin"
      ),
      stringsAsFactors = FALSE
    )

    final_status <- paste(
      "SUPPORTIVE_CSF_DETECTABILITY_ONLY",
      "CELL_YIELD_HETEROGENEITY_LIMITS_PREVALENCE_AND_ENRICHMENT_INFERENCE",
      sep = "__"
    )

    donor_specific_description <- if (length(donor_specific_clusters)) {
      details <- vapply(
        donor_specific_clusters,
        function(cluster_id) {
          d <- cluster_balance[
            as.character(cluster_balance$cluster) == cluster_id,
            ,
            drop = FALSE
          ]
          paste0(
            "cluster ", cluster_id, " (", d$cluster_cells_n,
            " cells, all from Patient ", d$largest_contributing_donor, ")"
          )
        },
        character(1)
      )
      paste0(
        " In addition, the marker rule selected a donor-specific cluster: ",
        paste(details, collapse = "; "),
        ". A separate sensitivity analysis therefore excludes that cluster."
      )
    } else {
      " No selected cluster was contributed by only one donor."
    }

    reviewer_response <- paste0(
      "In response to the reviewer, we transferred the frozen blood-derived ",
      "IFN-response-high scoring framework to independently selected ",
      "ventricular-CSF classical-monocyte candidates from GSE292141. The ",
      "program was detectable across donors; however, CSF monocyte yields were ",
      "highly unequal (range ", min(donor_audit$n_cells), "–",
      max(donor_audit$n_cells), " cells per donor), only ", primary_ge30_n,
      " of ", length(primary_patients),
      " primary matched donors contributed at least 30 candidate cells, and ",
      "Patient ", donor_audit$patient[largest_donor_index], " contributed ",
      sprintf("%.1f%%", largest_donor_share_pct),
      " of all selected CSF candidates.", donor_specific_description,
      " We therefore report this analysis as ",
      "supportive evidence that a corresponding IFN-response-high program can ",
      "be observed at the ventricular-CSF interface, not as evidence of CSF ",
      "enrichment, CNS origin, disease specificity, or association with ",
      "cognition. Donor-level values and prespecified minimum-cell sensitivity ",
      "analyses are provided, and the pooled cell-weighted percentage is not ",
      "used as a cohort-level burden estimate."
    )

    manuscript_result <- paste0(
      "Exploratory transfer of the frozen blood-derived IFN-response-high ",
      "scoring framework identified threshold-positive cells among ",
      "ventricular-CSF classical-monocyte candidates. CSF candidate recovery ",
      "was highly heterogeneous across donors (", min(donor_audit$n_cells),
      "–", max(donor_audit$n_cells), " cells), with only ", primary_ge30_n,
      " of ", length(primary_patients),
      " primary matched donors contributing at least 30 cells.",
      donor_specific_description,
      " Therefore, the CSF analysis was interpreted as supportive evidence of ",
      "cross-compartment detectability rather than a quantitative test of ",
      "enrichment."
    )

    manuscript_limitation <- paste0(
      "The ventricular-CSF analysis was constrained by markedly unequal and, ",
      "for some donors, very low monocyte recovery. These data support only ",
      "exploratory detection of a corresponding transcriptional program and ",
      "cannot establish its CSF prevalence, enrichment, CNS origin, ",
      "relationship to cognition, or disease specificity."
    )

    immediate_inputs <- input_path
    own_script <- "analysis/05C_audit_csf_cell_count_heterogeneity.R"
    if (file.exists(own_script)) immediate_inputs <- c(immediate_inputs, own_script)
    immediate_md5_before <- unname(tools::md5sum(immediate_inputs))
    check(!anyNA(immediate_md5_before), "Cannot hash an immediate 05C input.")

    parent <- "results/revision/csf_cell_count_heterogeneity_audit"
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    output_directory <- tempfile(
      paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
      tmpdir = parent
    )
    check(dir.create(output_directory), "Cannot create the 05C output directory.")

    result <- list(
      analysis_id = "revision_05C_csf_cell_count_heterogeneity_audit",
      build_id = BUILD,
      status = final_status,
      source_05B = input_path,
      source_05B_md5 = input_md5_before,
      donor_audit = donor_audit,
      donor_balanced_summary = donor_balanced_summary,
      selected_cluster_donor_counts = cluster_counts,
      selected_cluster_balance = cluster_balance,
      shared_cluster_donor_summary = shared_cluster_donor_summary,
      shared_cluster_paired_summary = shared_cluster_paired,
      shared_cluster_tests = shared_cluster_tests,
      minimum_cell_sensitivity = minimum_cell_sensitivity,
      inference_audit = inference_audit,
      summary = list(
        selected_cells_n = nrow(cells),
        primary_donors_n = length(primary_patients),
        primary_donors_ge30_n = primary_ge30_n,
        largest_donor = as.character(donor_audit$patient[largest_donor_index]),
        largest_donor_share_pct = largest_donor_share_pct,
        effective_donor_number = effective_donor_number,
        donor_specific_selected_clusters = donor_specific_clusters,
        shared_selected_clusters = shared_clusters,
        reporting_decision = final_status,
        pooled_cell_weighted_percentage_for_cohort_reporting = FALSE,
        definitive_CSF_enrichment_claim_supported = FALSE,
        supportive_CSF_detectability_statement_supported = TRUE
      ),
      manuscript_text = list(
        result = manuscript_result,
        limitation = manuscript_limitation,
        response_to_reviewer = reviewer_response
      ),
      input_provenance = data.frame(
        input = immediate_inputs,
        md5 = immediate_md5_before,
        stringsAsFactors = FALSE
      ),
      output_directory = output_directory,
      session_info = utils::sessionInfo()
    )

    utils::write.csv(
      donor_audit,
      file.path(output_directory, "csf_donor_cell_count_and_uncertainty.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      donor_balanced_summary,
      file.path(output_directory, "csf_donor_balanced_summary.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      cluster_counts,
      file.path(output_directory, "csf_selected_cluster_donor_counts.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      cluster_balance,
      file.path(output_directory, "csf_selected_cluster_balance.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      shared_cluster_donor_summary,
      file.path(output_directory, "csf_shared_cluster_donor_summary.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      shared_cluster_paired,
      file.path(output_directory, "csf_shared_cluster_paired_summary.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      shared_cluster_tests,
      file.path(output_directory, "csf_shared_cluster_paired_tests.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      minimum_cell_sensitivity,
      file.path(output_directory, "csf_minimum_cell_sensitivity.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      inference_audit,
      file.path(output_directory, "csf_inference_audit.csv"),
      row.names = FALSE,
      na = ""
    )
    writeLines(
      c(
        "MANUSCRIPT RESULT",
        manuscript_result,
        "",
        "MANUSCRIPT LIMITATION",
        manuscript_limitation,
        "",
        "RESPONSE TO REVIEWER",
        reviewer_response
      ),
      file.path(output_directory, "csf_reportable_text.txt")
    )

    result_path <- file.path(
      output_directory,
      "revision_csf_cell_count_heterogeneity_audit.rds"
    )
    saveRDS(result, result_path)
    check(identical(readRDS(result_path), result), "05C RDS read-back failed.")

    check(
      identical(
        immediate_md5_before,
        unname(tools::md5sum(immediate_inputs))
      ),
      "An immediate input changed during Analysis 05C."
    )

    exports <- list.files(output_directory, full.names = TRUE)
    check(all(file.info(exports)$size > 0), "An Analysis 05C export is empty.")
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
      c(final_status, BUILD),
      file.path(output_directory, "_SUCCESS.txt")
    )

    cat("\n============================================================\n")
    cat("  REVISION ANALYSIS 05C COMPLETE\n")
    cat("============================================================\n")
    cat("  Status                         :", final_status, "\n")
    cat("  Selected CSF candidates        :", nrow(cells), "\n")
    cat("  Cell-count range per donor     :", min(donor_audit$n_cells), "to",
        max(donor_audit$n_cells), "\n")
    cat("  Primary donors with >=30 cells :", primary_ge30_n, "/",
        length(primary_patients), "\n")
    cat("  Largest donor contribution     : Patient",
        donor_audit$patient[largest_donor_index], sprintf("(%.2f%%)", largest_donor_share_pct), "\n")
    cat("  Effective donor number         :", round(effective_donor_number, 2), "/ 10\n")
    cat("  Donor-specific selected cluster:",
        if (length(donor_specific_clusters)) paste(donor_specific_clusters, collapse = ", ") else "none",
        "\n")
    cat("  Pooled 15.21% as cohort burden : DO NOT REPORT\n")
    cat("  Supportive CSF detectability   : REPORT WITH LIMITATIONS\n")
    cat("  Definitive CSF enrichment      : NOT SUPPORTED\n")
    cat("\nMinimum-cell sensitivity for transfer score:\n")
    print(
      minimum_cell_sensitivity[
        minimum_cell_sensitivity$metric ==
          "median blood-derived transfer score",
        c(
          "minimum_CSF_candidates_per_donor",
          "eligible_donors_n",
          "eligible_patients",
          "median_paired_difference_CSF_minus_PBMC",
          "donors_with_positive_difference_n",
          "nominal_unadjusted_p"
        )
      ],
      row.names = FALSE,
      digits = 4
    )
    cat("\nSelected-cluster donor balance:\n")
    print(cluster_balance, row.names = FALSE, digits = 4)
    cat("\nSensitivity excluding donor-specific selected cluster(s):\n")
    print(shared_cluster_tests, row.names = FALSE, digits = 4)
    cat("\nDonor-balanced summaries:\n")
    print(donor_balanced_summary, row.names = FALSE, digits = 4)
    cat("\nOutput directory:", output_directory, "\n")
    cat("Elapsed minutes:",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), "\n")
    cat("============================================================\n")

    invisible(result)
  }

  run()
})
