# Reporting exports from saved CSF score tables. Run from the project root.
# No clustering, cell selection, AUCell scoring or LINCS queries are rerun.
# Original and selection-gene-removed scores are summarized separately.
# The cluster-13 exclusion is a new targeted annotation-sensitivity summary.

local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  source_dir <- file.path(root, "results/revision/csf_ifn_contextualization")
  inputs <- file.path(source_dir, c("csf_ifn_cell_scores.csv",
    "csf_pbmc_paired_donor_summary.csv", "csf_classical_monocyte_cluster_diagnostics.csv",
    "csf_contextualization_reporting_notes.txt"))
  if (!all(file.exists(inputs))) stop("Missing saved CSF input(s): ",
                                      paste(inputs[!file.exists(inputs)], collapse = ", "))
  read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
    check.names = FALSE, na.strings = c("", "NA"))
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE,
                                                na = "", fileEncoding = "UTF-8")
  cells <- read_csv(inputs[1]); paired <- read_csv(inputs[2])
  diagnostics <- read_csv(inputs[3]); notes <- readLines(inputs[4], warn = FALSE)
  required <- c("cell_id", "patient", "csf_cluster", "transfer_score_exact",
                "transfer_score_selection_independent")
  if (!all(required %in% names(cells)) || anyDuplicated(cells$cell_id)) {
    stop("Unexpected cell-score schema or duplicated cell identifiers.")
  }
  if (!all(c("patient", "n_cells_CSF", "median_transfer_score_CSF",
             "median_transfer_score_PBMC") %in% names(paired)) ||
      anyDuplicated(paired$patient)) stop("Unexpected paired-donor schema.")
  if (any(!is.finite(cells$transfer_score_exact)) ||
      any(!is.finite(cells$transfer_score_selection_independent))) {
    stop("Non-finite cell scores require inspection.")
  }
  groups <- split(cells, cells$patient)
  donor_scores <- do.call(rbind, lapply(groups, function(d) data.frame(
    patient = as.character(d$patient[1]), n_cells = nrow(d),
    original_median = stats::median(d$transfer_score_exact),
    selection_removed_median = stats::median(d$transfer_score_selection_independent),
    median_cellwise_change = stats::median(d$transfer_score_selection_independent -
                                           d$transfer_score_exact),
    max_absolute_cellwise_change = max(abs(d$transfer_score_selection_independent -
                                           d$transfer_score_exact)),
    stringsAsFactors = FALSE)))
  donor_scores$change_in_donor_median <- donor_scores$selection_removed_median -
    donor_scores$original_median
  donor_scores <- donor_scores[order(as.numeric(donor_scores$patient)), ]
  idx <- match(as.character(paired$patient), donor_scores$patient)
  if (anyNA(idx) || !all(donor_scores$n_cells[idx] == paired$n_cells_CSF) ||
      max(abs(donor_scores$original_median[idx] - paired$median_transfer_score_CSF)) > 1e-10) {
    stop("Saved paired summaries do not agree with the supplied cell scores.")
  }
  extract_overlap <- function(prefix) {
    line <- notes[startsWith(notes, prefix)]
    if (length(line) != 1L) stop("Selection-overlap information is missing or ambiguous.")
    value <- trimws(substring(line, nchar(prefix) + 1L))
    if (!nzchar(value)) return(character(0))
    trimws(strsplit(value, "|", fixed = TRUE)[[1]])
  }
  removed_up <- extract_overlap("Selection-panel overlap with transfer UP:")
  removed_down <- extract_overlap("Selection-panel overlap with transfer DOWN:")
  score_summary <- data.frame(n_cells = nrow(cells),
    removed_up_genes = paste(removed_up, collapse = " | "),
    removed_down_genes = paste(removed_down, collapse = " | "),
    spearman_original_vs_selection_removed = stats::cor(cells$transfer_score_exact,
      cells$transfer_score_selection_independent, method = "spearman"),
    max_absolute_cell_score_change = max(abs(cells$transfer_score_selection_independent -
                                            cells$transfer_score_exact)),
    scope = "CSF score sensitivity; no comparison of modified CSF scores with unmodified PBMC scores")
  summarize_cells <- function(x, label) {
    do.call(rbind, lapply(split(x, x$patient), function(d) data.frame(
      analysis = label, patient = as.character(d$patient[1]), n_cells_CSF = nrow(d),
      median_transfer_score_CSF = stats::median(d$transfer_score_exact),
      stringsAsFactors = FALSE)))
  }
  per_donor <- summarize_cells(cells, "original_selected_clusters")
  if ("13" %in% as.character(cells$csf_cluster)) {
    per_donor <- rbind(per_donor,
      summarize_cells(cells[as.character(cells$csf_cluster) != "13", ],
                      "exclude_cluster_13_annotation_sensitivity"))
  }
  per_donor <- merge(per_donor,
    data.frame(patient = as.character(paired$patient),
      median_transfer_score_PBMC = paired$median_transfer_score_PBMC),
    by = "patient", all = FALSE, sort = FALSE)
  per_donor$paired_difference_CSF_minus_PBMC <-
    per_donor$median_transfer_score_CSF - per_donor$median_transfer_score_PBMC
  per_donor$meets_30_CSF_cell_criterion <- per_donor$n_cells_CSF >= 30
  tests <- list(); k <- 0L
  for (label in unique(per_donor$analysis)) {
    for (minimum in c(1L, 30L)) {
      d <- per_donor[per_donor$analysis == label & per_donor$n_cells_CSF >= minimum, ]
      d <- d[order(as.numeric(d$patient)), ]
      delta <- d$paired_difference_CSF_minus_PBMC
      if (nrow(d) < 3L) stop("Too few paired donors for the planned summary.")
      test <- stats::wilcox.test(d$median_transfer_score_CSF,
        d$median_transfer_score_PBMC, paired = TRUE, alternative = "two.sided",
        exact = FALSE, correct = TRUE)
      k <- k + 1L
      tests[[k]] <- data.frame(analysis = label, minimum_CSF_candidates = minimum,
        donor_pairs_n = nrow(d), donors = paste(d$patient, collapse = " | "),
        median_paired_difference_CSF_minus_PBMC = stats::median(delta),
        minimum_paired_difference = min(delta), maximum_paired_difference = max(delta),
        positive_differences_n = sum(delta > 0), negative_differences_n = sum(delta < 0),
        zero_differences_n = sum(delta == 0), wilcoxon_V = unname(test$statistic),
        nominal_unadjusted_p = test$p.value, exact = FALSE, continuity_correction = TRUE,
        test_method = test$method, interpretation = "Exploratory detectability context only")
    }
  }
  tests <- do.call(rbind, tests)
  parent <- file.path(root, "results/revision/csf_reporting_details")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  if (!dir.create(out)) stop("Cannot create reporting output directory.")
  write_csv(score_summary, file.path(out, "csf_selection_gene_removal_summary.csv"))
  write_csv(donor_scores, file.path(out, "csf_selection_gene_removal_by_donor.csv"))
  write_csv(per_donor, file.path(out, "csf_paired_differences_and_cluster_sensitivity.csv"))
  write_csv(tests, file.path(out, "csf_paired_tests_and_cluster_sensitivity.csv"))
  write_csv(as.data.frame(table(patient = cells$patient, cluster = cells$csf_cluster)),
            file.path(out, "csf_selected_cluster_donor_counts.csv"))
  write_csv(diagnostics, file.path(out, "csf_cluster_annotation_diagnostics.csv"))
  state_path <- file.path(root, "data/processed/revision_csf_ifn_contextualization.rds")
  configuration_exported <- FALSE
  if (file.exists(state_path)) {
    state <- readRDS(state_path)
    cfg <- state$configuration
    gene_fields <- c("transfer_up_exact", "transfer_down_exact",
                     "transfer_up_selection_independent", "transfer_down_selection_independent")
    if (!all(gene_fields %in% names(cfg))) stop("Saved 05B gene definitions are incomplete.")
    if (!setequal(as.character(state$csf_cell_scores$cell_id), as.character(cells$cell_id))) {
      stop("Saved 05B configuration belongs to a different set of cells.")
    }
    saved_order <- match(cells$cell_id, state$csf_cell_scores$cell_id)
    if (max(abs(state$csf_cell_scores$transfer_score_exact[saved_order] -
                cells$transfer_score_exact)) > 1e-10 ||
        max(abs(state$csf_cell_scores$transfer_score_selection_independent[saved_order] -
                cells$transfer_score_selection_independent)) > 1e-10) {
      stop("Saved 05B scores differ from the supplied CSV.")
    }
    genes <- do.call(rbind, lapply(gene_fields, function(field) data.frame(
      gene_set = field, gene = as.character(cfg[[field]]), stringsAsFactors = FALSE)))
    write_csv(genes, file.path(out, "csf_transferred_signature_gene_sets.csv"))
    write_csv(data.frame(gene_set = gene_fields,
      gene_count = vapply(gene_fields, function(f) length(cfg[[f]]), integer(1))),
      file.path(out, "csf_transferred_signature_gene_counts.csv"))
    if (!is.null(state$session_info)) capture.output(state$session_info,
      file = file.path(out, "saved_05B_execution_session_info.txt"))
    inputs <- c(inputs, state_path)
    configuration_exported <- TRUE
  }
  write_csv(data.frame(saved_05B_configuration_exported = configuration_exported,
    remaining_gene_definition_gap = !configuration_exported),
    file.path(out, "csf_configuration_export_status.csv"))
  write_csv(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs))),
            file.path(out, "input_provenance.csv"))
  capture.output(utils::sessionInfo(), file = file.path(out, "csf_reporting_session_info.txt"))
  writeLines(c(
    "Original transfer-score summaries were checked against the saved paired donor table.",
    "Selection-gene-removal scores were already computed in 05B; only their reporting summaries were exported here.",
    "High Spearman correlation describes score ordering; it does not imply identical score magnitudes.",
    "Modified CSF scores were not compared with unmodified PBMC scores.",
    "Cluster-13 exclusion is a targeted additional sensitivity summary from saved cells, motivated by its mixed lineage-panel expression.",
    "Aggregate annotation-panel scores do not prove the identity or purity of every selected cell.",
    "Paired tests use a two-sided normal approximation with continuity correction, exact=FALSE and correct=TRUE.",
    "P values are nominal and exploratory; results support detectability and do not establish compartmental enrichment."
  ), file.path(out, "csf_reporting_notes.txt"))
  writeLines(if (configuration_exported) "PASS_REPORTING_EXPORTS" else
    "REPORTING_EXPORTS_COMPLETE; SAVED_05B_GENE_CONFIGURATION_NOT_AVAILABLE",
    file.path(out, "_SUCCESS.txt"))
  cat("\nCSF reporting exports completed:\n", out, "\n", sep = "")
  print(tests, row.names = FALSE)
})
