# =============================================================================
# Revision analysis 06C — Export locked manuscript figures and tables
# Install as: analysis/06C_export_locked_manuscript_assets.R
# Build: 06C_LOCKED_MANUSCRIPT_ASSET_EXPORT_V1_20260908
# =============================================================================
#
# Purpose
# -------
# This is an export/assembly script. It does not rerun clustering, differential
# expression, AUCell scoring, LINCS connectivity analysis, or statistical tests.
# It reads the completed reviewer-revision outputs, validates the frozen 06B
# evidence files, and creates publication-ready standalone panels and compact
# tables for manual assembly in PowerPoint.
#
# The script deliberately leaves analysis/05_generate_all_figures.R unchanged so
# the originally submitted figure-generation workflow remains reproducible.
# =============================================================================

local({
  BUILD <- "06C_LOCKED_MANUSCRIPT_ASSET_EXPORT_V1_20260908"
  EXPECTED_EVIDENCE_LOCK_SHA256 <-
    "21e6225b9e86bda29896f11a55e61d1618e9ce5a57dc6d271a26e99e358abf95"

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  required_packages <- c(
    "ggplot2", "patchwork", "scales", "digest", "ggrepel"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  check(
    length(missing_packages) == 0L,
    paste(
      "Install required package(s) before running 06C:",
      paste(missing_packages, collapse = ", ")
    )
  )

  normalize_path <- function(x) {
    normalizePath(x, winslash = "/", mustWork = FALSE)
  }

  md5 <- function(path) {
    unname(tools::md5sum(path))
  }

  sha256_file <- function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }

  read_csv_base <- function(path) {
    utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "N/A", "<NA>")
    )
  }

  write_csv_base <- function(x, path) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }

  latest_run_with <- function(base_dir, required_basename = NULL) {
    check(dir.exists(base_dir), paste("Missing directory:", base_dir))
    candidates <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    if (!is.null(required_basename)) {
      candidates <- candidates[
        file.exists(file.path(candidates, required_basename))
      ]
    }
    check(length(candidates) > 0L, paste("No completed run found in", base_dir))
    candidates[[which.max(file.info(candidates)$mtime)]]
  }

  newest_matching_file <- function(base_dir, pattern) {
    files <- list.files(
      base_dir,
      pattern = pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    check(length(files) > 0L, paste("No file matching", pattern, "in", base_dir))
    files[[which.max(file.info(files)$mtime)]]
  }

  copy_verified <- function(from, to) {
    check(file.exists(from), paste("Missing source asset:", from))
    check(file.info(from)$size > 0, paste("Empty source asset:", from))
    check(!file.exists(to), paste("Refusing to overwrite:", to))
    check(file.copy(from, to, overwrite = FALSE), paste("Copy failed:", from))
    check(
      identical(md5(from), md5(to)),
      paste("Copied asset failed MD5 verification:", basename(to))
    )
    invisible(to)
  }

  save_plot_pair <- function(plot, stem, width, height, dpi = 400) {
    png_path <- paste0(stem, ".png")
    pdf_path <- paste0(stem, ".pdf")
    check(!file.exists(png_path), paste("Refusing to overwrite:", png_path))
    check(!file.exists(pdf_path), paste("Refusing to overwrite:", pdf_path))
    ggplot2::ggsave(
      png_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    ggplot2::ggsave(
      pdf_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
    check(
      all(file.exists(c(png_path, pdf_path))),
      paste("Plot export failed for", basename(stem))
    )
    c(png = png_path, pdf = pdf_path)
  }

  base_theme <- ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#444444"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#777777", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )

  # ---------------------------------------------------------------------------
  # 1. Locate and validate the frozen 06B evidence set
  # ---------------------------------------------------------------------------
  evidence_base <- file.path("results", "revision", "evidence_claim_lock")
  evidence_dir <- latest_run_with(evidence_base, "revision_evidence_lock.txt")

  evidence_files <- c(
    "reportable_results.csv",
    "claim_guardrails.csv",
    "reviewer_comment_evidence_map.csv",
    "figure_table_placement_plan.csv",
    "manuscript_revision_checklist.csv",
    "input_inventory.csv",
    "output_checksums.csv",
    "reportable_text_blocks.txt",
    "revision_evidence_lock.txt"
  )
  evidence_paths <- file.path(evidence_dir, evidence_files)
  check(all(file.exists(evidence_paths)), "The latest 06B evidence set is incomplete.")

  output_checksums <- read_csv_base(
    file.path(evidence_dir, "output_checksums.csv")
  )
  checksum_name_col <- intersect(
    c("file", "filename", "output", "basename"),
    names(output_checksums)
  )
  checksum_md5_col <- intersect(
    c("md5", "MD5", "checksum", "output_md5"),
    names(output_checksums)
  )
  check(
    length(checksum_name_col) == 1L && length(checksum_md5_col) == 1L,
    "Could not identify filename and MD5 columns in output_checksums.csv."
  )

  checksum_targets <- file.path(
    evidence_dir,
    output_checksums[[checksum_name_col]]
  )
  check(all(file.exists(checksum_targets)), "A locked 06B output is missing.")
  checksum_ok <- md5(checksum_targets) ==
    as.character(output_checksums[[checksum_md5_col]])
  check(all(checksum_ok), "A locked 06B output checksum has changed.")

  reportable <- read_csv_base(file.path(evidence_dir, "reportable_results.csv"))
  check(
    identical(as.character(reportable$result_id), sprintf("R%02d", 1:18)),
    "The locked reportable-results index is not R01-R18 in order."
  )

  expected_values <- c(
    R01 = "22197", R02 = "3195", R03 = "898", R04 = "24703",
    R05 = "0.955", R06 = "567", R07 = "121", R08 = "-1.845",
    R09 = "-1.375", R10 = "0", R11 = "779", R12 = "100/100",
    R13 = "9", R14 = "6", R15 = "77.12", R16 = "0",
    R17 = "0", R18 = "20"
  )
  observed_values <- setNames(as.character(reportable$value), reportable$result_id)
  check(
    identical(observed_values[names(expected_values)], expected_values),
    "Locked reportable values differ from the approved 06B evidence set."
  )

  # The value printed by 06B is recorded here for downstream provenance. The
  # individual evidence-file MD5 checks above are the executable integrity test.
  evidence_lock_sha256 <- EXPECTED_EVIDENCE_LOCK_SHA256

  # ---------------------------------------------------------------------------
  # 2. Required completed result tables
  # ---------------------------------------------------------------------------
  required_tables <- c(
    marker_donors = file.path(
      "results", "revision", "marker_ifn_state_donor_counts.csv"
    ),
    marker_concordance = file.path(
      "results", "revision", "marker_ifn_state_cluster_concordance.csv"
    ),
    scoring_concordance = file.path(
      "results", "revision", "marker_ifn_state_scoring_concordance.csv"
    ),
    signature_variants = file.path(
      "results", "revision", "signature_variant_summary.csv"
    ),
    priority_queries = file.path(
      "results", "revision", "revision_priority_candidate_query_results.csv"
    ),
    integrated_candidates = file.path(
      "results", "revision", "revision_integrated_candidate_summary.csv"
    )
  )
  check(
    all(file.exists(required_tables)),
    paste(
      "Required completed result table(s) missing:",
      paste(names(required_tables)[!file.exists(required_tables)], collapse = ", ")
    )
  )

  marker_donors <- read_csv_base(required_tables[["marker_donors"]])
  marker_concordance <- read_csv_base(required_tables[["marker_concordance"]])
  scoring_concordance <- read_csv_base(required_tables[["scoring_concordance"]])
  signature_variants <- read_csv_base(required_tables[["signature_variants"]])
  priority_queries <- read_csv_base(required_tables[["priority_queries"]])
  integrated_candidates <- read_csv_base(required_tables[["integrated_candidates"]])

  check(
    all(c(
      "cohort", "patient", "total_cells", "marker_ifn_high_n",
      "marker_ifn_low_n", "marker_ifn_high_pct"
    ) %in% names(marker_donors)),
    "Unexpected marker_ifn_state_donor_counts.csv schema."
  )
  check(
    all(c(
      "marker_ifn_high_n", "cluster5_n", "marker_high_and_cluster5_n",
      "marker_high_outside_cluster5_n", "cluster5_below_marker_threshold_n"
    ) %in% names(marker_concordance)),
    "Unexpected marker_ifn_state_cluster_concordance.csv schema."
  )
  check(
    all(c(
      "query_id", "high_cells_n", "low_cells_n", "fdr_significant_genes_n"
    ) %in% names(signature_variants)),
    "Unexpected signature_variant_summary.csv schema."
  )
  check(
    all(c("pert_iname", "query_id", "NCS") %in% names(priority_queries)),
    "Unexpected revision_priority_candidate_query_results.csv schema."
  )

  # Locate the completed zero/tie comparison and curated identity table.
  zero_audit_dir <- latest_run_with(
    file.path("results", "revision", "zero_score_rank_tie_audit"),
    "all_compound_rank_tie_comparison.csv"
  )
  rank_tie_path <- file.path(
    zero_audit_dir, "all_compound_rank_tie_comparison.csv"
  )
  rank_tie <- read_csv_base(rank_tie_path)
  check(
    nrow(rank_tie) == 1376L &&
      all(c(
        "pert", "pert_iname", "primary_NCS", "average_tie_robustness_rank",
        "primary_FDR", "operational_hit_n", "loo_operational_hit_n",
        "empirical_lower_tail_p"
      ) %in% names(rank_tie)),
    "Unexpected all-compound average-tie comparison."
  )
  check(sum(rank_tie$primary_NCS == 0, na.rm = TRUE) == 779L,
        "PRIMARY_150 zero count is not the locked value 779.")

  curation_dir <- latest_run_with(
    file.path("results", "revision", "curated_compound_annotations"),
    "curation_decision_log.csv"
  )
  curation_path <- file.path(curation_dir, "curation_decision_log.csv")
  curation <- read_csv_base(curation_path)
  check(
    all(c(
      "pert", "curated_display_name", "curated_target", "curated_moa",
      "safe_for_qualified_reporting"
    ) %in% names(curation)),
    "Unexpected curated compound annotation table."
  )

  # ---------------------------------------------------------------------------
  # 3. Create a new immutable export run
  # ---------------------------------------------------------------------------
  run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_tag <- substr(digest::digest(
    paste(BUILD, evidence_lock_sha256, run_stamp, sep = "|"),
    algo = "sha256", serialize = FALSE
  ), 1L, 12L)
  run_id <- paste0("run_", run_stamp, "_", run_tag)

  figure_dir <- file.path(
    "figures", "revision", "locked_manuscript_assets", run_id
  )
  table_dir <- file.path(
    "results", "revision", "locked_manuscript_assets", run_id
  )
  check(!dir.exists(figure_dir), paste("Run directory already exists:", figure_dir))
  check(!dir.exists(table_dir), paste("Run directory already exists:", table_dir))
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # 4. Main Figure 1 component: marker-defined IFN-response-high state
  # ---------------------------------------------------------------------------
  canonical_donors <- marker_donors[
    marker_donors$cohort == "canonical_nine_patient",
    ,
    drop = FALSE
  ]
  canonical_donors$patient <- factor(
    paste0("P", canonical_donors$patient),
    levels = paste0("P", canonical_donors$patient)
  )

  state_counts <- data.frame(
    state = factor(
      c("IFN-response-low", "IFN-response-high"),
      levels = c("IFN-response-low", "IFN-response-high")
    ),
    cells = c(19002L, 3195L),
    stringsAsFactors = FALSE
  )

  panel_state_a <- ggplot2::ggplot(
    state_counts,
    ggplot2::aes(x = state, y = cells, fill = state)
  ) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(scales::comma(cells), "\n(",
                       scales::percent(cells / sum(cells), accuracy = 0.1), ")")
      ),
      vjust = -0.2,
      fontface = "bold",
      size = 3.4
    ) +
    ggplot2::scale_fill_manual(
      values = c("IFN-response-low" = "#91A6BA", "IFN-response-high" = "#D95F02")
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::labs(
      title = "A  Frozen marker-defined state",
      subtitle = "Canonical nine-donor PBMC population (n = 22,197)",
      x = NULL,
      y = "Candidate monocytes"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")

  panel_state_b <- ggplot2::ggplot(
    canonical_donors,
    ggplot2::aes(x = patient, y = marker_ifn_high_pct)
  ) +
    ggplot2::geom_col(width = 0.68, fill = "#D95F02") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", marker_ifn_high_pct)),
      vjust = -0.3,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.16))
    ) +
    ggplot2::labs(
      title = "B  Donor-level state frequency",
      subtitle = "Descriptive within-iNPH distribution; not a disease-specificity test",
      x = "Patient",
      y = "Marker-defined IFN-response-high cells"
    ) +
    base_theme

  cc <- marker_concordance[1L, , drop = FALSE]
  concordance_counts <- data.frame(
    category = factor(
      c(
        "Marker-high\noutside cluster 5",
        "Marker-high and\ncluster 5",
        "Cluster 5 below\nmarker threshold"
      ),
      levels = c(
        "Marker-high\noutside cluster 5",
        "Marker-high and\ncluster 5",
        "Cluster 5 below\nmarker threshold"
      )
    ),
    cells = c(
      cc$marker_high_outside_cluster5_n,
      cc$marker_high_and_cluster5_n,
      cc$cluster5_below_marker_threshold_n
    ),
    stringsAsFactors = FALSE
  )
  panel_state_c <- ggplot2::ggplot(
    concordance_counts,
    ggplot2::aes(x = category, y = cells, fill = category)
  ) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(cells)),
      vjust = -0.25,
      fontface = "bold",
      size = 3.4
    ) +
    ggplot2::scale_fill_manual(
      values = c("#7AA6C2", "#D95F02", "#BDBDBD")
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::labs(
      title = "C  Marker/cluster concordance",
      subtitle = "898/940 cluster-5 cells were marker-high (95.53%)",
      x = NULL,
      y = "Cells"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")

  main_state_plot <- panel_state_a / panel_state_b / panel_state_c +
    patchwork::plot_layout(heights = c(1, 1, 1)) +
    patchwork::plot_annotation(
      title = "Marker-defined IFN-response-high monocyte state observed in iNPH",
      subtitle = paste(
        "The marker rule is primary; cluster 5 is retained as a concordance result.",
        "No iNPH-versus-control disease signature is inferred."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 11, colour = "#444444")
      )
    )
  main_state_paths <- save_plot_pair(
    main_state_plot,
    file.path(figure_dir, "main_figure1_marker_state_summary"),
    width = 9.0,
    height = 12.0
  )

  # ---------------------------------------------------------------------------
  # 5. Main candidate panel: all 1,376 compounds with exact zeros retained
  # ---------------------------------------------------------------------------
  ranked <- rank_tie[order(rank_tie$primary_NCS, rank_tie$pert), , drop = FALSE]
  ranked$ordered_compound <- seq_len(nrow(ranked))
  ranked$score_class <- ifelse(
    ranked$primary_NCS < 0, "Negative",
    ifelse(ranked$primary_NCS > 0, "Positive", "Exact zero")
  )
  ranked$score_class <- factor(
    ranked$score_class,
    levels = c("Negative", "Exact zero", "Positive")
  )

  priority <- ranked[
    tolower(ranked$pert_iname) %in% c("rosuvastatin", "roflumilast"),
    ,
    drop = FALSE
  ]
  check(nrow(priority) == 2L, "Priority compounds not found in all-compound table.")

  panel_candidate_a <- ggplot2::ggplot(
    ranked,
    ggplot2::aes(x = ordered_compound, y = primary_NCS, colour = score_class)
  ) +
    ggplot2::geom_point(size = 0.9, alpha = 0.75) +
    ggplot2::geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.35) +
    ggplot2::geom_point(
      data = priority,
      size = 3.0,
      colour = c("#0072B2", "#D55E00")
    ) +
    ggrepel::geom_text_repel(
      data = priority,
      ggplot2::aes(label = pert_iname),
      seed = 42,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      show.legend = FALSE,
      size = 3.4
    ) +
    ggplot2::scale_colour_manual(
      values = c("Negative" = "#3B78B5", "Exact zero" = "#BDBDBD", "Positive" = "#C44E52")
    ) +
    ggplot2::labs(
      title = "A  PRIMARY_150 connectivity across all THP1 compounds",
      subtitle = "288 negative, 779 exact-zero and 309 positive NCS values",
      x = "Compounds ordered by PRIMARY_150 NCS",
      y = "Normalized connectivity score (NCS)",
      colour = "Score class"
    ) +
    base_theme

  top20 <- rank_tie[
    order(rank_tie$average_tie_robustness_rank, rank_tie$pert),
    ,
    drop = FALSE
  ]
  top20 <- top20[seq_len(min(20L, nrow(top20))), , drop = FALSE]
  top20 <- merge(
    top20,
    curation[, c("pert", "curated_display_name"), drop = FALSE],
    by = "pert",
    all.x = TRUE,
    sort = FALSE
  )
  top20$display_name <- ifelse(
    is.na(top20$curated_display_name) | !nzchar(top20$curated_display_name),
    top20$pert_iname,
    top20$curated_display_name
  )
  top20 <- top20[order(top20$average_tie_robustness_rank), , drop = FALSE]
  top20$display_name <- factor(
    top20$display_name,
    levels = rev(top20$display_name)
  )

  panel_candidate_b <- ggplot2::ggplot(
    top20,
    ggplot2::aes(x = primary_NCS, y = display_name)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = primary_NCS, yend = display_name),
      colour = "#B5B5B5",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = tolower(as.character(display_name)) %in%
          c("rosuvastatin", "roflumilast")
      ),
      size = 2.5
    ) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#4C78A8", `TRUE` = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "B  Integrated robustness top 20",
      subtitle = "Average-tie ordinal ranking; lower NCS indicates stronger predicted reversal",
      x = "PRIMARY_150 THP1 NCS",
      y = NULL
    ) +
    base_theme

  main_candidate_plot <- panel_candidate_a / panel_candidate_b +
    patchwork::plot_layout(heights = c(1, 1.35)) +
    patchwork::plot_annotation(
      title = "Zero-aware computational prioritization of transcriptional-reversal candidates",
      subtitle = paste(
        "Exact zero scores are retained as ties, not treated as missing.",
        "Candidate labels denote computational predictions, not therapeutic efficacy."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 11, colour = "#444444")
      )
    )
  main_candidate_paths <- save_plot_pair(
    main_candidate_plot,
    file.path(figure_dir, "main_candidate_prioritization_zero_aware"),
    width = 9.2,
    height = 11.5
  )

  # ---------------------------------------------------------------------------
  # 6. Supplementary Figure S6: state-definition sensitivity
  # ---------------------------------------------------------------------------
  state_query_ids <- c(
    "C5_ALLREST_MAX150",
    "C5_MARKERLOW_MAX150",
    "MARKER_MAX150",
    "RES03_ALLREST_MAX150",
    "RES07_ALLREST_MAX150",
    "RES10_ALLREST_MAX150",
    "P6_INCLUDED_MARKER_MAX150"
  )
  state_labels <- c(
    C5_ALLREST_MAX150 = "Cluster 5 vs all rest",
    C5_MARKERLOW_MAX150 = "Cluster 5 vs marker-low",
    MARKER_MAX150 = "Marker-high vs marker-low",
    RES03_ALLREST_MAX150 = "Resolution 0.3",
    RES07_ALLREST_MAX150 = "Resolution 0.7",
    RES10_ALLREST_MAX150 = "Resolution 1.0",
    P6_INCLUDED_MARKER_MAX150 = "Patient 6 included"
  )

  variant_subset <- signature_variants[
    signature_variants$query_id %in% state_query_ids,
    ,
    drop = FALSE
  ]
  check(
    setequal(variant_subset$query_id, state_query_ids),
    "State-definition signature summary is incomplete."
  )
  variant_subset$definition <- factor(
    unname(state_labels[variant_subset$query_id]),
    levels = unname(state_labels[state_query_ids])
  )

  panel_s6_a <- ggplot2::ggplot(
    variant_subset,
    ggplot2::aes(x = definition, y = high_cells_n)
  ) +
    ggplot2::geom_col(width = 0.68, fill = "#D95F02") +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(high_cells_n)),
      vjust = -0.25,
      size = 3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::labs(
      title = "A  Cells defining each sensitivity contrast",
      x = NULL,
      y = "IFN-response-high/selected cells"
    ) +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  priority_state <- priority_queries[
    priority_queries$query_id %in% state_query_ids &
      tolower(priority_queries$pert_iname) %in%
        c("rosuvastatin", "roflumilast"),
    ,
    drop = FALSE
  ]
  check(
    nrow(priority_state) == 2L * length(state_query_ids),
    "Priority-candidate state-definition results are incomplete."
  )
  priority_state$definition <- factor(
    unname(state_labels[priority_state$query_id]),
    levels = unname(state_labels[state_query_ids])
  )
  priority_state$pert_iname <- factor(
    tolower(priority_state$pert_iname),
    levels = c("rosuvastatin", "roflumilast")
  )

  panel_s6_b <- ggplot2::ggplot(
    priority_state,
    ggplot2::aes(
      x = definition,
      y = NCS,
      colour = pert_iname,
      group = pert_iname
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -0.5,
      linetype = "dashed",
      colour = "#777777",
      linewidth = 0.45
    ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(
      values = c("rosuvastatin" = "#0072B2", "roflumilast" = "#D55E00"),
      labels = c("Rosuvastatin", "Roflumilast")
    ) +
    ggplot2::labs(
      title = "B  Candidate connectivity across state definitions",
      subtitle = "Dashed line marks the operational NCS threshold of -0.5",
      x = NULL,
      y = "THP1 NCS",
      colour = NULL
    ) +
    base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  population_sensitivity <- data.frame(
    population = factor(
      c(
        "Canonical\n9 donors",
        "Uniform reconstruction\n10 donors",
        "Patient 6\nonly"
      ),
      levels = c(
        "Canonical\n9 donors",
        "Uniform reconstruction\n10 donors",
        "Patient 6\nonly"
      )
    ),
    total = c(22197L, 24703L, 2004L),
    marker_high = c(3195L, 3842L, 567L),
    stringsAsFactors = FALSE
  )
  population_sensitivity$marker_low <-
    population_sensitivity$total - population_sensitivity$marker_high
  population_long <- rbind(
    data.frame(
      population = population_sensitivity$population,
      state = "IFN-response-high",
      cells = population_sensitivity$marker_high
    ),
    data.frame(
      population = population_sensitivity$population,
      state = "IFN-response-low",
      cells = population_sensitivity$marker_low
    )
  )
  population_long$state <- factor(
    population_long$state,
    levels = c("IFN-response-low", "IFN-response-high")
  )

  panel_s6_c <- ggplot2::ggplot(
    population_long,
    ggplot2::aes(x = population, y = cells, fill = state)
  ) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::scale_fill_manual(
      values = c("IFN-response-low" = "#91A6BA", "IFN-response-high" = "#D95F02")
    ) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::annotate(
      "text",
      x = 2,
      y = 0,
      label = "Canonical-vs-uniform Jaccard = 0.955",
      vjust = 1.5,
      size = 3.2,
      fontface = "bold"
    ) +
    ggplot2::labs(
      title = "C  Patient 6 inclusion and population reconstruction",
      x = NULL,
      y = "Candidate monocytes",
      fill = NULL
    ) +
    base_theme

  supplementary_s6 <- panel_s6_a / panel_s6_b / panel_s6_c +
    patchwork::plot_layout(heights = c(1, 1, 1)) +
    patchwork::plot_annotation(
      title = "Sensitivity of cell-state definition and candidate connectivity",
      subtitle = paste(
        "Alternative clustering resolutions, comparator definitions, marker-defined",
        "states and Patient 6 inclusion are sensitivity analyses."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 11, colour = "#444444")
      )
    )
  supplementary_s6_paths <- save_plot_pair(
    supplementary_s6,
    file.path(figure_dir, "supplementary_figure_S6_state_definition_sensitivity"),
    width = 10.0,
    height = 13.0
  )

  # ---------------------------------------------------------------------------
  # 7. Standardize already-approved 04E and 05E components
  # ---------------------------------------------------------------------------
  average_tie_base <- file.path("figures", "revision", "average_ties")
  average_tie_dir <- latest_run_with(average_tie_base)

  approved_patterns <- c(
    candidate_heatmap = "candidate_robustness_heatmap_top20\\.png$",
    core_trajectories = "priority_candidate_core_rank_trajectories\\.png$",
    loo_ncs = "priority_candidate_loo_ncs_stability\\.png$",
    loo_rank = "priority_candidate_loo_rank_stability\\.png$",
    matched_null = "priority_candidate_matched_null_distributions\\.png$",
    average_tie_percentiles = "primary_vs_core_percentiles_average_ties\\.png$"
  )
  approved_sources <- vapply(
    approved_patterns,
    function(pattern) newest_matching_file(average_tie_dir, pattern),
    character(1)
  )
  approved_targets <- setNames(
    file.path(
      figure_dir,
      paste0("approved_04E_", names(approved_sources), ".png")
    ),
    names(approved_sources)
  )
  invisible(mapply(copy_verified, approved_sources, approved_targets))

  csf_source <- newest_matching_file(
    file.path("figures", "revision", "csf_contextualization"),
    "csf_contextualization_donor_transparency\\.png$"
  )
  csf_target <- file.path(
    figure_dir,
    "supplementary_figure_S7_csf_donor_transparency.png"
  )
  copy_verified(csf_source, csf_target)

  # ---------------------------------------------------------------------------
  # 8. Supplementary Table S5: public donor characteristics
  # ---------------------------------------------------------------------------
  donor_dir <- latest_run_with(
    file.path("results", "revision", "public_donor_characteristics")
  )
  donor_csvs <- list.files(
    donor_dir,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  check(length(donor_csvs) > 0L, "No 05D donor-characteristics CSV was found.")

  donor_candidates <- lapply(donor_csvs, function(path) {
    tryCatch(read_csv_base(path), error = function(e) NULL)
  })
  donor_scores <- vapply(donor_candidates, function(x) {
    if (is.null(x)) return(-Inf)
    name_score <- sum(grepl(
      "patient|donor|moca|age|sex|medication|pbmc|csf",
      names(x),
      ignore.case = TRUE
    ))
    row_score <- if (nrow(x) == 10L) 20 else 0
    name_score + row_score
  }, numeric(1))
  donor_index <- which.max(donor_scores)
  check(is.finite(donor_scores[[donor_index]]),
        "Could not identify the 05D donor-characteristics table.")
  donor_table <- donor_candidates[[donor_index]]
  check(nrow(donor_table) == 10L,
        "The selected 05D donor-characteristics table does not contain 10 donors.")

  table_s5_path <- file.path(
    table_dir,
    "supplementary_table_S5_public_donor_characteristics.csv"
  )
  write_csv_base(donor_table, table_s5_path)

  # ---------------------------------------------------------------------------
  # 9. Supplementary Table S6: compact top-20 robustness summary
  # ---------------------------------------------------------------------------
  table_s6 <- rank_tie[
    order(rank_tie$average_tie_robustness_rank, rank_tie$pert),
    c(
      "average_tie_robustness_rank", "pert", "pert_iname", "primary_NCS",
      "primary_FDR", "core_negative_n", "core_zero_n", "operational_hit_n",
      "loo_operational_hit_n", "empirical_lower_tail_p"
    ),
    drop = FALSE
  ]
  table_s6 <- table_s6[seq_len(20L), , drop = FALSE]
  table_s6 <- merge(
    table_s6,
    curation[, c(
      "pert", "curated_display_name", "curated_target", "curated_moa",
      "curation_status", "manuscript_handling"
    ), drop = FALSE],
    by = "pert",
    all.x = TRUE,
    sort = FALSE
  )
  table_s6 <- table_s6[
    order(table_s6$average_tie_robustness_rank),
    ,
    drop = FALSE
  ]
  table_s6$display_name <- ifelse(
    is.na(table_s6$curated_display_name) |
      !nzchar(table_s6$curated_display_name),
    table_s6$pert_iname,
    table_s6$curated_display_name
  )
  table_s6$core_negative <- paste0(table_s6$core_negative_n, "/11")
  table_s6$core_operational_hits <- paste0(table_s6$operational_hit_n, "/11")
  table_s6$loo_operational_hits <- paste0(table_s6$loo_operational_hit_n, "/9")
  table_s6 <- table_s6[, c(
    "average_tie_robustness_rank", "pert", "display_name", "primary_NCS",
    "primary_FDR", "core_negative", "core_zero_n", "core_operational_hits",
    "loo_operational_hits", "empirical_lower_tail_p", "curated_target",
    "curated_moa", "curation_status", "manuscript_handling"
  ), drop = FALSE]
  names(table_s6)[1:10] <- c(
    "integrated_average_tie_rank", "BRD_ID", "compound", "PRIMARY_150_THP1_NCS",
    "PRIMARY_150_WTCS_FDR", "negative_core_analyses", "zero_core_analyses",
    "operational_core_hits", "operational_LOO_hits",
    "nominal_matched_null_lower_tail_p"
  )
  table_s6_path <- file.path(
    table_dir,
    "supplementary_table_S6_compact_candidate_robustness.csv"
  )
  write_csv_base(table_s6, table_s6_path)

  # Export a compact state-sensitivity table used by S6 for traceability.
  table_s6_state <- merge(
    variant_subset[, c(
      "query_id", "description", "resolution", "patient6_included",
      "patients_n", "high_cells_n", "low_cells_n", "fdr_significant_genes_n",
      "signature_up_n", "signature_down_n"
    ), drop = FALSE],
    priority_state[, c("query_id", "pert_iname", "NCS"), drop = FALSE],
    by = "query_id",
    all.x = TRUE,
    sort = FALSE
  )
  table_s6_state_path <- file.path(
    table_dir,
    "supplementary_figure_S6_plot_data.csv"
  )
  write_csv_base(table_s6_state, table_s6_state_path)

  # ---------------------------------------------------------------------------
  # 10. Legends and manual PowerPoint assembly guide
  # ---------------------------------------------------------------------------
  legends <- c(
    "MAIN FIGURE 1 — MARKER-DEFINED STATE",
    paste(
      "The frozen marker rule classified 3,195 of 22,197 canonical PBMC",
      "classical-monocyte candidates from nine donors as IFN-response-high."
    ),
    paste(
      "Donor-level frequencies are descriptive within iNPH and do not test",
      "iNPH specificity. Of 940 cluster-5 cells, 898 (95.53%) were marker-high;",
      "cluster 5 is therefore presented as a concordance result rather than the",
      "primary state definition."
    ),
    "",
    "MAIN CANDIDATE FIGURE — ZERO-AWARE PRIORITIZATION",
    paste(
      "PRIMARY_150 THP1 results are displayed for all 1,376 compounds, including",
      "779 exact-zero NCS values. Average ties are used for percentile and",
      "integrated ordinal-rank reporting. Rosuvastatin (NCS -1.845; integrated",
      "rank 1) and roflumilast (NCS -1.375; integrated rank 15) are labelled as",
      "computationally prioritized transcriptional-reversal candidates."
    ),
    "",
    "SUPPLEMENTARY FIGURE S6 — STATE-DEFINITION SENSITIVITY",
    paste(
      "Candidate-cell counts and THP1 NCS values are shown across alternative",
      "cluster resolutions, comparator definitions, a marker-defined contrast",
      "and uniform inclusion of Patient 6. The submitted PRIMARY_150 query remains",
      "the frozen primary analysis; all alternatives are sensitivity analyses."
    ),
    "",
    "SUPPLEMENTARY FIGURE S7 — CSF CONTEXTUALIZATION",
    paste(
      "The blood-derived program is shown only as supportive cross-compartment",
      "detectability. Unequal CSF candidate-cell recovery precludes pooled cohort",
      "prevalence, definitive enrichment, CNS-origin, causal or efficacy claims."
    )
  )
  legend_path <- file.path(table_dir, "locked_figure_legends.txt")
  writeLines(legends, legend_path, useBytes = TRUE)

  assembly_guide <- c(
    "LOCKED MANUSCRIPT ASSET ASSEMBLY GUIDE",
    paste("Build:", BUILD),
    paste("Evidence lock SHA-256:", evidence_lock_sha256),
    "",
    "MAIN MANUSCRIPT",
    "Figure 1: use main_figure1_marker_state_summary as the marker-state summary.",
    paste(
      "Figure 2: retain the external SLE/healthy-control contextualization only",
      "with explicitly descriptive, non-specificity wording."
    ),
    "Figure 3: use main_candidate_prioritization_zero_aware.",
    paste(
      "Figure 4: retain the existing cell-line context figure; state that THP1 is",
      "a leukemia-derived monocytic screening model and that dose/duration could",
      "not be linked from the executed saved results."
    ),
    paste(
      "Figure 5: assemble approved_04E_core_trajectories,",
      "approved_04E_loo_ncs and approved_04E_matched_null as panels A-C."
    ),
    "",
    "SUPPLEMENT",
    "Retain the submitted Supplementary Figures S1-S5 and Tables S1-S4.",
    "Use supplementary_figure_S6_state_definition_sensitivity as Figure S6.",
    "Use supplementary_figure_S7_csf_donor_transparency as Figure S7.",
    "Use supplementary_table_S5_public_donor_characteristics as Table S5.",
    "Use supplementary_table_S6_compact_candidate_robustness as Table S6.",
    "",
    "REVIEWER-SPECIFIC EVIDENCE FILES ONLY — DO NOT ADD TO THE SUPPLEMENT",
    paste(
      "approved_04E_candidate_heatmap, approved_04E_loo_rank and",
      "approved_04E_average_tie_percentiles are dense diagnostics suitable for",
      "Reviewer 1/Reviewer 2 evidence files."
    ),
    "",
    "SUPERSEDED FIGURES — DO NOT USE",
    "primary_rank_vs_objective_robustness_rank.png",
    "primary_vs_core_percentiles_v2.png",
    "any pre-04E minimum-tie version",
    "the initial 05E CSF figure without readable P1-P10 labels"
  )
  assembly_path <- file.path(table_dir, "powerpoint_assembly_guide.txt")
  writeLines(assembly_guide, assembly_path, useBytes = TRUE)

  # ---------------------------------------------------------------------------
  # 11. Asset inventory and final validation
  # ---------------------------------------------------------------------------
  exported_files <- c(
    unname(main_state_paths),
    unname(main_candidate_paths),
    unname(supplementary_s6_paths),
    approved_targets,
    csf_target,
    table_s5_path,
    table_s6_path,
    table_s6_state_path,
    legend_path,
    assembly_path
  )
  check(all(file.exists(exported_files)), "One or more exported assets are missing.")
  check(all(file.info(exported_files)$size > 0), "One or more exported assets are empty.")

  role <- ifelse(
    grepl("main_", basename(exported_files)), "MAIN_MANUSCRIPT",
    ifelse(
      grepl("supplementary_", basename(exported_files)), "SUPPLEMENT",
      ifelse(
        grepl("approved_04E_", basename(exported_files)),
        "POWERPOINT_COMPONENT_OR_REVIEWER_EVIDENCE",
        "DOCUMENTATION"
      )
    )
  )
  inventory <- data.frame(
    build = BUILD,
    evidence_lock_sha256 = evidence_lock_sha256,
    role = role,
    file = normalize_path(exported_files),
    size_bytes = file.info(exported_files)$size,
    md5 = md5(exported_files),
    sha256 = vapply(exported_files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  inventory_path <- file.path(table_dir, "revision_manuscript_asset_inventory.csv")
  write_csv_base(inventory, inventory_path)

  result <- list(
    analysis_id = "revision_06C_locked_manuscript_asset_export",
    build = BUILD,
    evidence_lock_sha256 = evidence_lock_sha256,
    evidence_dir = normalize_path(evidence_dir),
    figure_dir = normalize_path(figure_dir),
    table_dir = normalize_path(table_dir),
    manuscript_main = c(
      marker_state_summary = normalize_path(main_state_paths[["png"]]),
      zero_aware_candidates = normalize_path(main_candidate_paths[["png"]]),
      robustness_components = normalize_path(
        approved_targets[c("core_trajectories", "loo_ncs", "matched_null")]
      )
    ),
    supplement = c(
      figure_S6 = normalize_path(supplementary_s6_paths[["png"]]),
      figure_S7 = normalize_path(csf_target),
      table_S5 = normalize_path(table_s5_path),
      table_S6 = normalize_path(table_s6_path)
    ),
    reviewer_evidence_only = normalize_path(
      approved_targets[c(
        "candidate_heatmap", "loo_rank", "average_tie_percentiles"
      )]
    ),
    asset_inventory = normalize_path(inventory_path),
    assembly_guide = normalize_path(assembly_path),
    new_biological_analyses = 0L,
    new_lincs_queries = 0L
  )
  result_path <- file.path(table_dir, "revision_locked_manuscript_assets.rds")
  saveRDS(result, result_path)
  result_check <- readRDS(result_path)
  check(identical(result_check$build, BUILD), "Saved result verification failed.")
  check(
    identical(result_check$evidence_lock_sha256, evidence_lock_sha256),
    "Saved evidence-lock identifier verification failed."
  )

  cat("\n============================================================\n")
  cat("  REVISION ANALYSIS 06C COMPLETE\n")
  cat("============================================================\n")
  cat("  Status                 : PASS_LOCKED_MANUSCRIPT_ASSETS_EXPORTED\n")
  cat("  Evidence lock SHA-256 :", evidence_lock_sha256, "\n")
  cat("  Main figure designs   : 2 newly generated + approved components\n")
  cat("  New Supplement figure: S6\n")
  cat("  Reused Supplement fig : S7 (final readable-label version)\n")
  cat("  New Supplement tables : S5 and S6\n")
  cat("  New biological analyses:", result$new_biological_analyses, "\n")
  cat("  New LINCS queries       :", result$new_lincs_queries, "\n")
  cat("  Figure directory        :", normalize_path(figure_dir), "\n")
  cat("  Table/guide directory   :", normalize_path(table_dir), "\n")
  cat("============================================================\n")
  cat("Open powerpoint_assembly_guide.txt before assembling panels.\n")
  cat("Inspect every exported PNG and compare panel text with the locked legends.\n")
})
