# =============================================================================
# Revision analysis 05E — Transparent CSF contextualization figure
# Install as analysis/05E_visualize_csf_contextualization.R
# Build: 05E_CSF_CONTEXTUALIZATION_FIGURE_V1_20260907
# =============================================================================
#
# Creates a supplementary diagnostic figure from the completed 05B and 05C
# results. No cells are reclustered, no scores are recalculated, and no LINCS
# query is run. The figure explicitly displays donor-level CSF cell yields and
# uncertainty so the contextualization cannot be mistaken for a CSF prevalence
# or enrichment analysis.
# =============================================================================

local({
  BUILD <- "05E_CSF_CONTEXTUALIZATION_FIGURE_V1_20260907"

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  required_packages <- c("ggplot2", "patchwork")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  check(
    length(missing_packages) == 0L,
    paste("Missing package(s):", paste(missing_packages, collapse = ", "))
  )

  csf_path <- "data/processed/revision_csf_ifn_contextualization.rds"
  check(file.exists(csf_path), paste("Missing input:", csf_path))

  audit_parent <- "results/revision/csf_cell_count_heterogeneity_audit"
  audit_paths <- if (dir.exists(audit_parent)) {
    list.files(
      audit_parent,
      pattern = "^revision_csf_cell_count_heterogeneity_audit[.]rds$",
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    character()
  }
  audit_paths <- audit_paths[
    file.exists(file.path(dirname(audit_paths), "_SUCCESS.txt"))
  ]
  check(length(audit_paths) > 0L, "No completed 05C audit could be found.")
  audit_path <- audit_paths[which.max(file.info(audit_paths)$mtime)]

  input_paths <- c(csf_path, audit_path)
  input_md5_before <- unname(tools::md5sum(input_paths))
  check(!anyNA(input_md5_before), "An input checksum could not be calculated.")

  csf <- readRDS(csf_path)
  audit <- readRDS(audit_path)
  check(
    identical(csf$analysis_id, "revision_05B_csf_monocyte_ifn_contextualization"),
    "Unexpected 05B result."
  )
  check(
    identical(audit$analysis_id, "revision_05C_csf_cell_count_heterogeneity_audit"),
    "Unexpected 05C audit."
  )
  check(
    identical(audit$source_05B_md5, input_md5_before[[1]]),
    "The completed 05C audit does not correspond to the current 05B result."
  )

  paired <- as.data.frame(csf$paired_donor_summary)
  donor <- as.data.frame(audit$donor_audit)
  sensitivity <- as.data.frame(audit$minimum_cell_sensitivity)

  required_paired <- c(
    "patient", "n_cells_CSF", "median_transfer_score_PBMC",
    "median_transfer_score_CSF"
  )
  required_donor <- c(
    "patient", "n_cells", "marker_ifn_high_pct",
    "marker_ifn_high_wilson_lower_pct",
    "marker_ifn_high_wilson_upper_pct"
  )
  check(
    all(required_paired %in% names(paired)) &&
      all(required_donor %in% names(donor)),
    "Unexpected 05B/05C donor schema."
  )

  paired$patient <- as.character(paired$patient)
  donor$patient <- as.character(donor$patient)
  paired <- paired[order(as.numeric(paired$patient)), , drop = FALSE]
  donor <- donor[order(as.numeric(donor$patient)), , drop = FALSE]
  paired$yield_group <- ifelse(
    paired$n_cells_CSF >= 30L,
    ">=30 CSF candidates",
    "<30 CSF candidates"
  )
  donor$yield_group <- ifelse(
    donor$n_cells >= 30L,
    ">=30 CSF candidates",
    "<30 CSF candidates"
  )

  paired_long <- rbind(
    data.frame(
      patient = paired$patient,
      compartment = "PBMC",
      compartment_order = 1,
      median_transfer_score = paired$median_transfer_score_PBMC,
      n_CSF_candidates = paired$n_cells_CSF,
      yield_group = paired$yield_group,
      stringsAsFactors = FALSE
    ),
    data.frame(
      patient = paired$patient,
      compartment = "Ventricular CSF",
      compartment_order = 2,
      median_transfer_score = paired$median_transfer_score_CSF,
      n_CSF_candidates = paired$n_cells_CSF,
      yield_group = paired$yield_group,
      stringsAsFactors = FALSE
    )
  )
  paired_long$compartment <- factor(
    paired_long$compartment,
    levels = c("PBMC", "Ventricular CSF")
  )
  paired_long$yield_group <- factor(
    paired_long$yield_group,
    levels = c(">=30 CSF candidates", "<30 CSF candidates")
  )
  donor$yield_group <- factor(
    donor$yield_group,
    levels = c(">=30 CSF candidates", "<30 CSF candidates")
  )

  transfer_sensitivity <- sensitivity[
    sensitivity$metric == "median blood-derived transfer score",
    ,
    drop = FALSE
  ]
  check(
    all(c(1L, 30L) %in%
      transfer_sensitivity$minimum_CSF_candidates_per_donor),
    "05C transfer-score sensitivity lacks the 1- or 30-cell threshold."
  )
  row_all <- transfer_sensitivity[
    transfer_sensitivity$minimum_CSF_candidates_per_donor == 1L,
    ,
    drop = FALSE
  ]
  row_30 <- transfer_sensitivity[
    transfer_sensitivity$minimum_CSF_candidates_per_donor == 30L,
    ,
    drop = FALSE
  ]

  format_p <- function(x) {
    if (x < 0.001) format(x, scientific = TRUE, digits = 2) else sprintf("%.3f", x)
  }

  palette <- c(
    ">=30 CSF candidates" = "#0072B2",
    "<30 CSF candidates" = "#D55E00"
  )

  base_theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#444444"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold")
    )

  panel_a <- ggplot2::ggplot(
    paired_long,
    ggplot2::aes(
      x = compartment,
      y = median_transfer_score,
      group = patient,
      colour = yield_group
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.35) +
    ggplot2::geom_line(alpha = 0.75, linewidth = 0.65) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(
      data = paired_long[paired_long$compartment == "Ventricular CSF", ],
      ggplot2::aes(label = patient),
      hjust = -0.65,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.12, 0.24))) +
    ggplot2::labs(
      title = "Paired transfer of the blood-derived program",
      subtitle = paste0(
        "All paired donors: n=", row_all$eligible_donors_n,
        ", nominal P=", format_p(row_all$nominal_unadjusted_p),
        "; >=30-cell sensitivity: n=", row_30$eligible_donors_n,
        ", nominal P=", format_p(row_30$nominal_unadjusted_p)
      ),
      x = NULL,
      y = "Median transfer score"
    ) +
    base_theme

  panel_b <- ggplot2::ggplot(
    donor,
    ggplot2::aes(
      x = factor(patient, levels = patient),
      y = n_cells,
      fill = yield_group
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 30,
      linetype = "dashed",
      colour = "#555555",
      linewidth = 0.5
    ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_cells),
      vjust = -0.35,
      size = 3
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_log10(
      breaks = c(1, 10, 30, 100, 1000, 3000),
      labels = c("1", "10", "30", "100", "1,000", "3,000"),
      expand = ggplot2::expansion(mult = c(0.02, 0.14))
    ) +
    ggplot2::labs(
      title = "Unequal CSF candidate-cell recovery",
      subtitle = "Dashed line marks 30 candidates; logarithmic y-axis",
      x = "Patient",
      y = "CSF candidates"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")

  panel_c <- ggplot2::ggplot(
    donor,
    ggplot2::aes(
      x = factor(patient, levels = patient),
      y = marker_ifn_high_pct,
      colour = yield_group
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = marker_ifn_high_wilson_lower_pct,
        ymax = marker_ifn_high_wilson_upper_pct
      ),
      width = 0.15,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::labs(
      title = "Descriptive IFN-high fractions and uncertainty",
      subtitle = "Points are donor estimates; bars are 95% Wilson intervals",
      x = "Patient",
      y = "CSF IFN-high cells"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")

  combined <- panel_a / (panel_b + panel_c) +
    patchwork::plot_annotation(
      title = "Ventricular-CSF contextualization of the blood-derived IFN-response-high program",
      subtitle = paste(
        "Supportive cross-compartment analysis only; unequal cell recovery",
        "precludes cohort prevalence or definitive enrichment inference"
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 11, colour = "#444444")
      )
    ) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))

  run_id <- paste0(
    "run_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    substr(input_md5_before[[2]], 1, 12)
  )
  figure_directory <- file.path(
    "figures/revision/csf_contextualization",
    run_id
  )
  result_directory <- file.path(
    "results/revision/csf_contextualization_figures",
    run_id
  )
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

  png_path <- file.path(
    figure_directory,
    "csf_contextualization_donor_transparency.png"
  )
  pdf_path <- file.path(
    figure_directory,
    "csf_contextualization_donor_transparency.pdf"
  )
  ggplot2::ggsave(
    png_path,
    combined,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_path,
    combined,
    width = 12,
    height = 10,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  utils::write.csv(
    paired_long,
    file.path(result_directory, "csf_paired_transfer_plot_data.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    donor,
    file.path(result_directory, "csf_donor_yield_uncertainty_plot_data.csv"),
    row.names = FALSE
  )

  figure_legend <- paste(
    "Exploratory ventricular-CSF contextualization of the frozen",
    "blood-derived IFN-response-high monocyte program.",
    "(A) Paired donor-level median transfer scores in canonical PBMC",
    "classical monocytes and independently selected ventricular-CSF",
    "classical-monocyte candidates. Colors indicate whether each donor",
    "contributed at least 30 CSF candidates. P values are nominal,",
    "two-sided paired Wilcoxon tests. (B) CSF candidate-cell recovery per",
    "donor on a logarithmic scale. The dashed line marks 30 cells.",
    "(C) Donor-level percentages of CSF candidates above the frozen",
    "blood-derived IFN threshold with 95% Wilson intervals. Because CSF",
    "cell yields were markedly heterogeneous, these data support",
    "cross-compartment detectability only and do not establish cohort",
    "prevalence, CSF enrichment, CNS origin, disease specificity, cognition",
    "association, or therapeutic relevance."
  )
  writeLines(
    figure_legend,
    file.path(result_directory, "csf_contextualization_figure_legend.txt")
  )

  result <- list(
    analysis_id = "revision_05E_visualize_csf_contextualization",
    build_id = BUILD,
    status = "PASS_SUPPORTIVE_SUPPLEMENTARY_FIGURE_CREATED",
    source_05B = csf_path,
    source_05C = audit_path,
    input_md5 = data.frame(
      input = input_paths,
      md5 = input_md5_before,
      stringsAsFactors = FALSE
    ),
    paired_plot_data = paired_long,
    donor_plot_data = donor,
    transfer_sensitivity = transfer_sensitivity,
    figure_legend = figure_legend,
    figure_paths = c(PNG = png_path, PDF = pdf_path),
    result_directory = result_directory,
    session_info = utils::sessionInfo()
  )
  result_path <- file.path(
    result_directory,
    "revision_csf_contextualization_figures.rds"
  )
  saveRDS(result, result_path)
  check(identical(readRDS(result_path), result), "05E RDS read-back failed.")

  check(
    identical(input_md5_before, unname(tools::md5sum(input_paths))),
    "An input changed during Analysis 05E."
  )
  outputs <- c(
    list.files(figure_directory, full.names = TRUE),
    list.files(result_directory, full.names = TRUE)
  )
  check(all(file.info(outputs)$size > 0), "An Analysis 05E output is empty.")
  utils::write.csv(
    data.frame(
      file = outputs,
      md5 = unname(tools::md5sum(outputs)),
      stringsAsFactors = FALSE
    ),
    file.path(result_directory, "output_checksums.csv"),
    row.names = FALSE
  )
  writeLines(
    c(result$status, BUILD),
    file.path(result_directory, "_SUCCESS.txt")
  )

  cat("\n============================================================\n")
  cat("  REVISION ANALYSIS 05E COMPLETE\n")
  cat("============================================================\n")
  cat("  Status                      :", result$status, "\n")
  cat("  Paired PBMC/CSF donors      :", nrow(paired), "\n")
  cat("  Donors with >=30 CSF cells :", sum(paired$n_cells_CSF >= 30L), "\n")
  cat("  Figure designs             : 1\n")
  cat("  PNG/PDF files              : 2\n")
  cat("  PNG                         :", png_path, "\n")
  cat("  Figure legend               :",
      file.path(result_directory, "csf_contextualization_figure_legend.txt"),
      "\n")
  cat("  New scoring or LINCS query : 0\n")
  cat("============================================================\n")
})
