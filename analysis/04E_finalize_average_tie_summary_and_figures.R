# =============================================================================
# Build: 04E_AVERAGE_TIES_V1_20260907
# Install as analysis/04E_finalize_average_tie_summary_and_figures.R
# =============================================================================
# Authorized downstream change: average ranks for all exact NCS ties, including
# zeros. Original submitted nonzero ranks, scores, FDRs, and null probabilities
# remain unchanged. Minimum-tie values are retained in explicit *_minimum fields.
# All outputs go to new run folders. No existing 04B/04C output is overwritten.
# No new LINCS query. Requires a completed 04D audit matching current input files.

local({
  SCRIPT_BUILD <- "04E_AVERAGE_TIES_V1_20260907"
  EXPECTED_LOCK <- paste0("1e3a2b827a1429e406303871d547a9dec",
                          "07ee4068078c666df69ca82e5428380")
  check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  same <- function(x, y) identical(as.numeric(x), as.numeric(y))
  key <- function(x) paste(x$query_id, x$pert, sep = "::")
  update_columns <- function(old, revised, cols) {
    check(!anyDuplicated(revised$pert) && setequal(old$pert, revised$pert),
          "Compound universe differs during summary update.")
    z <- match(old$pert, revised$pert)
    for (nm in cols) {
      old[[paste0(nm, "_minimum")]] <- old[[nm]]
      old[[nm]] <- revised[[nm]][z]
    }
    old
  }
  build_summary <- function(original, audit) {
    observed <- as.data.frame(original$observed_thp1)
    ar <- as.data.frame(audit$all_query_rank_comparison)
    check(nrow(observed) == 21L * 1376L && !anyDuplicated(key(observed)) &&
            !anyDuplicated(key(ar)), "Observed/audit query keys are not unique.")
    ix <- match(key(observed), key(ar))
    check(!anyNA(ix), "04D lacks an observed query/compound pair.")
    check(same(observed$NCS, ar$NCS[ix]) && same(observed$WTCS, ar$WTCS[ix]) &&
            same(observed$rank_all, ar$rank_min[ix]) &&
            same(observed$percentile_all, ar$percentile_min[ix]),
          "04B scores or minimum ranks differ from 04D.")
    check(all(is.finite(observed$NCS)), "Observed scores must all be finite.")
    # Independently reproduce each average rank before adopting the audit values.
    for (q in unique(observed$query_id)) {
      j <- which(observed$query_id == q)
      check(length(j) == 1376L, paste("Wrong THP1 population for", q))
      r <- as.numeric(rank(observed$NCS[j], ties.method = "average"))
      check(same(r, ar$rank_average[ix[j]]) &&
              same((r - 1) / 1375, ar$percentile_average[ix[j]]),
            paste("Average rank reproduction failed:", q))
    }
    observed$rank_all_minimum <- observed$rank_all
    observed$percentile_all_minimum <- observed$percentile_all
    observed$rank_all <- ar$rank_average[ix]
    observed$percentile_all <- ar$percentile_average[ix]
    observed$connectivity_direction <- ifelse(observed$NCS < 0, "negative",
                                             ifelse(observed$NCS > 0, "positive", "zero"))
    ids <- original$core_query_ids
    check(length(ids) == 11L && identical(ids, audit$core_query_ids), "Core query set mismatch.")
    perts <- sort(as.character(original$candidate_summary$pert), method = "radix")
    check(length(perts) == 1376L && !anyDuplicated(perts), "Expected 1,376 unique compounds.")
    cmp <- as.data.frame(audit$all_compound_comparison)
    check(nrow(cmp) == 1376L && !anyDuplicated(cmp$pert) && setequal(cmp$pert, perts),
          "Invalid all-compound audit comparison.")
    cmp <- cmp[match(perts, cmp$pert), ]
    primary <- observed[observed$query_id == "PRIMARY_150", ]
    primary <- primary[match(perts, primary$pert), ]
    core <- observed[observed$query_id %in% ids, ]
    loo <- observed[observed$query_role == "donor_sensitivity", ]
    core_new <- do.call(rbind, lapply(perts, function(p) {
      d <- core[core$pert == p, ]
      check(nrow(d) == 11L, "Incomplete core query coverage.")
      data.frame(pert = p, median_rank_all = median(d$rank_all),
                 min_rank_all = min(d$rank_all), max_rank_all = max(d$rank_all),
                 median_percentile_all = median(d$percentile_all),
                 min_percentile_all = min(d$percentile_all),
                 max_percentile_all = max(d$percentile_all), inverse_proportion = mean(d$NCS < 0))
    }))
    core_new$primary_rank_all <- primary$rank_all
    core_new$primary_percentile_all <- primary$percentile_all
    ordering <- order(core_new$median_percentile_all, -core_new$inverse_proportion,
                      core_new$primary_percentile_all, perts, method = "radix")
    core_new$objective_robustness_rank <- integer(length(perts))
    core_new$objective_robustness_rank[ordering] <- seq_along(perts)
    check(same(core_new$objective_robustness_rank, cmp$average_tie_robustness_rank) &&
            same(core_new$median_percentile_all, cmp$median_core_percentile_average),
          "Integrated average-tie ranking does not reproduce the completed audit.")
    old_cand <- as.data.frame(original$candidate_summary)
    check(same(old_cand$objective_robustness_rank,
               cmp$original_robustness_rank[match(old_cand$pert, cmp$pert)]),
          "Original integrated ranking differs from 04D.")
    rank_cols <- setdiff(names(core_new), c("pert", "inverse_proportion"))
    core_summary <- update_columns(as.data.frame(original$core_robustness), core_new, rank_cols)
    candidate_summary <- update_columns(old_cand, core_new, rank_cols)
    loo_new <- do.call(rbind, lapply(perts, function(p) {
      d <- loo[loo$pert == p, ]
      check(nrow(d) == 9L, "Incomplete LOO query coverage.")
      data.frame(pert = p, loo_median_rank_all = median(d$rank_all),
                 loo_min_rank_all = min(d$rank_all), loo_max_rank_all = max(d$rank_all))
    }))
    loo_cols <- setdiff(names(loo_new), "pert")
    loo_summary <- update_columns(as.data.frame(original$loo_summary), loo_new, loo_cols)
    candidate_summary <- update_columns(candidate_summary, loo_new, loo_cols)
    untouched <- setdiff(names(old_cand), c(rank_cols, loo_cols))
    check(identical(old_cand[untouched], candidate_summary[untouched]),
          "A non-ranking candidate field changed.")
    check(identical(as.data.frame(original$observed_thp1)[setdiff(names(original$observed_thp1),
                     c("rank_all", "percentile_all"))],
                   observed[setdiff(names(original$observed_thp1), c("rank_all", "percentile_all"))]),
          "A non-ranking observed field changed.")
    result <- original
    result$analysis_id <- "revision_candidate_robustness_average_ties_v1"
    result$script_build <- SCRIPT_BUILD
    result$ranking_method <- "average"
    result$observed_thp1 <- observed
    result$core_robustness <- core_summary[order(core_summary$objective_robustness_rank), ]
    result$loo_summary <- loo_summary
    result$candidate_summary <- candidate_summary[order(candidate_summary$objective_robustness_rank), ]
    focus <- unique(c(head(result$candidate_summary$pert, 30),
                       "BRD-K82941592", "BRD-K03194791"))
    result$priority_candidate_summary <- result$candidate_summary[result$candidate_summary$pert %in% focus, ]
    result$minimum_tie_comparison <- as.data.frame(audit$all_compound_comparison)
    result$ranking_definitions$percentile_all <-
      "Average rank of ascending NCS among 1,376 THP1 compounds; percentile=(rank-1)/1375; all exact ties included."
    result$ranking_definitions$objective_ranking <- paste(
      "Median average-tie percentile across the same 11 core queries; decreasing negative-NCS proportion,",
      "primary average-tie percentile, then BRD ID break remaining ties."
    )
    result$ranking_definitions$minimum_sensitivity <-
      "Original 04B minimum-tie rank fields are retained with the _minimum suffix. Submitted nonzero ranks are unchanged."
    result$baseline_input_provenance <- original$input_provenance
    result$output_paths <- NULL
    result
  }
  run_main <- function() {
    started <- Sys.time()
    for (p in c("dplyr", "ggplot2", "readr", "scales"))
      check(requireNamespace(p, quietly = TRUE), paste("Required package missing:", p))
    ROBUSTNESS_RDS <- "data/processed/revision_candidate_robustness_and_null_summary.rds"
    BATCH_INDEX_RDS <- "data/processed/revision_lincs_batch_index.rds"
    check(all(file.exists(c(ROBUSTNESS_RDS, BATCH_INDEX_RDS))), "Open the project root; required inputs are missing.")
    original <- readRDS(ROBUSTNESS_RDS)
    batch_index <- readRDS(BATCH_INDEX_RDS)
    check(identical(original$manifest_content_sha256, EXPECTED_LOCK) &&
            identical(batch_index$manifest_content_sha256, EXPECTED_LOCK) &&
            isTRUE(batch_index$all_locked_queries_complete), "Signature lock or batch completion mismatch.")
    chosen <- Sys.getenv("INPH_RANK_TIE_AUDIT_RDS", unset = "")
    audit_paths <- if (nzchar(chosen)) chosen else list.files(
      "results/revision/zero_score_rank_tie_audit", pattern = "^revision_zero_score_rank_tie_audit\\.rds$",
      recursive = TRUE, full.names = TRUE)
    audit_paths <- audit_paths[file.exists(audit_paths) &
                                 file.exists(file.path(dirname(audit_paths), "_SUCCESS.txt"))]
    check(length(audit_paths) > 0L, "No completed 04D audit found.")
    audit_paths <- audit_paths[order(file.info(audit_paths)$mtime, decreasing = TRUE)]
    audit_path <- NULL
    for (p in audit_paths) {
      a <- readRDS(p)
      if (identical(a$status, "PASS_AUDIT_COMPLETE_NOT_A_RANKING_REPLACEMENT") &&
          identical(a$manifest_content_sha256, EXPECTED_LOCK)) {
        audit <- a; audit_path <- p; break
      }
    }
    check(!is.null(audit_path), "No passing audit matches the locked signature manifest.")
    h <- as.data.frame(audit$input_hashes)
    check(all(file.exists(h$path)) && identical(unname(tools::md5sum(h$path)), as.character(h$md5_after)),
          "An input recorded by 04D changed. Investigate before finalizing ranks.")
    required_files <- unique(c(h$path, audit_path, ROBUSTNESS_RDS, BATCH_INDEX_RDS))
    own <- "analysis/04E_finalize_average_tie_summary_and_figures.R"
    if (file.exists(own)) required_files <- unique(c(required_files, own))
    hashes_before <- unname(tools::md5sum(required_files))
    cat("Using completed audit:", audit_path, "\n")
    robustness <- build_summary(original, audit)
    candidate_summary <- robustness$candidate_summary
    observed_thp1 <- robustness$observed_thp1
    core_query_ids <- robustness$core_query_ids
    parent <- "results/revision/average_tie_finalization"
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    RESULT_DIR <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
    check(dir.create(RESULT_DIR), "Cannot create a fresh finalization run directory.")
    FIGURE_DIR <- file.path("figures/revision/average_ties", basename(RESULT_DIR))
    check(dir.create(FIGURE_DIR, recursive = TRUE), "Cannot create the figure directory.")
    complete <- FALSE
    on.exit(if (!complete) message("Finalization incomplete; retain original outputs. Partial run: ", RESULT_DIR), add = TRUE)
    priority_perts <- candidate_summary |>
      dplyr::filter(tolower(pert_iname) %in% c("rosuvastatin", "roflumilast")) |>
      dplyr::pull(pert)

    if (length(priority_perts) != 2L) {
      stop("Rosuvastatin and roflumilast were not uniquely identified.")
    }

    display_name <- function(pert_iname, pert) {
      ifelse(
        is.na(pert_iname) | pert_iname == "" | pert_iname == "UNIDENTIFIED",
        pert,
        pert_iname
      )
    }

    theme_revision <- function(base_size = 11) {
      ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
          plot.title = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(color = "grey30"),
          legend.position = "right"
        )
    }

    save_plot_pair <- function(plot, stem, width, height) {
      png_path <- file.path(FIGURE_DIR, paste0(stem, ".png"))
      pdf_path <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))

      ggplot2::ggsave(
        png_path,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = 400,
        bg = "white"
      )
      ggplot2::ggsave(
        pdf_path,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
        bg = "white"
      )

      c(png = png_path, pdf = pdf_path)
    }

    query_labels <- c(
      PRIMARY_150 = "Primary\n150/150",
      PRIMARY_50 = "Primary\n50/50",
      PRIMARY_100 = "Primary\n100/100",
      PRIMARY_200 = "Primary\n200/200",
      C5_ALLREST_MAX150 = "C5 vs\nall rest",
      C5_MARKERLOW_MAX150 = "C5 vs\nmarker-low",
      MARKER_MAX150 = "Marker-high\nvs low",
      RES03_ALLREST_MAX150 = "Resolution\n0.3",
      RES07_ALLREST_MAX150 = "Resolution\n0.7",
      RES10_ALLREST_MAX150 = "Resolution\n1.0",
      P6_INCLUDED_MARKER_MAX150 = "Patient 6\nincluded"
    )

    if (!setequal(names(query_labels), core_query_ids)) {
      stop("Core-query labels do not match the locked core query set.")
    }


    # -----------------------------------------------------------------------------
    # Figure 1. Core-query NCS heatmap for objective top 20 candidates
    # -----------------------------------------------------------------------------

    top20 <- candidate_summary |>
      dplyr::arrange(objective_robustness_rank) |>
      dplyr::slice_head(n = 20L) |>
      dplyr::mutate(
        compound_label = display_name(pert_iname, pert),
        compound_label = make.unique(compound_label)
      )

    heatmap_data <- observed_thp1 |>
      dplyr::filter(
        pert %in% top20$pert,
        query_id %in% core_query_ids
      ) |>
      dplyr::left_join(
        top20 |>
          dplyr::select(pert, objective_robustness_rank, compound_label),
        by = "pert"
      ) |>
      dplyr::mutate(
        query_label = factor(
          query_labels[query_id],
          levels = unname(query_labels[core_query_ids])
        ),
        compound_label = factor(
          compound_label,
          levels = rev(top20$compound_label)
        ),
        zero_connectivity = !is.na(NCS) & NCS == 0
      )

    if (nrow(heatmap_data) != 20L * 11L) {
      stop("Top-20 heatmap does not contain exactly 220 observations.")
    }

    p_heatmap <- ggplot2::ggplot(
      heatmap_data,
      ggplot2::aes(x = query_label, y = compound_label, fill = NCS)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.35) +
      ggplot2::geom_point(
        data = heatmap_data |> dplyr::filter(zero_connectivity),
        shape = 4,
        size = 2.2,
        stroke = 0.8,
        color = "black"
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-2.25, 2.25),
        oob = scales::squish,
        name = "THP1 NCS"
      ) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        title = "Connectivity robustness of the top 20 candidates",
        subtitle = "Ordered by integrated average-tie rank; negative NCS indicates predicted reversal; x = zero"
      ) +
      theme_revision(10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 8.5)
      )

    heatmap_paths <- save_plot_pair(
      p_heatmap,
      "candidate_robustness_heatmap_top20",
      width = 12,
      height = 8.2
    )


    # -----------------------------------------------------------------------------
    # Figure 2. Priority-candidate trajectories across core definitions
    # -----------------------------------------------------------------------------

    priority_core <- observed_thp1 |>
      dplyr::filter(
        pert %in% priority_perts,
        query_id %in% core_query_ids
      ) |>
      dplyr::mutate(
        compound = display_name(pert_iname, pert),
        query_label = factor(
          query_labels[query_id],
          levels = unname(query_labels[core_query_ids])
        )
      )

    p_trajectory <- ggplot2::ggplot(
      priority_core,
      ggplot2::aes(
        x = query_label,
        y = rank_all,
        color = compound,
        group = compound
      )
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(shape = NCS == 0), size = 2.8, stroke = 0.8) +
      ggplot2::scale_y_reverse(
        breaks = c(1, 100, 300, 500, 700),
        expand = ggplot2::expansion(mult = c(0.04, 0.08))
      ) +
      ggplot2::scale_color_manual(
        values = c(rosuvastatin = "#0072B2", roflumilast = "#D55E00"),
        name = NULL
      ) +
      ggplot2::scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 4),
        labels = c(`FALSE` = "Non-zero NCS", `TRUE` = "NCS = 0"),
        name = NULL
      ) +
      ggplot2::labs(
        x = NULL,
        y = "THP1 average-tie rank among 1,376 compounds",
        title = "Priority-candidate rank stability across core analyses",
        subtitle = "Lower ranks indicate stronger predicted reversal; zeros retain their average tied rank"
      ) +
      theme_revision(10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    trajectory_paths <- save_plot_pair(
      p_trajectory,
      "priority_candidate_core_rank_trajectories",
      width = 12,
      height = 5.8
    )


    # -----------------------------------------------------------------------------
    # Figure 3. Leave-one-donor-out rank and NCS stability
    # -----------------------------------------------------------------------------

    loo_order <- c("1", "2", "3", "4", "5", "7", "8", "9", "10")

    priority_loo <- observed_thp1 |>
      dplyr::filter(
        pert %in% priority_perts,
        query_role == "donor_sensitivity"
      ) |>
      dplyr::mutate(
        compound = display_name(pert_iname, pert),
        excluded_patient = factor(excluded_patient, levels = loo_order)
      )

    if (nrow(priority_loo) != 18L) {
      stop("Priority LOO table does not contain 18 observations.")
    }

    p_loo_rank <- ggplot2::ggplot(
      priority_loo,
      ggplot2::aes(
        x = excluded_patient,
        y = rank_all,
        color = compound,
        group = compound
      )
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
      ggplot2::scale_color_manual(
        values = c(rosuvastatin = "#0072B2", roflumilast = "#D55E00"),
        name = NULL
      ) +
      ggplot2::labs(
        x = "Excluded patient",
        y = "THP1 average-tie rank among 1,376 compounds",
        title = "Leave-one-donor-out rank stability (average ties)",
        subtitle = "Each point represents a signature recalculated after excluding one donor"
      ) +
      theme_revision(11) +
      ggplot2::theme(legend.position = "bottom")

    p_loo_ncs <- ggplot2::ggplot(
      priority_loo,
      ggplot2::aes(
        x = excluded_patient,
        y = NCS,
        color = compound,
        group = compound
      )
    ) +
      ggplot2::geom_hline(
        yintercept = -0.5,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey40"
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::scale_color_manual(
        values = c(rosuvastatin = "#0072B2", roflumilast = "#D55E00"),
        name = NULL
      ) +
      ggplot2::labs(
        x = "Excluded patient",
        y = "THP1 NCS",
        title = "Leave-one-donor-out connectivity stability",
        subtitle = "Dashed line shows the operational NCS threshold of -0.5"
      ) +
      theme_revision(11) +
      ggplot2::theme(legend.position = "bottom")

    loo_rank_paths <- save_plot_pair(
      p_loo_rank,
      "priority_candidate_loo_rank_stability",
      width = 7.2,
      height = 5.4
    )

    loo_ncs_paths <- save_plot_pair(
      p_loo_ncs,
      "priority_candidate_loo_ncs_stability",
      width = 7.2,
      height = 5.4
    )


    # -----------------------------------------------------------------------------
    # Figure 4. Compound-specific matched-null distributions
    # -----------------------------------------------------------------------------

    null_status <- batch_index$status |>
      dplyr::filter(
        query_role == "specificity_null",
        status == "complete"
      ) |>
      dplyr::arrange(manifest_order)

    if (
      nrow(null_status) != 100L ||
        any(!file.exists(null_status$result_path))
    ) {
      stop("Expected 100 complete matched-null checkpoints.")
    }

    extract_priority_null <- function(result_path, query_id) {
      checkpoint <- readRDS(result_path)
      table <- if (
        is.list(checkpoint) && "result_table" %in% names(checkpoint)
      ) {
        as.data.frame(checkpoint$result_table)
      } else if (isS4(checkpoint)) {
        as.data.frame(signatureSearch::result(checkpoint))
      } else {
        stop("Unsupported null checkpoint structure for ", query_id, ".")
      }

      required <- c("pert", "cell", "type", "NCS")
      missing <- setdiff(required, colnames(table))
      if (length(missing) > 0L) {
        stop(query_id, " lacks required result columns.")
      }

      extracted <- table |>
        dplyr::filter(
          toupper(as.character(cell)) == "THP1",
          as.character(type) == "trt_cp",
          pert %in% priority_perts
        ) |>
        dplyr::transmute(
          query_id = query_id,
          pert = as.character(pert),
          null_NCS = as.numeric(NCS)
        )

      if (nrow(extracted) != 2L || anyDuplicated(extracted$pert)) {
        stop(query_id, " does not contain both unique priority compounds.")
      }

      extracted
    }

    priority_null <- dplyr::bind_rows(lapply(
      seq_len(nrow(null_status)),
      function(i) {
        extract_priority_null(
          null_status$result_path[[i]],
          null_status$query_id[[i]]
        )
      }
    )) |>
      dplyr::left_join(
        candidate_summary |>
          dplyr::filter(pert %in% priority_perts) |>
          dplyr::select(pert, pert_iname, primary_NCS),
        by = "pert"
      ) |>
      dplyr::mutate(compound = display_name(pert_iname, pert))

    if (
      nrow(priority_null) != 200L ||
        any(table(priority_null$pert) != 100L)
    ) {
      stop("Priority matched-null extraction failed.")
    }

    observed_lines <- priority_null |>
      dplyr::distinct(pert, compound, primary_NCS)

    p_null <- ggplot2::ggplot(priority_null, ggplot2::aes(x = null_NCS)) +
      ggplot2::geom_histogram(
        bins = 20,
        fill = "grey72",
        color = "white",
        linewidth = 0.3
      ) +
      ggplot2::geom_vline(
        data = observed_lines,
        ggplot2::aes(xintercept = primary_NCS),
        color = "#D55E00",
        linewidth = 1
      ) +
      ggplot2::facet_wrap(~compound, ncol = 1, scales = "free_y") +
      ggplot2::labs(
        x = "THP1 NCS",
        y = "Matched random signatures",
        title = "Compound-specific matched-null distributions",
        subtitle = "Orange line: observed PRIMARY_150 NCS\n100 matched-null queries per compound; nominal comparison"
      ) +
      theme_revision(11)

    null_paths <- save_plot_pair(
      p_null,
      "priority_candidate_matched_null_distributions",
      width = 7.2,
      height = 6.8
    )


    # -----------------------------------------------------------------------------
    # Figure 5. Primary percentile versus integrated robustness percentile
    # -----------------------------------------------------------------------------

    rank_scatter <- candidate_summary |>
      dplyr::filter(
        !is.na(primary_percentile_all),
        !is.na(median_percentile_all)
      ) |>
      dplyr::mutate(
        compound = display_name(pert_iname, pert),
        priority = tolower(pert_iname) %in% c("rosuvastatin", "roflumilast")
      )

    priority_scatter <- rank_scatter |>
      dplyr::filter(priority) |>
      dplyr::mutate(label_x = 0.06,
                    label_y = ifelse(tolower(pert_iname) == "rosuvastatin", 0.045, 0.10))

    p_rank_scatter <- ggplot2::ggplot(
      rank_scatter,
      ggplot2::aes(x = primary_percentile_all, y = median_percentile_all)
    ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey45"
      ) +
      ggplot2::geom_point(color = "grey65", alpha = 0.55, size = 1.4) +
      ggplot2::geom_point(
        data = priority_scatter,
        ggplot2::aes(color = compound),
        size = 3
      ) +
      ggplot2::geom_segment(
        data = priority_scatter,
        ggplot2::aes(xend = label_x, yend = label_y, color = compound),
        linewidth = 0.3, show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = priority_scatter,
        ggplot2::aes(x = label_x, y = label_y, label = compound, color = compound),
        hjust = 0,
        vjust = 0.5,
        size = 3.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::scale_color_manual(
        values = c(rosuvastatin = "#0072B2", roflumilast = "#D55E00"),
        name = NULL
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        x = "PRIMARY_150 average-tie percentile (1,376 compounds)",
        y = "Median average-tie percentile across 11 core analyses",
        title = "Primary connectivity versus integrated robustness",
        subtitle = paste(
          "Lower percentiles indicate stronger predicted reversal.\n",
          "Average ties; dashed line indicates equal percentiles."
        )
      ) +
      theme_revision(11) +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(8, 35, 8, 8)
      )

    rank_scatter_paths <- save_plot_pair(
      p_rank_scatter,
      "primary_vs_core_percentiles_average_ties",
      width = 9.2,
      height = 6.3
    )



    summary_path <- file.path(RESULT_DIR, "revision_candidate_robustness_average_ties.rds")
    tables <- list(
      candidate_summary_average_ties = candidate_summary,
      core_robustness_average_ties = robustness$core_robustness,
      loo_summary_average_ties = robustness$loo_summary,
      observed_query_results_average_ties = observed_thp1,
      minimum_vs_average_rank_comparison = robustness$minimum_tie_comparison,
      plot_heatmap = heatmap_data, plot_priority_core = priority_core,
      plot_priority_loo = priority_loo, plot_priority_null = priority_null,
      plot_primary_vs_core_percentiles = rank_scatter
    )
    csv_paths <- setNames(file.path(RESULT_DIR, paste0(names(tables), ".csv")), names(tables))
    for (nm in names(tables)) readr::write_csv(tables[[nm]], csv_paths[[nm]], na = "")
    figure_paths <- unname(c(heatmap_paths, trajectory_paths, loo_rank_paths,
                             loo_ncs_paths, null_paths, rank_scatter_paths))
    check(length(figure_paths) == 12L && all(file.exists(figure_paths)) &&
            all(file.info(figure_paths)$size > 0), "One or more figures were not saved.")
    figure_manifest <- data.frame(file = figure_paths, format = tools::file_ext(figure_paths),
                                   md5 = unname(tools::md5sum(figure_paths)))
    readr::write_csv(figure_manifest, file.path(RESULT_DIR, "figure_manifest.csv"))
    check(identical(hashes_before, unname(tools::md5sum(required_files))),
          "An input changed during finalization.")
    cmp <- robustness$minimum_tie_comparison
    changed <- cmp$rank_change_average_minus_min != 0
    first_changed <- if (any(changed)) min(cmp$original_robustness_rank[changed]) else NA_integer_
    decision <- c(
      paste("Build:", SCRIPT_BUILD), paste("Created:", format(Sys.time())),
      "Decision: use average-tie percentiles for the revision-wide integrated robustness summary and figures.",
      "Rationale: assign every exact tied block its midpoint rank rather than its best position.",
      "This rule applies to all 1,376 compounds, all observed queries, and all exact ties, including nonzero ties.",
      "Formula: percentile = (average ascending NCS rank - 1) / 1375; zeros are retained as observations.",
      "The 11 core query definitions, ordering keys after median percentile, and compound universe are unchanged.",
      "All updated rank columns retain their original 04B values under the _minimum suffix.",
      "Submitted nonzero ranks, NCS, WTCS, FDRs, operational-hit counts, and empirical null probabilities are unchanged.",
      "The original 04B summary and 04C figures are preserved; use this run for subsequent revision reporting.",
      paste("Compounds with changed integrated ordinal ranks:", sum(changed)),
      paste("First original rank affected:", first_changed),
      "Source audit: same-sign UP/DOWN enrichment sets bidirectional WTCS to zero; normalization preserves zero.",
      "The enrichment helper can also return zero for zero total matched-gene weight; individual ES components were not saved.",
      "The package trend label puts zero in down; use NCS < 0 to identify negative connectivity.",
      "Average ties can leave visible bands; these represent tied scores, not evidence of biological subgroups.",
      "Operational NCS/FDR thresholds are analysis rules. Empirical null probabilities remain nominal compound-specific values.",
      "These correlated sensitivity analyses support computational prioritization, not clinical efficacy or independent validation.",
      paste("Source audit RDS:", audit_path), paste("Frozen manifest SHA-256:", EXPECTED_LOCK)
    )
    writeLines(decision, file.path(RESULT_DIR, "ranking_method_decision.txt"))
    robustness$status <- "PASS"
    robustness$ranking_decision <- decision
    robustness$audit_rds <- audit_path
    robustness$output_directory <- RESULT_DIR
    robustness$output_paths <- csv_paths
    robustness$figure_directory <- FIGURE_DIR
    robustness$figure_manifest <- figure_manifest
    robustness$input_provenance <- data.frame(input = required_files, md5 = hashes_before)
    robustness$session_info <- utils::sessionInfo()
    saveRDS(robustness, summary_path)
    check(identical(readRDS(summary_path), robustness), "Average-tie summary read-back failed.")
    exported <- list.files(RESULT_DIR, full.names = TRUE)
    check(all(file.info(exported)$size > 0), "An export is empty.")
    readr::write_csv(data.frame(file = basename(exported), md5 = unname(tools::md5sum(exported))),
                     file.path(RESULT_DIR, "output_checksums.csv"))
    writeLines(c("PASS", SCRIPT_BUILD, paste("Summary:", summary_path)), file.path(RESULT_DIR, "_SUCCESS.txt"))
    complete <- TRUE
    cat("\nREVISION ANALYSIS 04E COMPLETE\n")
    cat("Average ties adopted for 1,376 compounds across 11 core analyses.\n")
    cat("Scores, FDRs, null probabilities, and submitted nonzero ranks unchanged: PASS\n")
    cat("Integrated ordinal ranks changed:", sum(changed), "| First original rank affected:", first_changed, "\n")
    print(as.data.frame(candidate_summary[candidate_summary$pert %in% priority_perts,
      c("pert_iname", "objective_robustness_rank_minimum", "objective_robustness_rank", "primary_NCS")]), row.names = FALSE)
    cat("\nSummary RDS:", summary_path, "\nTables:", RESULT_DIR, "\nFigures:", FIGURE_DIR, "\n")
    cat("Figures: 6 designs / 12 PNG-PDF files. Inspect the PNG files.\n")
    cat("Elapsed minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), "\n")
    cat("Next: 04F compound/perturbation metadata audit. No new LINCS queries.\n")
  }
  run_main()
})
