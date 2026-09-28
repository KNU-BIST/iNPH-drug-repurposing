# =============================================================================
# Script: analysis/04D_audit_lincs_zero_scores_and_rank_ties.R
# Build: 04D_ZERO_TIE_AUDIT_V1_20260906
# =============================================================================
# Read-only audit of the 121 locked LINCS queries and the existing 04B summary.
# No scoring engine is executed. No reference database is opened. No upstream
# output, signature, score, p-value, FDR, or existing ranking is overwritten.
#
# Outputs: a NEW run directory inside results/revision/zero_score_rank_tie_audit/
# containing source evidence, score counts, all-compound rank comparisons,
# checksums, session information, and a self-contained audit RDS.
#
# Average ties are evaluated as a sensitivity, NOT adopted as a replacement.
# Both methods use (rank - 1)/(n - 1), the same 11 core queries, and the existing
# secondary ordering: inverse proportion, primary percentile, then BRD ID.
# Exact zeros remain observations; they are not missing or untested compounds.
# A zero WTCS/NCS alone cannot identify its cause without the underlying ES
# components. Source keyword hits are evidence for inspection, not causal proof.
#
# References (interpret the CURRENT installed implementation alongside these):
# https://www.bioconductor.org/packages/release/bioc/vignettes/signatureSearch/inst/doc/signatureSearch.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html
#
# Optional developer check: INPH_04D_SELF_TEST_ONLY=1 runs numerical self-tests
# only, without reading project files or creating output. Normally leave unset.

local({
  BUILD_ID <- "04D_ZERO_TIE_AUDIT_V1_20260906"
  EXPECTED_LOCK <- paste0(
    "1e3a2b827a1429e406303871d547a9dec",
    "07ee4068078c666df69ca82e5428380"
  )
  EXPECTED_CORE <- c(
    "PRIMARY_150", "PRIMARY_50", "PRIMARY_100", "PRIMARY_200",
    "C5_ALLREST_MAX150", "C5_MARKERLOW_MAX150", "MARKER_MAX150",
    "RES03_ALLREST_MAX150", "RES07_ALLREST_MAX150",
    "RES10_ALLREST_MAX150", "P6_INCLUDED_MARKER_MAX150"
  )
  PRIORITY <- c(rosuvastatin = "BRD-K82941592", roflumilast = "BRD-K03194791")

  require_fields <- function(x, fields, label) {
    missing <- setdiff(fields, names(x))
    if (length(missing)) {
      stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  check <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
  }
  same_numeric <- function(x, y) identical(as.numeric(x), as.numeric(y))
  safe_median <- function(x) {
    if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  }
  rank_pair <- function(x) {
    check(is.numeric(x), "NCS must be numeric, not a factor or character.")
    # Inf is not a valid score: fail instead of silently treating it as extreme.
    check(!any(is.infinite(x)), "Infinite NCS encountered in rank calculation.")
    n <- sum(!is.na(x))
    r_min <- as.numeric(base::rank(x, ties.method = "min", na.last = "keep"))
    r_avg <- as.numeric(base::rank(x, ties.method = "average", na.last = "keep"))
    data.frame(
      rank_min = r_min, rank_average = r_avg,
      percentile_min = if (n > 1L) (r_min - 1) / (n - 1) else rep(NA_real_, length(x)),
      percentile_average = if (n > 1L) (r_avg - 1) / (n - 1) else rep(NA_real_, length(x)),
      stringsAsFactors = FALSE
    )
  }
  score_counts <- function(tab, query_id, scope) {
    out <- list(query_id = query_id, population = scope, profiles_n = nrow(tab))
    for (column in c("NCS", "WTCS")) {
      x <- tab[[column]]
      check(is.numeric(x), paste(query_id, column, "is not numeric."))
      finite <- is.finite(x)
      counts <- c(
        negative_n = sum(finite & x < 0),
        zero_n = sum(finite & x == 0),
        positive_n = sum(finite & x > 0),
        missing_n = sum(is.na(x)),
        infinite_n = sum(is.infinite(x)),
        tiny_nonzero_n = sum(finite & x != 0 & abs(x) <= 1e-12)
      )
      check(sum(counts[1:5]) == length(x), "Score-count partition failed.")
      for (nm in names(counts)) out[[paste0(column, "_", nm)]] <- unname(counts[[nm]])
    }
    both_finite <- is.finite(tab$NCS) & is.finite(tab$WTCS)
    out$NCS_zero_WTCS_zero_n <- sum(both_finite & tab$NCS == 0 & tab$WTCS == 0)
    out$NCS_zero_WTCS_nonzero_n <- sum(both_finite & tab$NCS == 0 & tab$WTCS != 0)
    out$NCS_nonzero_WTCS_zero_n <- sum(both_finite & tab$NCS != 0 & tab$WTCS == 0)
    # Placement of the zero block, without sorting or ranking the full database.
    n <- sum(is.finite(tab$NCS))
    z <- out$NCS_zero_n
    neg <- out$NCS_negative_n
    usable <- z > 0L && n > 1L && out$NCS_infinite_n == 0L
    out$zero_min_rank <- if (usable) neg + 1 else NA_real_
    out$zero_average_rank <- if (usable) neg + (z + 1) / 2 else NA_real_
    out$zero_min_percentile <- if (usable) neg / (n - 1) else NA_real_
    out$zero_average_percentile <- if (usable) (neg + (z - 1) / 2) / (n - 1) else NA_real_
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  integrated_order <- function(median_percentile, inverse_proportion, primary_percentile, pert) {
    check(!anyNA(c(median_percentile, inverse_proportion, primary_percentile)),
          "Integrated ranking requires complete finite core scores.")
    o <- base::order(median_percentile, -inverse_proportion, primary_percentile,
                     pert, method = "radix")
    ans <- integer(length(pert))
    ans[o] <- seq_along(o)
    ans
  }
  self_test <- function() {
    r <- rank_pair(c(-2, -1, 0, 0, 0, 1, NA_real_))
    stopifnot(
      same_numeric(r$rank_min, c(1, 2, 3, 3, 3, 6, NA)),
      same_numeric(r$rank_average, c(1, 2, 4, 4, 4, 6, NA)),
      same_numeric(r$percentile_min, c(0, .2, .4, .4, .4, 1, NA)),
      same_numeric(r$percentile_average, c(0, .2, .6, .6, .6, 1, NA))
    )
    z <- rank_pair(rep(0, 5))
    stopifnot(all(z$percentile_min == 0), all(z$percentile_average == .5))
    r <- rank_pair(c(-2, -2, -1, 1, 1))
    stopifnot(same_numeric(r$rank_min, c(1, 1, 3, 4, 4)),
              same_numeric(r$rank_average, c(1.5, 1.5, 3, 4.5, 4.5)))
    r <- rank_pair(c(-2, 0, 1))
    stopifnot(same_numeric(r$rank_min, r$rank_average))
    stopifnot(all(is.na(rank_pair(c(NA_real_, NaN))$percentile_min)),
              is.na(rank_pair(0)$percentile_average))
    cts <- score_counts(data.frame(NCS = c(-1, 0, 1, NA, Inf, -Inf, NaN),
                                  WTCS = c(-1, 0, 1, NA, Inf, -Inf, NaN)), "test", "test")
    stopifnot(cts$NCS_zero_n == 1, cts$NCS_missing_n == 2,
              cts$NCS_infinite_n == 2, is.na(cts$zero_average_rank))
    stopifnot(identical(integrated_order(c(.1, .1, .1), c(1, 1, 1),
                                       c(.2, .2, .1), c("B", "A", "C")), c(3L, 2L, 1L)))
    stopifnot(inherits(try(rank_pair(c(0, Inf)), silent = TRUE), "try-error"))
    invisible(TRUE)
  }

  inspect_installed_source <- function() {
    # Deparse functions; NEVER call the scoring functions themselves.
    ns <- asNamespace("signatureSearch")
    symbols <- ls(ns, all.names = TRUE)
    functions <- list()
    for (nm in symbols) {
      if (bindingIsActive(nm, ns)) next
      obj <- get(nm, envir = ns, inherits = FALSE)
      if (is.function(obj) && identical(environment(obj), ns)) functions[[nm]] <- obj
    }
    code <- lapply(functions, function(f) deparse(f, width.cutoff = 120L))
    roots <- intersect(c("gess_lincs", "lincsEnrich"), names(code))
    check(length(roots) == 2L, "Installed scoring entry points could not be inspected.")
    relevant <- names(code)[vapply(seq_along(code), function(i) {
      grepl("lincs|wtcs|ncs|ksstat", names(code)[i], ignore.case = TRUE) ||
        any(grepl("\\bWTCS\\b|\\bNCS\\b", code[[i]], perl = TRUE))
    }, logical(1))]
    selected <- unique(c(roots, relevant))
    # Follow direct function references, including helper names supplied to apply.
    frontier <- selected
    edges <- list()
    for (depth in seq_len(4L)) {
      next_names <- character()
      for (nm in frontier) {
        called <- intersect(all.names(body(functions[[nm]]), functions = TRUE,
                                      unique = TRUE), names(functions))
        if (length(called)) {
          edges[[length(edges) + 1L]] <- data.frame(caller = nm, callee = called)
          next_names <- c(next_names, called)
        }
      }
      frontier <- setdiff(unique(next_names), selected)
      if (!length(frontier)) break
      selected <- unique(c(selected, frontier))
    }
    evidence_pattern <- paste(
      "WTCS", "NCS", "sign\\(", "== *0", "!= *0", "ifelse", "pmax", "pmin",
      "ESup", "ESdown", "es_up", "es_down", "\\.Call", "\\.C\\(", sep = "|"
    )
    snapshot <- character()
    inventory <- list()
    hits <- list()
    for (nm in selected) {
      lines <- code[[nm]]
      snapshot <- c(snapshot, paste0("### signatureSearch::: ", nm),
                    sprintf("%5d | %s", seq_along(lines), lines), "")
      compiled <- any(grepl("\\.Call|\\.C\\(|\\.External", lines))
      inventory[[nm]] <- data.frame(
        function_name = nm, lines_n = length(lines),
        compiled_call_present = compiled,
        source_sha256 = digest::digest(paste(lines, collapse = "\n"),
                                       algo = "sha256", serialize = FALSE)
      )
      hit <- grep(evidence_pattern, lines, ignore.case = TRUE, perl = TRUE)
      for (i in hit) {
        context <- seq.int(max(1L, i - 2L), min(length(lines), i + 2L))
        hits[[length(hits) + 1L]] <- data.frame(
          function_name = nm, line = i, text = lines[i],
          context = paste(sprintf("%d: %s", context, lines[context]), collapse = " | "),
          interpretation = "Keyword evidence only; review full function and helpers."
        )
      }
    }
    source_hits <- if (length(hits)) do.call(rbind, hits) else data.frame(
      function_name = character(), line = integer(), text = character(),
      context = character(), interpretation = character()
    )
    dll_paths <- list.files(system.file("libs", package = "signatureSearch"),
                           pattern = "\\.(so|dll|dylib)$", recursive = TRUE, full.names = TRUE)
    dlls <- data.frame(path = dll_paths, md5 = unname(tools::md5sum(dll_paths)))
    list(
      package_version = as.character(utils::packageVersion("signatureSearch")),
      package_description = utils::packageDescription("signatureSearch"),
      source_snapshot = snapshot,
      function_inventory = do.call(rbind, inventory),
      keyword_evidence = source_hits,
      call_edges = if (length(edges)) unique(do.call(rbind, edges)) else data.frame(
        caller = character(), callee = character()
      ),
      compiled_libraries = dlls,
      limitation = paste(
        "This records the currently installed R implementation and reachable helpers (four levels).",
        "Compiled internals are not decompiled. Checkpoints do not certify this was their historical version.",
        "ESup/ESdown are not required result fields; individual zero causes may remain unresolved."
      )
    )
  }

  run_audit <- function() {
    started <- Sys.time()
    root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    check(dir.exists(file.path(root, "analysis")) && dir.exists(file.path(root, "data/processed")),
          "Open the project root in RStudio before running 04D.")
    for (pkg in c("digest", "signatureSearch")) {
      check(requireNamespace(pkg, quietly = TRUE), paste("Required installed package missing:", pkg))
    }
    paths <- c(
      manifest = "data/processed/revision_signature_manifest.rds",
      lock = "results/revision/revision_signature_lock.txt",
      batch = "data/processed/revision_lincs_batch_index.rds",
      summary = "data/processed/revision_candidate_robustness_and_null_summary.rds",
      primary = "data/processed/lincs2_gess_result.rds",
      hallmark = "data/processed/hallmark_ifn_sensitivity_result.rds"
    )
    paths <- setNames(file.path(root, paths), names(paths))
    check(all(file.exists(paths)), paste("Missing audit input(s):",
                                        paste(paths[!file.exists(paths)], collapse = ", ")))
    initial_md5 <- unname(tools::md5sum(paths))
    check(!anyNA(initial_md5), "Could not hash the initial audit inputs.")
    manifest_obj <- readRDS(paths[["manifest"]])
    batch <- readRDS(paths[["batch"]])
    summary_obj <- readRDS(paths[["summary"]])
    for (obj in list(manifest_obj, batch, summary_obj)) {
      check(identical(obj$manifest_content_sha256, EXPECTED_LOCK), "An input has the wrong signature lock.")
    }
    require_fields(manifest_obj, c("manifest", "queries", "query_order", "input_provenance"), "Manifest")
    require_fields(summary_obj, c("observed_thp1", "candidate_summary", "core_query_ids"), "04B summary")
    m <- as.data.frame(manifest_obj$manifest)
    require_fields(m, c("query_id", "manifest_order", "query_role", "execution_action",
                        "reference_scope", "signature_sha256", "execution_sha256", "lincs_down_n"), "Manifest table")
    check(identical(manifest_obj$manifest_status, "LOCKED") &&
            identical(manifest_obj$manifest_version, 1L) && nrow(m) == 121L &&
            !anyNA(m$query_id) && !anyDuplicated(m$query_id), "Invalid locked manifest structure.")
    check(identical(as.character(m$query_id), manifest_obj$query_order) &&
            identical(names(manifest_obj$queries), manifest_obj$query_order) &&
            identical(as.integer(m$manifest_order), seq_len(121L)), "Manifest query order mismatch.")
    p <- as.data.frame(manifest_obj$input_provenance)
    require_fields(p, c("input_name", "path", "md5"), "Manifest provenance")
    content <- c("revision_signature_manifest_version=1", "query_count=121",
                 paste0("input_md5:", p$input_name, "=", p$md5),
                 paste0(m$query_id, "=", m$execution_sha256))
    check(identical(digest::digest(paste(content, collapse = "\n"), algo = "sha256",
                                  serialize = FALSE), EXPECTED_LOCK), "Manifest content hash failed.")
    lock_text <- readLines(paths[["lock"]], warn = FALSE)
    check(paste0("manifest_content_sha256=", EXPECTED_LOCK) %in% lock_text,
          "Lock text SHA-256 does not match.")
    md5_line <- grep("^manifest_rds_md5=", lock_text, value = TRUE)
    check(length(md5_line) == 1L && identical(sub("^manifest_rds_md5=", "", md5_line),
          unname(tools::md5sum(paths[["manifest"]]))), "Manifest RDS checksum failed.")
    check(isTRUE(batch$all_locked_queries_complete) && batch$full_reference_complete_n == 19L &&
            batch$null_complete_n == 100L, "All 19 full and 100 null checkpoints must be complete.")
    sv <- batch$subset_validation
    check(is.list(sv) && identical(sv$status, "PASS") &&
            identical(sv$manifest_content_sha256, EXPECTED_LOCK) &&
            length(sv$maximum_absolute_difference[["NCS"]]) == 1L &&
            is.finite(sv$maximum_absolute_difference[["NCS"]]) &&
            sv$maximum_absolute_difference[["NCS"]] <= 1e-12,
          "Saved primary THP1 subset equivalence is not PASS. No test will be rerun by this script.")
    check(identical(sv$primary_full_result_md5, unname(tools::md5sum(paths[["primary"]]))),
          "The saved subset validation refers to a different primary result file.")
    status <- as.data.frame(batch$status)
    require_fields(status, c("query_id", "status", "result_path", "signature_sha256",
                             "execution_sha256", "reference_scope"), "Batch status")
    check(nrow(status) == 121L && !anyDuplicated(status$query_id) &&
            setequal(status$query_id, m$query_id), "Batch/manifest query accounting failed.")
    status <- status[match(m$query_id, status$query_id), , drop = FALSE]
    check(identical(as.character(status$signature_sha256), as.character(m$signature_sha256)) &&
            identical(as.character(status$execution_sha256), as.character(m$execution_sha256)) &&
            identical(as.character(status$reference_scope), as.character(m$reference_scope)),
          "Batch status provenance differs from the locked manifest.")
    reused <- m$query_id %in% c("PRIMARY_150", "HALLMARK_IFNA_NONOVERLAP")
    check(all(status$status[!reused] == "complete") &&
            all(status$status[reused] == "validated_existing"), "Unexpected completion status.")
    core_ids <- as.character(m$query_id[m$query_role %in% c("primary", "biological_sensitivity")])
    check(identical(core_ids, EXPECTED_CORE) && identical(summary_obj$core_query_ids, core_ids),
          "Core query definitions/order differ from the expected 11 analyses.")
    result_paths <- as.character(status$result_path)
    result_paths[m$query_id == "PRIMARY_150"] <- paths[["primary"]]
    result_paths[m$query_id == "HALLMARK_IFNA_NONOVERLAP"] <- paths[["hallmark"]]
    # No guessed path substitutions: stop if a recorded checkpoint cannot be read.
    check(!anyNA(result_paths) && all(file.exists(result_paths)),
          "A recorded checkpoint path is missing. Restore/locate it; do not rerun LINCS.")
    optional <- file.path(root, c(
      "renv.lock", "analysis/03B_run_revision_lincs_batch.R",
      "analysis/04B_candidate_robustness_and_null_summary.R",
      "analysis/04C_visualize_candidate_robustness.R",
      "analysis/04D_audit_lincs_zero_scores_and_rank_ties.R"
    ))
    input_paths <- unique(c(unname(paths), result_paths, optional[file.exists(optional)]))
    before_md5 <- unname(tools::md5sum(input_paths))
    check(!anyNA(before_md5), "Could not hash an audit input.")
    check(identical(initial_md5, before_md5[match(unname(paths), input_paths)]),
          "An initial input changed while it was being loaded.")
    cat("Signature lock validated. Inspecting installed scoring source (no scoring call)...\n")
    source_audit <- inspect_installed_source()

    saved <- as.data.frame(summary_obj$observed_thp1)
    candidates <- as.data.frame(summary_obj$candidate_summary)
    require_fields(saved, c("query_id", "pert", "NCS", "WTCS", "WTCS_FDR",
                            "rank_all", "percentile_all"), "04B observed table")
    require_fields(candidates, c("pert", "pert_iname", "objective_robustness_rank",
                                 "median_percentile_all", "primary_percentile_all",
                                 "inverse_proportion"), "04B candidate table")
    check(nrow(candidates) == 1376L && !anyNA(candidates$pert) && !anyDuplicated(candidates$pert),
          "04B candidate universe is not 1,376 unique identifiers.")
    perts <- sort(as.character(candidates$pert), method = "radix")
    observed_ids <- m$query_id[!grepl("^NULL_MATCHED_", m$query_id)]
    check(nrow(saved) == 21L * 1376L && setequal(unique(saved$query_id), observed_ids),
          "04B observed query accounting failed.")
    candidates <- candidates[match(perts, candidates$pert), , drop = FALSE]
    names_by_pert <- setNames(as.character(candidates$pert_iname), perts)
    counts <- list()
    long <- list()
    reproduction <- list()
    infinite_n <- 0L
    for (i in seq_len(nrow(m))) {
      q <- as.character(m$query_id[i])
      obj <- readRDS(result_paths[i])
      if (!reused[i]) {
        require_fields(obj, c("query_id", "manifest_content_sha256", "signature_sha256",
                              "execution_sha256", "reference_scope", "result_table"), q)
        check(identical(obj$query_id, q) && identical(obj$manifest_content_sha256, EXPECTED_LOCK) &&
                identical(obj$signature_sha256, as.character(m$signature_sha256[i])) &&
                identical(obj$execution_sha256, as.character(m$execution_sha256[i])) &&
                identical(obj$reference_scope, as.character(m$reference_scope[i])),
              paste("Checkpoint provenance mismatch:", q))
        tab <- as.data.frame(obj$result_table)
      } else {
        check(isS4(obj) && all(c("result", "query") %in% methods::slotNames(obj)),
              paste("Expected saved gessResult for", q))
        executed <- methods::slot(obj, "query")
        locked_query <- manifest_obj$queries[[q]]
        check(identical(unname(as.character(executed$upset)), locked_query$up_entrez_in_lincs) &&
                identical(unname(as.character(executed$downset)), locked_query$down_entrez_in_lincs),
              paste("Reused executed query differs from the manifest:", q))
        tab <- as.data.frame(methods::slot(obj, "result"))
      }
      require_fields(tab, c("pert", "cell", "type", "NCS", "WTCS", "WTCS_FDR"), q)
      expected_n <- if (m$reference_scope[i] == "full_lincs2_reference") 136460L else 1376L
      check(nrow(tab) == expected_n, paste("Unexpected reference row count:", q))
      selected <- !is.na(tab$cell) & toupper(as.character(tab$cell)) == "THP1" &
        !is.na(tab$type) & as.character(tab$type) == "trt_cp" &
        !is.na(tab$pert) & grepl("^BRD-", as.character(tab$pert))
      thp <- tab[selected, , drop = FALSE]
      check(nrow(thp) == 1376L && !anyDuplicated(thp$pert) && setequal(thp$pert, perts),
            paste("THP1 universe mismatch:", q))
      thp <- thp[match(perts, thp$pert), , drop = FALSE]
      for (pop in c("saved_reference", "THP1_trt_cp")) {
        ct <- score_counts(if (pop == "saved_reference") tab else thp, q, pop)
        ct$manifest_order <- m$manifest_order[i]
        ct$query_role <- as.character(m$query_role[i])
        ct$reference_scope <- as.character(m$reference_scope[i])
        ct$query_direction <- if (m$lincs_down_n[i] == 0L) "UP_only" else "bidirectional"
        counts[[length(counts) + 1L]] <- ct
      }
      infinite_n <- infinite_n + sum(is.infinite(thp$NCS))
      if (any(is.infinite(thp$NCS))) {
        # Preserve score counts, but do not invent ranks for invalid scores.
        rp <- data.frame(rank_min = rep(NA_real_, 1376L), rank_average = NA_real_,
                         percentile_min = NA_real_, percentile_average = NA_real_)
      } else {
        rp <- rank_pair(thp$NCS)
      }
      long[[q]] <- cbind(data.frame(
        query_id = q, manifest_order = m$manifest_order[i],
        query_role = as.character(m$query_role[i]), pert = perts,
        pert_iname = unname(names_by_pert[perts]), NCS = thp$NCS, WTCS = thp$WTCS,
        is_core = q %in% core_ids, stringsAsFactors = FALSE
      ), rp)
      if (q %in% observed_ids) {
        old <- saved[saved$query_id == q, , drop = FALSE]
        check(nrow(old) == 1376L && !anyDuplicated(old$pert) && setequal(old$pert, perts),
              paste("04B observed compound accounting failed:", q))
        old <- old[match(perts, old$pert), , drop = FALSE]
        score_match <- all(vapply(c("NCS", "WTCS", "WTCS_FDR"), function(nm) {
          same_numeric(old[[nm]], thp[[nm]])
        }, logical(1)))
        rank_match <- same_numeric(old$rank_all, rp$rank_min) &&
          same_numeric(old$percentile_all, rp$percentile_min)
        reproduction[[q]] <- data.frame(query_id = q, scores_equal = score_match,
                                        minimum_ranks_equal = rank_match)
      }
      rm(obj, tab, thp)
      if (i %% 20L == 0L || i == nrow(m)) cat("  Saved queries audited:", i, "/", nrow(m), "\n")
    }
    count_table <- do.call(rbind, counts)
    rank_long <- do.call(rbind, long)
    rownames(rank_long) <- NULL
    rank_long$rank_change_average_minus_min <- rank_long$rank_average - rank_long$rank_min
    rank_long$percentile_change_average_minus_min <- rank_long$percentile_average - rank_long$percentile_min
    reproduction_table <- do.call(rbind, reproduction)
    score_reproduced <- all(reproduction_table$scores_equal)
    check(score_reproduced, "04B scores differ from checkpoints. No audit outputs were finalized; investigate inputs.")
    core <- rank_long[rank_long$is_core, , drop = FALSE]
    complete_core <- all(is.finite(core$NCS))
    rank_reproduced <- all(reproduction_table$minimum_ranks_equal)
    integrated_eligible <- complete_core && rank_reproduced && infinite_n == 0L
    comparison <- data.frame()
    priority <- data.frame()
    top_overlap <- data.frame()
    per_query <- do.call(rbind, lapply(m$query_id, function(q) {
      x <- long[[q]]
      ok <- is.finite(x$rank_min) & is.finite(x$rank_average)
      data.frame(query_id = q, comparable_n = sum(ok),
                 rank_changed_n = sum(x$rank_min[ok] != x$rank_average[ok]),
                 zero_score_n = sum(x$NCS == 0, na.rm = TRUE),
                 nonzero_tie_changed_n = sum(ok & !is.na(x$NCS) & x$NCS != 0 &
                                              x$rank_min != x$rank_average, na.rm = TRUE),
                 maximum_rank_shift = if (any(ok)) max(x$rank_average[ok] - x$rank_min[ok]) else NA_real_)
    }))
    original_order_reproduced <- FALSE
    if (integrated_eligible) {
      # All compound vectors are aligned by BRD ID, not row position in source.
      min_mat <- vapply(core_ids, function(q) long[[q]]$percentile_min, numeric(1376L))
      avg_mat <- vapply(core_ids, function(q) long[[q]]$percentile_average, numeric(1376L))
      ncs_mat <- vapply(core_ids, function(q) long[[q]]$NCS, numeric(1376L))
      med_min <- apply(min_mat, 1L, stats::median)
      med_avg <- apply(avg_mat, 1L, stats::median)
      inv <- rowMeans(ncs_mat < 0)
      orig_rank <- integrated_order(med_min, inv, min_mat[, 1L], perts)
      avg_rank <- integrated_order(med_avg, inv, avg_mat[, 1L], perts)
      original_order_reproduced <- same_numeric(orig_rank, candidates$objective_robustness_rank) &&
        same_numeric(med_min, candidates$median_percentile_all) &&
        same_numeric(min_mat[, 1L], candidates$primary_percentile_all) &&
        same_numeric(inv, candidates$inverse_proportion)
      check(original_order_reproduced,
            "Existing 04B integrated ranking was not exactly reproduced. Stop before interpreting a new ranking.")
      comparison <- data.frame(
        pert = perts, pert_iname = as.character(candidates$pert_iname),
        primary_NCS = ncs_mat[, 1L], core_negative_n = rowSums(ncs_mat < 0),
        core_zero_n = rowSums(ncs_mat == 0), core_positive_n = rowSums(ncs_mat > 0),
        inverse_proportion = inv,
        primary_percentile_min = min_mat[, 1L], primary_percentile_average = avg_mat[, 1L],
        median_core_percentile_min = med_min, median_core_percentile_average = med_avg,
        original_robustness_rank = orig_rank, average_tie_robustness_rank = avg_rank,
        rank_change_average_minus_min = avg_rank - orig_rank,
        stringsAsFactors = FALSE
      )
      # Values are copied for context, never recalculated or reinterpreted as FDR.
      context_cols <- intersect(c("primary_FDR", "operational_hit_n", "loo_operational_hit_n",
                                 "empirical_lower_tail_p"), names(candidates))
      for (nm in context_cols) comparison[[nm]] <- candidates[[nm]]
      comparison <- comparison[order(comparison$original_robustness_rank), , drop = FALSE]
      priority <- comparison[match(unname(PRIORITY), comparison$pert), , drop = FALSE]
      check(!anyNA(priority$pert), "Expected priority compounds not found in the full comparison.")
      top_overlap <- do.call(rbind, lapply(c(10L, 20L, 50L, 100L), function(k) {
        a <- comparison$pert[comparison$original_robustness_rank <= k]
        b <- comparison$pert[comparison$average_tie_robustness_rank <= k]
        overlap <- length(intersect(a, b))
        data.frame(top_k = k, overlap_n = overlap, retained_proportion = overlap / k,
                   jaccard = overlap / length(union(a, b)))
      }))
    }
    audit_status <- if (integrated_eligible && original_order_reproduced) {
      "PASS_AUDIT_COMPLETE_NOT_A_RANKING_REPLACEMENT"
    } else {
      "COUNTS_ONLY_RANK_COMPARISON_WITHHELD"
    }
    after_md5 <- unname(tools::md5sum(input_paths))
    check(identical(before_md5, after_md5), "An input changed during the audit. No outputs finalized.")
    hashes <- data.frame(path = input_paths, md5_before = before_md5, md5_after = after_md5,
                         unchanged = before_md5 == after_md5, stringsAsFactors = FALSE)
    primary_counts <- count_table[count_table$query_id == "PRIMARY_150" &
                                    count_table$population == "THP1_trt_cp", , drop = FALSE]
    notes <- c(
      paste("Build:", BUILD_ID), paste("Status:", audit_status), paste("Signature lock:", EXPECTED_LOCK),
      paste("Current signatureSearch version:", source_audit$package_version),
      "",
      "WHAT THIS AUDIT DOES",
      "Reads all 121 saved queries. Score counts include full saved references and THP1 separately.",
      "Ranks every query under minimum and average ties. Uses all 1,376 THP1 profiles; exact zeros are retained.",
      "Reproduces 04B scores and ranks before comparing integrated ranks across the SAME 11 core analyses.",
      "All-compound ordering uses median core percentile, decreasing inverse proportion, primary percentile, BRD ID.",
      "Positive rank_change_average_minus_min means a WORSE ordinal position under average ties.",
      "Percentiles use (rank-1)/(n-1). Average ties apply to ALL exact ties, not just zeros.",
      "The average percentile of a zero block need not be 50%; it depends on counts on both sides.",
      "An all-zero query gives minimum-tie percentiles 0 and average-tie percentiles 0.5; neither is evidence of reversal.",
      "Integrated comparison is withheld if core scores are missing/nonfinite or original rank checks fail.",
      "",
      "SOURCE EVIDENCE AND LIMITS",
      "installed_scoring_source.txt contains currently installed R functions and reachable helpers, with line numbers.",
      "source_keyword_evidence.csv points to zero/sign/normalization/compiled-call code for manual review.",
      "The documented bidirectional WTCS convention gates discordant enrichment; same-sign ES components can give zero.",
      "Verify the actual branch in the captured implementation; UP-only Hallmark scoring is a separate case.",
      "NCS=0 is not evidence of a missing perturbation, a failed query, or a biologically inert compound.",
      "WTCS/NCS zero concordance alone does NOT prove which branch caused each zero.",
      source_audit$limitation,
      "",
      "INTERPRETATION",
      "This is a ranking sensitivity audit, not independent replication or proof of treatment efficacy.",
      "No scores, FDRs, null probabilities, signatures, or existing summaries have been changed.",
      "Minimum ties put the zero block immediately after negative scores; average ties put it at the block midpoint.",
      "This can move the scatter bands without changing any measured score or providing new biological evidence.",
      "A decision to revise the reportable ranking must be documented separately and applied to every compound.",
      "Empirical null probabilities retained from 04B remain nominal compound-specific values, not adjusted FDR.",
      "The 11 analyses are correlated sensitivity analyses, not 11 independent biological replications.",
      "",
      "REFERENCES",
      "https://www.bioconductor.org/packages/release/bioc/vignettes/signatureSearch/inst/doc/signatureSearch.html",
      "https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html"
    )
    result <- list(
      analysis_id = "revision_04D_zero_score_rank_tie_audit_v1", build_id = BUILD_ID,
      status = audit_status, manifest_content_sha256 = EXPECTED_LOCK,
      created_at = started, core_query_ids = core_ids, source_audit = source_audit,
      score_counts = count_table, per_query_tie_summary = per_query,
      all_query_rank_comparison = rank_long, all_compound_comparison = comparison,
      priority_comparison = priority, top_k_overlap = top_overlap,
      reproduction_checks = reproduction_table,
      original_integrated_order_reproduced = original_order_reproduced,
      complete_finite_core = complete_core,
      input_hashes = hashes, notes = notes, session_info = utils::sessionInfo()
    )
    base_out <- file.path(root, "results/revision/zero_score_rank_tie_audit")
    if (!dir.exists(base_out)) check(dir.create(base_out, recursive = TRUE), "Cannot create audit output parent.")
    out <- tempfile(pattern = paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = base_out)
    check(!dir.exists(out) && dir.create(out), "Cannot create new audit run directory.")
    saved_ok <- FALSE
    on.exit({
      if (!saved_ok) message("Audit export did not finish. Partial files, if any, are in: ", out,
                             "\nDo not treat a directory without _SUCCESS.txt as complete.")
    }, add = TRUE)
    result$output_directory <- out
    csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
    csv(count_table, "query_zero_score_counts.csv")
    csv(per_query, "query_tie_rank_summary.csv")
    csv(rank_long, "all_query_compound_tie_ranks.csv")
    csv(reproduction_table, "original_04B_reproduction_checks.csv")
    if (nrow(comparison)) {
      csv(comparison, "all_compound_rank_tie_comparison.csv")
      csv(priority, "priority_candidate_rank_tie_comparison.csv")
      csv(top_overlap, "top_k_overlap.csv")
      top_union <- comparison$original_robustness_rank <= 20L | comparison$average_tie_robustness_rank <= 20L
      csv(comparison[top_union, , drop = FALSE], "top20_union_rank_tie_comparison.csv")
    }
    csv(source_audit$function_inventory, "source_function_inventory.csv")
    csv(source_audit$keyword_evidence, "source_keyword_evidence.csv")
    csv(source_audit$call_edges, "source_function_calls.csv")
    csv(source_audit$compiled_libraries, "source_compiled_libraries.csv")
    csv(hashes, "input_hashes.csv")
    writeLines(source_audit$source_snapshot, file.path(out, "installed_scoring_source.txt"), useBytes = TRUE)
    writeLines(notes, file.path(out, "audit_notes.txt"), useBytes = TRUE)
    writeLines(capture.output(utils::sessionInfo()), file.path(out, "session_info.txt"), useBytes = TRUE)
    rds_path <- file.path(out, "revision_zero_score_rank_tie_audit.rds")
    saveRDS(result, rds_path, compress = TRUE)
    read_back <- readRDS(rds_path)
    check(identical(read_back, result), "Saved audit RDS read-back failed.")
    exports <- list.files(out, full.names = TRUE)
    check(all(file.info(exports)$size > 0), "An audit export is empty.")
    export_md5 <- unname(tools::md5sum(exports))
    check(!anyNA(export_md5), "Could not checksum the audit exports.")
    csv(data.frame(file = basename(exports), md5 = export_md5), "output_checksums.csv")
    writeLines(c(audit_status, paste("Completed:", format(Sys.time())),
                 "File export complete; see audit_notes.txt for interpretation limits."),
               file.path(out, "_SUCCESS.txt"))
    saved_ok <- TRUE
    cat("\n============================================================\n")
    cat("  REVISION ANALYSIS 04D COMPLETE\n")
    cat("============================================================\n")
    cat("Status:", audit_status, "\n")
    cat("Queries audited:", nrow(m), "| THP1 compounds per query: 1376\n")
    cat("Existing 04B score reproduction:", score_reproduced, "\n")
    cat("Existing 04B minimum-rank reproduction:", rank_reproduced, "\n")
    cat("Existing integrated ranking reproduced:", original_order_reproduced, "\n")
    cat("All input hashes unchanged:", all(hashes$unchanged), "\n")
    cat("Installed scoring functions captured:", nrow(source_audit$function_inventory), "\n")
    cat("\nPRIMARY_150 THP1 zero block:\n")
    print(primary_counts[, c("NCS_negative_n", "NCS_zero_n", "NCS_positive_n", "NCS_missing_n",
                              "zero_min_percentile", "zero_average_percentile")], row.names = FALSE)
    if (nrow(comparison)) {
      cat("\nCompounds changing integrated ordinal rank:",
          sum(comparison$rank_change_average_minus_min != 0), "/", nrow(comparison), "\n")
      print(top_overlap, row.names = FALSE)
      cat("\nPriority candidates (also included in the all-compound comparison):\n")
      print(priority[, c("pert_iname", "primary_NCS", "core_negative_n", "core_zero_n",
                          "original_robustness_rank", "average_tie_robustness_rank",
                          "rank_change_average_minus_min")], row.names = FALSE)
    }
    cat("\nOutput directory:\n", out, "\n", sep = "")
    cat("Audit RDS:\n", rds_path, "\n", sep = "")
    cat("Elapsed minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), "\n")
    cat("No new LINCS query was run. Existing summaries remain unchanged.\n")
    cat("Next: inspect the captured scoring code and all-compound comparison before revising figures.\n")
    invisible(result)
  }

  self_test()
  cat("04D numerical self-tests: PASS\n")
  if (identical(Sys.getenv("INPH_04D_SELF_TEST_ONLY", unset = ""), "1")) {
    cat("Self-test-only mode: no project inputs read and no output files written.\n")
  } else {
    run_audit()
  }
})
