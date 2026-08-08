# =============================================================================
# Full Pipeline Reproducibility Test
# Script: run_full_pipeline.R
# =============================================================================
#
# Project
# -------
# Systems biology prioritizes drug candidates for interferon-activated
# monocytes in idiopathic normal pressure hydrocephalus
#
# Author
# ------
# Md. Masudul Haque
#
# Purpose
# -------
# Execute the complete analysis workflow as an unattended reproducibility test
# and verify that the public repository can regenerate the expected analysis
# objects, figures, and manuscript checks from a clean project state.
#
# Scope
# -----
# - Verifies that all canonical analysis scripts and helper files exist.
# - Verifies all staged external inputs required by the finalized workflow.
# - In strict mode, requires data/processed/ to contain only the released
#   canonical duy_monocytes_classical.rds input before execution.
# - Parses every canonical R script before starting the expensive computation,
#   so syntax errors are detected immediately.
# - Runs each analysis script in a separate R process, preventing objects,
#   attached packages, options, or memory state from leaking between scripts.
# - Stops at the first failed script or failed output check.
# - Verifies key output files after each generating script.
# - Optionally requires key outputs to have been regenerated during the current
#   run rather than silently accepting stale files from an earlier session.
# - Writes a master run log, one log per script, a CSV summary, and session
#   information for the runner.
#
# Canonical analysis scripts
# --------------------------
# 01_derive_inph_ifn_signature.R
# 02_cross_disease_contextualization.R
# 03_run_lincs2_connectivity.R
# 04_leave_one_out_robustness.R
# 05_generate_all_figures.R
# 06_audit_manuscript_results.R
# 07_hallmark_ifn_sensitivity.R
# 08_export_full_thp1_compound_ranking.R
# 09_verify_manuscript_tables.R
#
# Execution order
# ---------------
# The execution order below is dependency-aware rather than strictly numeric:
# analytical/sensitivity outputs are completed first, then figures, the
# manuscript audit, and final table verification are run against the complete
# set of generated results.
#
# Usage
# -----
# From the repository root:
#
#   Rscript run_full_pipeline.R
#
# or, from an interactive R session:
#
#   source("run_full_pipeline.R")
#
# Quick preflight only (file names + R syntax; no analysis):
#
#   Rscript run_full_pipeline.R --preflight-only
#
# By default, key generated outputs must be newer than the start time of the
# script that produced them. To permit already-existing outputs:
#
#   Rscript run_full_pipeline.R --allow-existing-outputs
#
# Important
# ---------
# Do not run with --vanilla if the repository uses renv activation through
# .Rprofile. Each child R process intentionally starts normally so that the
# repository's declared R environment can be activated.
# =============================================================================


# =============================================================================
# CONFIGURATION
# =============================================================================

get_runner_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)

  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )

  if (length(file_arg) > 0L) {
    return(
      normalizePath(
        sub("^--file=", "", file_arg[[1]]),
        winslash = "/",
        mustWork = FALSE
      )
    )
  }

  ofile <- tryCatch(
    sys.frame(1)$ofile,
    error = function(e) NULL
  )

  if (!is.null(ofile)) {
    return(
      normalizePath(
        ofile,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }

  NULL
}


RUNNER_PATH <- get_runner_path()

PROJECT_ROOT <- if (!is.null(RUNNER_PATH)) {
  dirname(RUNNER_PATH)
} else {
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}


PIPELINE_STEPS <- list(

  list(
    label = "iNPH IFN signature derivation",
    script = "01_derive_inph_ifn_signature.R",
    expected_outputs = c(
      "data/processed/inph_ifn_signature.rds",
      "data/processed/de_ifn_subtype_results.rds",
      "data/processed/seu_monocyte_subclustered.rds",
      "data/processed/isg_reference_panel.rds"
    )
  ),

  list(
    label = "Cross-disease contextualization",
    script = "02_cross_disease_contextualization.R",
    expected_outputs = c(
      "data/processed/phase2_contextualization.rds"
    )
  ),

  list(
    label = "LINCS2 connectivity analysis",
    script = "03_run_lincs2_connectivity.R",
    expected_outputs = c(
      "data/processed/lincs2_gess_result.rds",
      "data/processed/lincs2_compound_summary.rds",
      "data/processed/lincs2_significant_thp1_reversers.rds",
      "data/processed/lincs2_final_candidates.rds",
      "data/processed/lincs2_top200_thp1_reversal_candidates.csv",
      "data/processed/lincs2_significant_thp1_reversers.csv"
    )
  ),

  list(
    label = "Leave-one-out robustness",
    script = "04_leave_one_out_robustness.R",
    expected_outputs = c(
      "data/processed/loo_signature_stability.rds",
      "data/processed/loo_up_gene_persistence.csv"
    )
  ),

  list(
    label = "Hallmark IFN-alpha sensitivity analysis",
    script = "07_hallmark_ifn_sensitivity.R",
    expected_outputs = c(
      "data/processed/hallmark_ifn_sensitivity_result.rds",
      "data/processed/hallmark_ifn_alpha_pathway_sensitivity_analysis.csv"
    )
  ),

  list(
    label = "Full THP1 compound-ranking export",
    script = "08_export_full_thp1_compound_ranking.R",
    expected_outputs = c(
      "data/processed/lincs2_thp1_nonzero_597_compounds.csv"
    )
  ),

  list(
    label = "Main and supplementary figure generation",
    script = "05_generate_all_figures.R",
    expected_outputs = c(
      "data/processed/pbmc_reconstructed_historical_qc.rds",
      "figures/clean/fig1_combined.png",
      "figures/clean/fig2_combined.png",
      "figures/clean/fig3_combined.png",
      "figures/clean/fig4_combined.png",
      "figures/clean/suppfig_S1_combined.png",
      "figures/clean/suppfig_S4_combined.png",
      "figures/clean/suppfig_S6_mechanism.png"
    )
  ),

  list(
    label = "Manuscript-results audit",
    script = "06_audit_manuscript_results.R",
    expected_outputs = character(0)
  ),

  list(
    label = "Manuscript-table verification",
    script = "09_verify_manuscript_tables.R",
    expected_outputs = character(0)
  )
)


HELPER_FILES <- c(
  "R/data_acquisition_and_loading.R",
  "R/pbmc_clustering_and_annotation.R"
)


# External/canonical inputs staged before the strict fresh-clone test.
#
# EH7297 is intentionally not listed here because Scripts 03 and 07 obtain it
# through ExperimentHub. The GSE292141 helper can download missing GSM files,
# but the strict fresh-clone reproducibility test requires all 20 sample
# triplets to be staged in advance so the test is not confounded by network
# availability during the expensive pipeline run.
CORE_REQUIRED_INPUTS <- c(
  "data/processed/duy_monocytes_classical.rds",
  "data/raw/perez2022_lupus/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad",
  "data/raw/GSE70138_Broad_LINCS_pert_info.txt.gz",
  "data/raw/compoundinfo_beta.txt"
)

DUY_SAMPLE_PREFIXES <- c(
  "GSM8851048_1C",
  "GSM8851049_1B",
  "GSM8851050_2C",
  "GSM8851051_2B",
  "GSM8851052_3C",
  "GSM8851053_3B",
  "GSM8851054_4C",
  "GSM8851055_4B",
  "GSM8851056_5C",
  "GSM8851057_5B",
  "GSM8851058_6C",
  "GSM8851059_6B",
  "GSM8851060_101822C",
  "GSM8851061_101822B",
  "GSM8851062_102522C2",
  "GSM8851063_102522B2",
  "GSM8851064_102522C3",
  "GSM8851065_102522B3",
  "GSM8851066_110122C",
  "GSM8851067_110122B"
)

DUY_REQUIRED_FILES <- unlist(
  lapply(
    DUY_SAMPLE_PREFIXES,
    function(prefix) {
      file.path(
        "data",
        "raw",
        "duy",
        "extracted",
        paste0(
          prefix,
          c(
            "_barcodes.tsv.gz",
            "_features.tsv.gz",
            "_matrix.mtx.gz"
          )
        )
      )
    }
  ),
  use.names = FALSE
)

REQUIRED_INPUTS <- c(
  CORE_REQUIRED_INPUTS,
  DUY_REQUIRED_FILES
)

stopifnot(
  length(DUY_SAMPLE_PREFIXES) == 20L,
  length(DUY_REQUIRED_FILES) == 60L
)


LEGACY_SCRIPT_MAP <- c(
  "01_phase1_signature_derivation.R" = "01_derive_inph_ifn_signature.R",
  "02_lupus_benchmark.R" = "02_cross_disease_contextualization.R",
  "03_lincs_connectivity.R" = "03_run_lincs2_connectivity.R",
  "05_robustness.R" = "04_leave_one_out_robustness.R",
  "06_figures_all.R" = "05_generate_all_figures.R",
  "08_numbers_audit.R" = "06_audit_manuscript_results.R",
  "09_hallmark_ifn_positive_control.R" = "07_hallmark_ifn_sensitivity.R",
  "build_full_597_csv.R" = "08_export_full_thp1_compound_ranking.R",
  "verify_manuscript_tables.R" = "09_verify_manuscript_tables.R"
)


# =============================================================================
# HELPERS
# =============================================================================

project_path <- function(...) {
  file.path(
    PROJECT_ROOT,
    ...
  )
}


timestamp_string <- function() {
  format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )
}


safe_git <- function(args) {
  git_bin <- Sys.which("git")

  if (!nzchar(git_bin)) {
    return(character(0))
  }

  tryCatch(
    system2(
      git_bin,
      c("-C", shQuote(PROJECT_ROOT), args),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0)
  )
}


verify_generated_outputs <- function(
    relative_paths,
    step_start,
    require_fresh_outputs = TRUE,
    freshness_tolerance_sec = 5
) {

  if (length(relative_paths) == 0L) {
    return(
      list(
        ok = TRUE,
        problems = character(0)
      )
    )
  }

  full_paths <- project_path(relative_paths)

  missing <- relative_paths[
    !file.exists(full_paths)
  ]

  if (length(missing) > 0L) {
    return(
      list(
        ok = FALSE,
        problems = paste0(
          "Missing expected output: ",
          missing
        )
      )
    )
  }

  info <- file.info(full_paths)

  empty <- relative_paths[
    is.na(info$size) |
      info$size <= 0
  ]

  problems <- character(0)

  if (length(empty) > 0L) {
    problems <- c(
      problems,
      paste0(
        "Expected output is empty: ",
        empty
      )
    )
  }

  if (require_fresh_outputs) {

    cutoff <- step_start -
      as.difftime(
        freshness_tolerance_sec,
        units = "secs"
      )

    stale <- relative_paths[
      is.na(info$mtime) |
        info$mtime < cutoff
    ]

    if (length(stale) > 0L) {
      problems <- c(
        problems,
        paste0(
          "Expected output was not regenerated during this run: ",
          stale
        )
      )
    }
  }

  list(
    ok = length(problems) == 0L,
    problems = problems
  )
}


parse_r_file <- function(path) {
  tryCatch(
    {
      parse(
        file = path,
        keep.source = FALSE
      )
      NULL
    },
    error = function(e) {
      conditionMessage(e)
    }
  )
}


make_step_wrapper <- function(
    script_path,
    wrapper_path
) {

  source_call <- paste0(
    "source(",
    deparse(
      normalizePath(
        script_path,
        winslash = "/",
        mustWork = TRUE
      )
    ),
    ", chdir = FALSE)"
  )

  wrapper_lines <- c(
    "options(warn = 1)",
    "cat('Child R process started:', format(Sys.time()), '\\n')",
    source_call,
    "cat('\\n--- CHILD SESSION INFO ---\\n')",
    "print(sessionInfo())",
    "cat('\\nChild R process completed:', format(Sys.time()), '\\n')"
  )

  writeLines(
    wrapper_lines,
    wrapper_path,
    useBytes = TRUE
  )
}


# =============================================================================
# MAIN RUNNER
# =============================================================================

run_full_pipeline <- function(
    preflight_only = FALSE,
    require_fresh_outputs = TRUE
) {

  old_wd <- getwd()

  on.exit(
    setwd(old_wd),
    add = TRUE
  )

  setwd(PROJECT_ROOT)


  # ---------------------------------------------------------------------------
  # Logging setup
  # ---------------------------------------------------------------------------

  run_id <- timestamp_string()

  log_dir <- project_path(
    "logs",
    "reproducibility",
    paste0(
      "run_",
      run_id
    )
  )

  dir.create(
    log_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  master_log <- file.path(
    log_dir,
    "pipeline_master.log"
  )

  summary_csv <- file.path(
    log_dir,
    "pipeline_summary.csv"
  )

  runner_session_file <- file.path(
    log_dir,
    "runner_sessionInfo.txt"
  )

  log_con <- file(
    master_log,
    open = "wt"
  )

  sink(
    log_con,
    split = TRUE
  )

  sink(
    log_con,
    type = "message"
  )

  on.exit(
    {
      while (sink.number(type = "message") > 0L) {
        sink(type = "message")
      }

      while (sink.number() > 0L) {
        sink()
      }

      if (isOpen(log_con)) {
        close(log_con)
      }
    },
    add = TRUE
  )


  # ---------------------------------------------------------------------------
  # Run metadata
  # ---------------------------------------------------------------------------

  cat(
    "===============================================================================\n"
  )

  cat(
    "FULL PIPELINE REPRODUCIBILITY TEST\n"
  )

  cat(
    "===============================================================================\n"
  )

  cat(
    "Project root: ",
    PROJECT_ROOT,
    "\n",
    sep = ""
  )

  cat(
    "Start time: ",
    format(Sys.time()),
    "\n",
    sep = ""
  )

  cat(
    "Fresh-output verification: ",
    if (require_fresh_outputs) "ENABLED" else "DISABLED",
    "\n",
    sep = ""
  )

  cat(
    "Preflight-only mode: ",
    if (preflight_only) "YES" else "NO",
    "\n\n",
    sep = ""
  )


  git_commit <- safe_git(
    c(
      "rev-parse",
      "HEAD"
    )
  )

  if (length(git_commit) > 0L) {
    cat(
      "Git commit: ",
      git_commit[[1]],
      "\n",
      sep = ""
    )
  } else {
    cat(
      "Git commit: unavailable\n"
    )
  }


  git_status <- safe_git(
    c(
      "status",
      "--porcelain"
    )
  )

  if (length(git_status) == 0L) {
    cat(
      "Git working tree: clean or Git status unavailable\n"
    )
  } else {
    cat(
      "Git working tree contains uncommitted changes:\n"
    )

    cat(
      paste(
        git_status,
        collapse = "\n"
      ),
      "\n"
    )
  }


  if (file.exists(project_path("renv.lock"))) {
    cat(
      "renv.lock: found\n"
    )
  } else {
    cat(
      "WARNING: renv.lock was not found. Package-version reproducibility ",
      "cannot be verified from this clone.\n",
      sep = ""
    )
  }


  if (file.exists(project_path(".Rprofile"))) {
    cat(
      ".Rprofile: found\n"
    )
  } else {
    cat(
      "NOTE: .Rprofile was not found. If the project relies on renv automatic ",
      "activation, restore/activation must be handled explicitly.\n",
      sep = ""
    )
  }


  if (
    file.exists(project_path("renv.lock")) &&
    requireNamespace("renv", quietly = TRUE)
  ) {

    cat(
      "\n--- renv::status() ---\n"
    )

    tryCatch(
      renv::status(),
      error = function(e) {
        cat(
          "WARNING: renv::status() could not be completed: ",
          conditionMessage(e),
          "\n",
          sep = ""
        )
      }
    )

    cat(
      "--- end renv::status() ---\n"
    )
  }


  writeLines(
    capture.output(
      sessionInfo()
    ),
    runner_session_file,
    useBytes = TRUE
  )


  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: repository structure
  # ---------------------------------------------------------------------------

  cat(
    "\n===============================================================================\n"
  )

  cat(
    "PRE-FLIGHT CHECKS\n"
  )

  cat(
    "===============================================================================\n\n"
  )


  analysis_scripts <- vapply(
    PIPELINE_STEPS,
    function(x) {
      file.path(
        "analysis",
        x$script
      )
    },
    character(1)
  )


  required_code_files <- c(
    analysis_scripts,
    HELPER_FILES
  )


  missing_code <- required_code_files[
    !file.exists(
      project_path(
        required_code_files
      )
    )
  ]


  if (length(missing_code) > 0L) {

    cat(
      "FAILED: required code file(s) are missing:\n"
    )

    cat(
      paste0(
        "  - ",
        missing_code,
        collapse = "\n"
      ),
      "\n\n"
    )


    # Give a useful rename hint when an old file is still present.
    legacy_present <- names(
      LEGACY_SCRIPT_MAP
    )[
      file.exists(
        project_path(
          "analysis",
          names(
            LEGACY_SCRIPT_MAP
          )
        )
      )
    ]


    if (length(legacy_present) > 0L) {

      cat(
        "Legacy script name(s) detected. Rename them before testing:\n"
      )

      for (old_name in legacy_present) {
        cat(
          "  ",
          old_name,
          " -> ",
          LEGACY_SCRIPT_MAP[[old_name]],
          "\n",
          sep = ""
        )
      }

      cat(
        "\n"
      )
    }


    stop(
      "Pre-flight failed because canonical code files are missing."
    )
  }


  cat(
    "Canonical analysis scripts: all present\n"
  )

  cat(
    "Required helper files: all present\n"
  )


  # Warn about duplicate legacy script names even when the canonical names
  # are also present. Leaving both versions in analysis/ can confuse users
  # about which script is authoritative.
  legacy_present <- names(
    LEGACY_SCRIPT_MAP
  )[
    file.exists(
      project_path(
        "analysis",
        names(
          LEGACY_SCRIPT_MAP
        )
      )
    )
  ]

  if (length(legacy_present) > 0L) {

    cat(
      "\nWARNING: legacy script name(s) are still present in analysis/:\n"
    )

    for (old_name in legacy_present) {
      cat(
        "  ",
        old_name,
        " -> canonical name: ",
        LEGACY_SCRIPT_MAP[[old_name]],
        "\n",
        sep = ""
      )
    }

    cat(
      "Remove or archive these duplicate legacy files before the public GitHub push.\n"
    )
  }


  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: syntax parsing
  # ---------------------------------------------------------------------------

  syntax_failures <- list()


  for (rel_path in required_code_files) {

    full_path <- project_path(
      rel_path
    )

    parse_error <- parse_r_file(
      full_path
    )

    if (!is.null(parse_error)) {
      syntax_failures[[rel_path]] <- parse_error
    }
  }


  if (length(syntax_failures) > 0L) {

    cat(
      "\nFAILED: R syntax error(s) detected before execution:\n"
    )

    for (nm in names(syntax_failures)) {
      cat(
        "\n",
        nm,
        "\n  ",
        syntax_failures[[nm]],
        "\n",
        sep = ""
      )
    }

    stop(
      "Pre-flight failed because one or more R files could not be parsed."
    )
  }


  cat(
    "R syntax parse check: all canonical scripts/helpers passed\n"
  )


  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: required non-automated inputs
  # ---------------------------------------------------------------------------

  missing_inputs <- REQUIRED_INPUTS[
    !file.exists(
      project_path(
        REQUIRED_INPUTS
      )
    )
  ]


  if (length(missing_inputs) > 0L) {

    cat(
      "\nFAILED: required input file(s) are missing:\n"
    )

    cat(
      paste0(
        "  - ",
        missing_inputs,
        collapse = "\n"
      ),
      "\n\n"
    )

    cat(
      "Stage the canonical 22,197-cell monocyte object, the Perez comparator ",
      "h5ad, the GSE70138 perturbation metadata, compoundinfo_beta.txt, and ",
      "all 20 GSE292141 GSM triplets (60 files) before the strict fresh-clone ",
      "run. See README.md for the documented input locations.\n",
      sep = ""
    )

    stop(
      "Pre-flight failed because required non-automated inputs are missing."
    )
  }


  input_info <- file.info(
    project_path(
      REQUIRED_INPUTS
    )
  )

  bad_inputs <- REQUIRED_INPUTS[
    is.na(input_info$size) |
      input_info$size <= 0
  ]


  if (length(bad_inputs) > 0L) {

    stop(
      "Required input file(s) are empty: ",
      paste(
        bad_inputs,
        collapse = ", "
      )
    )
  }


  cat(
    "Required canonical/raw inputs: present and non-empty\n"
  )

  core_input_md5 <- tools::md5sum(
    project_path(
      CORE_REQUIRED_INPUTS
    )
  )

  cat(
    "Core required-input MD5 checksums:\n"
  )

  for (j in seq_along(CORE_REQUIRED_INPUTS)) {
    cat(
      "  ",
      CORE_REQUIRED_INPUTS[[j]],
      "  ",
      unname(core_input_md5[[j]]),
      "\n",
      sep = ""
    )
  }

  cat(
    "GSE292141 staged sample sets: ",
    length(DUY_SAMPLE_PREFIXES),
    " | required non-empty triplet files: ",
    length(DUY_REQUIRED_FILES),
    "\n",
    sep = ""
  )


  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: strict fresh-clone state
  # ---------------------------------------------------------------------------

  if (require_fresh_outputs) {

    allowed_processed_basename <- "duy_monocytes_classical.rds"

    processed_dir <- project_path(
      "data",
      "processed"
    )

    if (dir.exists(processed_dir)) {
      processed_entries <- list.files(
        processed_dir,
        all.files = FALSE,
        full.names = FALSE,
        recursive = FALSE
      )

      unexpected_processed <- setdiff(
        processed_entries,
        allowed_processed_basename
      )
    } else {
      unexpected_processed <- character(0)
    }

    if (length(unexpected_processed) > 0L) {
      cat(
        "\nFAILED: strict fresh-clone mode found generated file(s) already ",
        "present in data/processed/:\n",
        sep = ""
      )

      cat(
        paste0(
          "  - data/processed/",
          unexpected_processed,
          collapse = "\n"
        ),
        "\n\n",
        sep = ""
      )

      cat(
        "For the strict reproducibility test, data/processed/ must contain ",
        "only duy_monocytes_classical.rds before the pipeline starts.\n",
        sep = ""
      )

      stop(
        "Pre-flight failed because generated processed outputs already exist."
      )
    }

    cat(
      "Strict fresh-clone processed state: PASS ",
      "(only duy_monocytes_classical.rds is present)\n",
      sep = ""
    )
  }


  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: Rscript executable
  # ---------------------------------------------------------------------------

  rscript_bin <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") {
      "Rscript.exe"
    } else {
      "Rscript"
    }
  )


  if (!file.exists(rscript_bin)) {
    stop(
      "Rscript executable not found at: ",
      rscript_bin
    )
  }


  cat(
    "Child-process executable: ",
    rscript_bin,
    "\n",
    sep = ""
  )


  if (preflight_only) {

    cat(
      "\nPRE-FLIGHT PASSED. No analysis scripts were executed.\n"
    )

    cat(
      "Log directory: ",
      log_dir,
      "\n",
      sep = ""
    )

    return(
      invisible(
        list(
          status = "PREFLIGHT_OK",
          log_dir = log_dir
        )
      )
    )
  }


  # ---------------------------------------------------------------------------
  # RESULTS TABLE
  # ---------------------------------------------------------------------------

  results <- data.frame(
    execution_order = seq_along(
      PIPELINE_STEPS
    ),
    script = vapply(
      PIPELINE_STEPS,
      `[[`,
      character(1),
      "script"
    ),
    label = vapply(
      PIPELINE_STEPS,
      `[[`,
      character(1),
      "label"
    ),
    status = "NOT_RUN",
    exit_code = NA_integer_,
    elapsed_min = NA_real_,
    output_check = NA_character_,
    step_log = NA_character_,
    stringsAsFactors = FALSE
  )


  overall_start <- Sys.time()


  # ---------------------------------------------------------------------------
  # EXECUTE EACH SCRIPT IN A CLEAN R PROCESS
  # ---------------------------------------------------------------------------

  for (i in seq_along(PIPELINE_STEPS)) {

    step <- PIPELINE_STEPS[[i]]

    script_name <- step$script

    script_path <- project_path(
      "analysis",
      script_name
    )

    step_log <- file.path(
      log_dir,
      paste0(
        sprintf(
          "%02d",
          i
        ),
        "_",
        tools::file_path_sans_ext(
          script_name
        ),
        ".log"
      )
    )

    wrapper_path <- tempfile(
      pattern = paste0(
        "pipeline_step_",
        i,
        "_"
      ),
      fileext = ".R"
    )

    on.exit(
      unlink(
        wrapper_path,
        force = TRUE
      ),
      add = TRUE
    )


    make_step_wrapper(
      script_path = script_path,
      wrapper_path = wrapper_path
    )


    cat(
      "\n",
      strrep(
        "=",
        79
      ),
      "\n",
      sep = ""
    )

    cat(
      "STARTING [",
      i,
      "/",
      length(PIPELINE_STEPS),
      "] ",
      script_name,
      "\n",
      sep = ""
    )

    cat(
      "Purpose: ",
      step$label,
      "\n",
      sep = ""
    )

    cat(
      "Start time: ",
      format(Sys.time()),
      "\n",
      sep = ""
    )

    cat(
      "Step log: ",
      step_log,
      "\n",
      sep = ""
    )

    cat(
      strrep(
        "=",
        79
      ),
      "\n\n",
      sep = ""
    )


    step_start <- Sys.time()


    exit_code <- tryCatch(
      system2(
        command = rscript_bin,
        args = shQuote(wrapper_path),
        stdout = step_log,
        stderr = step_log,
        wait = TRUE
      ),
      error = function(e) {

        cat(
          "Runner-level error while launching ",
          script_name,
          ": ",
          conditionMessage(e),
          "\n",
          sep = ""
        )

        999L
      }
    )


    step_elapsed <- as.numeric(
      difftime(
        Sys.time(),
        step_start,
        units = "mins"
      )
    )


    results$exit_code[i] <- as.integer(
      exit_code
    )

    results$elapsed_min[i] <- round(
      step_elapsed,
      2
    )

    results$step_log[i] <- normalizePath(
      step_log,
      winslash = "/",
      mustWork = FALSE
    )


    if (!identical(
      as.integer(exit_code),
      0L
    )) {

      results$status[i] <- "FAILED"

      results$output_check[i] <- "NOT_CHECKED"


      cat(
        "\nFAILED: ",
        script_name,
        " exited with status ",
        exit_code,
        ".\n",
        sep = ""
      )


      if (file.exists(step_log)) {

        failure_lines <- readLines(
          step_log,
          warn = FALSE
        )

        cat(
          "\nLast lines from the failed step log:\n"
        )

        cat(
          paste(
            utils::tail(
              failure_lines,
              40
            ),
            collapse = "\n"
          ),
          "\n"
        )
      }


      cat(
        "\nSTOPPING PIPELINE. Later steps were not attempted.\n"
      )

      break
    }


    output_check <- verify_generated_outputs(
      relative_paths = step$expected_outputs,
      step_start = step_start,
      require_fresh_outputs = require_fresh_outputs
    )


    if (!output_check$ok) {

      results$status[i] <- "OUTPUT_CHECK_FAILED"

      results$output_check[i] <- paste(
        output_check$problems,
        collapse = " | "
      )


      cat(
        "\nSCRIPT COMPLETED BUT OUTPUT VERIFICATION FAILED:\n"
      )

      cat(
        paste0(
          "  - ",
          output_check$problems,
          collapse = "\n"
        ),
        "\n"
      )

      cat(
        "\nSTOPPING PIPELINE. This usually indicates a stale-file dependency, ",
        "an output-name mismatch, or a script that did not write the expected ",
        "artifact.\n",
        sep = ""
      )

      break
    }


    results$status[i] <- "OK"

    results$output_check[i] <- if (
      length(step$expected_outputs) > 0L
    ) {
      "PASSED"
    } else {
      "NOT_APPLICABLE"
    }


    cat(
      "\nCOMPLETE: ",
      script_name,
      " (",
      round(
        step_elapsed,
        2
      ),
      " min)\n",
      sep = ""
    )


    if (length(step$expected_outputs) > 0L) {

      cat(
        "Expected-output verification: PASSED\n"
      )

      for (output_path in step$expected_outputs) {
        cat(
          "  OK: ",
          output_path,
          "\n",
          sep = ""
        )
      }
    }
  }


  # ---------------------------------------------------------------------------
  # FINAL SUMMARY
  # ---------------------------------------------------------------------------

  overall_elapsed <- as.numeric(
    difftime(
      Sys.time(),
      overall_start,
      units = "hours"
    )
  )


  utils::write.csv(
    results,
    summary_csv,
    row.names = FALSE
  )


  cat(
    "\n\n",
    strrep(
      "=",
      79
    ),
    "\n",
    sep = ""
  )

  cat(
    "PIPELINE RUN SUMMARY\n"
  )

  cat(
    strrep(
      "=",
      79
    ),
    "\n\n",
    sep = ""
  )

  print(
    results,
    row.names = FALSE
  )

  cat(
    sprintf(
      "\nTotal elapsed: %.2f hours\n",
      overall_elapsed
    )
  )

  cat(
    "End time: ",
    format(Sys.time()),
    "\n",
    sep = ""
  )

  cat(
    "Master log: ",
    master_log,
    "\n",
    sep = ""
  )

  cat(
    "Summary CSV: ",
    summary_csv,
    "\n",
    sep = ""
  )

  cat(
    "Per-script logs: ",
    log_dir,
    "\n",
    sep = ""
  )


  all_ok <- all(
    results$status == "OK"
  )


  if (all_ok) {

    cat(
      "\nALL CANONICAL SCRIPTS COMPLETED SUCCESSFULLY.\n"
    )

    cat(
      "Next: review 06_audit_manuscript_results.R and ",
      "09_verify_manuscript_tables.R logs for the numerical comparison ",
      "against the locked manuscript values.\n",
      sep = ""
    )

  } else {

    cat(
      "\nPIPELINE DID NOT COMPLETE SUCCESSFULLY.\n"
    )

    cat(
      "Review the first FAILED / OUTPUT_CHECK_FAILED / NOT_RUN row and its ",
      "corresponding step log.\n",
      sep = ""
    )
  }


  invisible(
    list(
      results = results,
      log_dir = log_dir,
      master_log = master_log,
      summary_csv = summary_csv,
      all_ok = all_ok
    )
  )
}


# =============================================================================
# ENTRY POINT
# =============================================================================

cli_args <- commandArgs(
  trailingOnly = TRUE
)

PRE_FLIGHT_ONLY <- "--preflight-only" %in%
  cli_args

REQUIRE_FRESH_OUTPUTS <- !(
  "--allow-existing-outputs" %in%
    cli_args
)


pipeline_result <- run_full_pipeline(
  preflight_only = PRE_FLIGHT_ONLY,
  require_fresh_outputs = REQUIRE_FRESH_OUTPUTS
)
