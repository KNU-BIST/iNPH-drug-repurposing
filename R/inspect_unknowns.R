# =============================================================================
# INSPECT UNKNOWNS - answers every open question from the repo audit in one run
# =============================================================================
PROJECT_ROOT <- here::here()

# -----------------------------------------------------------------------------
# 1. Orphaned / unrecognized data/processed files - structure only, not full data
# -----------------------------------------------------------------------------
cat("=================== 1. UNKNOWN PROCESSED FILES - STRUCTURE ===================\n\n")

rds_unknowns <- c(
  "inph_mono_auc_scores.rds", "inph_mono_scored.rds", "lupus_mono_panel_expr.rds",
  "phase5_moa_class_counts.rds", "phase6_jak_inhibitor_recovery.rds",
  "phase6_positive_control_ra.rds"
)
for (f in rds_unknowns) {
  path <- file.path(PROJECT_ROOT, "data", "processed", f)
  cat("---", f, "---\n")
  if (!file.exists(path)) { cat("  (file not found)\n\n"); next }
  obj <- tryCatch(readRDS(path), error = function(e) { cat("  READ ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(obj)) {
    cat("  class:", paste(class(obj), collapse = ", "), "\n")
    if (!is.null(dim(obj))) cat("  dim:", paste(dim(obj), collapse = " x "), "\n")
    else if (!is.null(length(obj))) cat("  length:", length(obj), "\n")
    if (is.data.frame(obj) || is.matrix(obj)) {
      cat("  colnames (first 10):", paste(head(colnames(obj), 10), collapse = ", "), "\n")
      cat("  head:\n"); print(utils::head(obj, 3))
    } else if (is.list(obj)) {
      cat("  names:", paste(names(obj), collapse = ", "), "\n")
      cat("  str (1 level):\n"); str(obj, max.level = 1)
    } else {
      cat("  str:\n"); str(obj)
    }
    rm(obj); gc(verbose = FALSE)
  }
  cat("\n")
}

csv_unknowns <- c("phase3_compounds_named.csv", "phase5_thp1_candidates.csv",
                  "supplementary_table1_full_DE_genes.csv")
for (f in csv_unknowns) {
  path <- file.path(PROJECT_ROOT, "data", "processed", f)
  cat("---", f, "---\n")
  if (!file.exists(path)) { cat("  (file not found)\n\n"); next }
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE, n_max = 5),
                 error = function(e) { cat("  READ ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(df)) {
    cat("  columns:", paste(colnames(df), collapse = ", "), "\n")
    n_total <- length(readr::read_lines(path)) - 1
    cat("  approx row count:", n_total, "\n")
    print(df)
  }
  cat("\n")
}

xlsx_path <- file.path(PROJECT_ROOT, "data", "processed",
                       "Supplementary_Table_2_threshold_sensitivity.xlsx")
cat("--- Supplementary_Table_2_threshold_sensitivity.xlsx ---\n")
if (file.exists(xlsx_path)) {
  sheets <- readxl::excel_sheets(xlsx_path)
  cat("  sheets:", paste(sheets, collapse = ", "), "\n")
  print(readxl::read_excel(xlsx_path, sheet = sheets[1], n_max = 5))
} else cat("  (file not found)\n")
cat("\n")

# -----------------------------------------------------------------------------
# 2. R/ stub files + Diagnostic Script.R + tests - full content (all tiny)
# -----------------------------------------------------------------------------
cat("=================== 2. R/ STUB FILE CONTENTS (full text, all near-empty) ===================\n\n")

stub_files <- c("R/lincs.R", "R/mechanism.R", "R/plotting.R", "R/robustness.R",
                "R/signature.R", "R/specificity.R", "R/translation.R", "R/utils.R",
                "analysis/Diagnostic Script.R",
                "tests/testthat.R", "tests/testthat/test-utils.R")
for (f in stub_files) {
  path <- file.path(PROJECT_ROOT, f)
  cat("----------------------------------------------------------------\n")
  cat("FILE:", f, "\n")
  cat("----------------------------------------------------------------\n")
  if (!file.exists(path)) { cat("(not found)\n\n"); next }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) cat("(empty file, 0 lines)\n")
  else cat(paste(lines, collapse = "\n"), "\n")
  cat("\n")
}

# -----------------------------------------------------------------------------
# 3. Locate every .qmd file with full path + first 15 lines each
# -----------------------------------------------------------------------------
cat("=================== 3. .qmd FILES - LOCATION + PREVIEW ===================\n\n")

qmd_files <- list.files(PROJECT_ROOT, pattern = "\\.qmd$", recursive = TRUE,
                        full.names = TRUE)
qmd_files <- qmd_files[!grepl("renv/library|\\.git", qmd_files)]

if (length(qmd_files) == 0) {
  cat("No .qmd files found (unexpected given the audit count - check exclude patterns).\n")
} else {
  for (f in qmd_files) {
    relf <- sub(paste0("^", PROJECT_ROOT, "/?"), "", f)
    cat("----------------------------------------------------------------\n")
    cat("FILE:", relf, " (", round(file.info(f)$size / 1024, 1), " KB, modified ",
        as.character(file.info(f)$mtime), ")\n", sep = "")
    cat("----------------------------------------------------------------\n")
    lines <- readLines(f, n = 15, warn = FALSE)
    cat(paste(lines, collapse = "\n"), "\n\n")
  }
}

cat("=================== INSPECTION COMPLETE ===================\n")
