# =============================================================================
# REPO AUDIT - run before any curation/renaming decisions
# =============================================================================
# Purpose: full inventory of the project folder, with automatic flags for the
# specific hazard patterns already found in this project this session:
#   - multiple versions of the same script sitting side by side (v1/v2/v3,
#     V2/V3/V4) with no clear signal which one is canonical
#   - retired-pipeline artifacts (targets, old CLUEio annotation script,
#     Book1.xlsx) that should be archived or clearly marked non-canonical
#   - scripts not yet matching the master file's 01-08 naming scheme
#   - large files sitting in a location that would bloat a git repo if
#     accidentally tracked (raw data, caches)
#   - no assumptions are made about what to DELETE - this only reports
# =============================================================================

PROJECT_ROOT <- here::here()
cat("Project root:", PROJECT_ROOT, "\n\n")

# Exclude noisy/irrelevant dirs from the walk
EXCLUDE_DIRS <- c("renv/library", "\\.git", "\\.Rproj\\.user", "node_modules",
                  "experiment_hub_cache")
exclude_pattern <- paste(EXCLUDE_DIRS, collapse = "|")

all_files <- list.files(PROJECT_ROOT, recursive = TRUE, full.names = TRUE,
                        all.files = FALSE, no.. = TRUE)
all_files <- all_files[!grepl(exclude_pattern, all_files)]

info <- file.info(all_files)
df <- data.frame(
  path    = all_files,
  relpath = sub(paste0("^", PROJECT_ROOT, "/?"), "", all_files),
  ext     = tolower(tools::file_ext(all_files)),
  size_mb = round(info$size / 1e6, 2),
  mtime   = info$mtime,
  stringsAsFactors = FALSE
)
df <- df[order(df$relpath), ]

cat("=================== 1. FULL FILE COUNT BY TYPE ===================\n")
print(df |> dplyr::count(ext, sort = TRUE))
cat("\nTotal files scanned:", nrow(df), "\n\n")

cat("=================== 2. R/ HELPER FUNCTIONS (sourced, not run) ===================\n")
r_helpers <- df[grepl("^R/", df$relpath) & df$ext == "r", ]
print(r_helpers[, c("relpath", "size_mb", "mtime")], row.names = FALSE)
cat("\n")

cat("=================== 3. TOP-LEVEL / analysis SCRIPTS ===================\n")
# Any .R file NOT under R/, renv/, or figures/ - i.e. the actual pipeline scripts
top_scripts <- df[df$ext == "r" & !grepl("^R/", df$relpath), ]
print(top_scripts[, c("relpath", "size_mb", "mtime")], row.names = FALSE)
cat("\n")

cat("=================== 4. VERSIONED-DUPLICATE DETECTION ===================\n")
cat("Files whose name matches a versioning pattern (v1/v2/_V3/etc). Groups\n")
cat("with >1 file are exactly the 'which one is canonical' hazard that broke\n")
cat("script 02/03 earlier this project. Manually confirm before deleting -\n")
cat("this only flags the pattern, it does not know which version is real.\n\n")

ver_pattern <- "_[vV][0-9]+\\.[Rr]$|_[vV][0-9]+\\.txt$"
versioned <- df[grepl(ver_pattern, df$relpath), ]
if (nrow(versioned) > 0) {
  # Strip version suffix to find the "family" each file belongs to
  versioned$family <- sub(ver_pattern, "", versioned$relpath)
  # Also catch bare (unversioned) files that share the same base name/stem
  base_names <- unique(versioned$family)
  for (fam in base_names) {
    cat("--- Family:", basename(fam), "---\n")
    related <- df[grepl(paste0("^", tools::file_path_sans_ext(fam)), 
                        tools::file_path_sans_ext(df$relpath)) &
                  df$ext %in% c("r"), ]
    print(related[, c("relpath", "size_mb", "mtime")], row.names = FALSE)
    cat("\n")
  }
} else {
  cat("None found by pattern match - check section 3 above manually for\n")
  cat("near-duplicate names without a clean version suffix.\n\n")
}

cat("=================== 5. RETIRED-PIPELINE ARTIFACTS ===================\n")
cat("Known-retired items from this session. Recommend archiving (not\n")
cat("deleting outright, in case content needs referencing) rather than\n")
cat("leaving in an ambiguous state next to the active pipeline.\n\n")

retired_patterns <- c(
  "Book1\\.xlsx", "CLUEio", "_targets", "targets\\.R", "supp_table3_thp1_ranked_moa"
)
for (pat in retired_patterns) {
  hits <- df[grepl(pat, df$relpath, ignore.case = TRUE), ]
  if (nrow(hits) > 0) {
    cat("Pattern '", pat, "':\n", sep = "")
    print(hits[, c("relpath", "size_mb", "mtime")], row.names = FALSE)
    cat("\n")
  }
}

cat("=================== 6. NAMING vs THE 01-08 CANONICAL SCHEME ===================\n")
canonical <- c(
  "01" = "signature derivation",
  "02" = "lupus cross-disease benchmark",
  "03" = "LINCS connectivity",
  "04" = "cell-line profiles",
  "05" = "robustness (LOO)",
  "06" = "main figures",
  "07" = "supplementary figures",
  "08" = "numbers audit"
)
cat("Expected canonical files (per project master tracking doc):\n")
for (n in names(canonical))
  cat("  ", n, "_*.R  ->  ", canonical[n], "\n", sep = "")
cat("\nActual top-level .R files matching a leading NN_ pattern:\n")
numbered <- top_scripts[grepl("^[0-9]{2}_", basename(top_scripts$relpath)), ]
print(numbered[, c("relpath", "size_mb")], row.names = FALSE)
cat("\nTop-level .R files NOT matching NN_ naming (candidates to rename or archive):\n")
unnumbered <- top_scripts[!grepl("^[0-9]{2}_", basename(top_scripts$relpath)), ]
print(unnumbered[, c("relpath", "size_mb")], row.names = FALSE)
cat("\n")

cat("=================== 7. LARGE FILES (>50 MB) - CHECK GIT TRACKING ===================\n")
large <- df[df$size_mb > 50, ]
large <- large[order(-large$size_mb), ]
if (nrow(large) > 0) {
  print(large[, c("relpath", "size_mb")], row.names = FALSE)
  cat("\nNone of these should be committed to git directly. Confirm .gitignore\n")
  cat("covers their parent directories (see section 9).\n\n")
} else {
  cat("None found.\n\n")
}

cat("=================== 8. data/raw AND data/processed CONTENTS ===================\n")
for (sub in c("data/raw", "data/processed", "data/raw/duy",
              "data/raw/perez2022_lupus", "figures")) {
  full <- file.path(PROJECT_ROOT, sub)
  if (dir.exists(full)) {
    items <- list.files(full, recursive = FALSE)
    cat("--- ", sub, " (", length(items), " items) ---\n", sep = "")
    if (length(items) <= 25) print(items) else {
      print(items[1:25]); cat("  ... and", length(items) - 25, "more\n")
    }
    cat("\n")
  }
}

cat("=================== 9. .gitignore / .Renviron STATUS ===================\n")
gi_path <- file.path(PROJECT_ROOT, ".gitignore")
if (file.exists(gi_path)) {
  cat(".gitignore exists. Contents:\n")
  cat(paste(readLines(gi_path), collapse = "\n"), "\n\n")
} else {
  cat(".gitignore NOT FOUND. Needed before any git push - data/, figures/,\n")
  cat("experiment_hub_cache/, and .Renviron should not be tracked.\n\n")
}

renv_path <- file.path(PROJECT_ROOT, ".Renviron")
cat(".Renviron exists:", file.exists(renv_path),
    if (file.exists(renv_path)) paste0(" (", length(readLines(renv_path)), " lines - confirm empty or minimal per earlier cleanup)") else "", "\n\n")

cat("=================== 10. README STATUS ===================\n")
readme_hits <- df[grepl("^readme", df$relpath, ignore.case = TRUE), ]
if (nrow(readme_hits) > 0) {
  print(readme_hits[, c("relpath", "size_mb", "mtime")], row.names = FALSE)
} else {
  cat("No README found at project root. Needed before publishing: run order,\n")
  cat("external data accessions (GSE292141, Perez CELLxGENE collection,\n")
  cat("LINCS EH7297, GSE70138), and the environment setup notes.\n")
}

cat("\n=================== AUDIT COMPLETE ===================\n")
cat("This report only DESCRIBES the current state. No files were moved,\n")
cat("renamed, or deleted by this script.\n")
