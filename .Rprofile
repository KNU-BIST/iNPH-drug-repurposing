# =============================================================================
# iNPH Drug Repurposing Project - R startup
# =============================================================================

# Activate the project-specific renv environment.
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else if (interactive()) {
  message(
    "\nrenv activation script was not found.\n",
    "Install renv if needed, then run `renv::restore()` from the project root.\n"
  )
}


# -----------------------------------------------------------------------------
# Project-wide options
# -----------------------------------------------------------------------------

n_cores <- parallel::detectCores(logical = TRUE)

if (is.na(n_cores)) {
  n_cores <- 1L
}

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  BioC_mirror = "https://bioconductor.org",
  Ncpus = max(1L, n_cores - 1L),
  download.file.method = "libcurl",
  scipen = 100,
  digits = 4,
  warn = 1
)


# -----------------------------------------------------------------------------
# Friendly interactive startup message
# -----------------------------------------------------------------------------

if (interactive()) {
  cat(
    "\n",
    "iNPH Drug Repurposing Project\n",
    "-----------------------------\n",
    "Project root: ", getwd(), "\n",
    "Environment: renv-managed; run `renv::restore()` after a fresh clone.\n",
    "Preflight:   Rscript run_full_pipeline.R --preflight-only\n",
    "Pipeline:    Rscript run_full_pipeline.R\n\n",
    sep = ""
  )
}
