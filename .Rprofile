# Activate renv if it has been initialised; otherwise prompt to run setup.R.
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else if (interactive()) {
  message(
    "\nrenv has not been initialised yet.\n",
    "Run `source(\"setup.R\")` once to install all dependencies.\n"
  )
}

# Project-wide options
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  BioC_mirror = "https://bioconductor.org",
  Ncpus = max(1L, parallel::detectCores() - 1L),
  download.file.method = "libcurl",
  scipen = 100,
  digits = 4,
  warn = 1
)

# Friendly startup message
if (interactive()) {
  cat("\n",
      "iNPH Drug Repurposing Project\n",
      "-----------------------------\n",
      "Project root: ", getwd(), "\n",
      "Run setup.R once on first clone to install dependencies.\n",
      "Pipeline: targets::tar_make() — see README for full workflow.\n\n",
      sep = "")
}
