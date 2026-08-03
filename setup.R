# =============================================================================
# setup.R — one-time bootstrap for the iNPH drug repurposing project
# =============================================================================
# Run this ONCE on first clone. It will:
#   1. Install renv if missing
#   2. Initialise renv and pin packages from CRAN + Bioconductor + GitHub
#   3. Create renv.lock so collaborators get the same versions
#
# Usage:  source("setup.R")
# =============================================================================

message("--- iNPH Drug Repurposing: bootstrap ---")

# 1. Ensure renv is available -------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# 2. Initialise renv if not already initialised -------------------------------
if (!file.exists("renv.lock")) {
  message("Initialising renv (this creates renv/ and renv.lock)...")
  renv::init(bare = TRUE, restart = FALSE)
}

# 3. Define the package set ---------------------------------------------------
# Grouped by purpose so it's easy to add/remove later.

cran_pkgs <- c(
  # Core workflow
  "renv", "here", "targets", "tarchetypes", "config",
  # Data wrangling and I/O
  "tidyverse", "data.table", "fs", "glue", "janitor", "readr", "vroom",
  # Logging and reporting
  "logger", "cli", "quarto", "knitr", "rmarkdown",
  # Plotting
  "ggplot2", "patchwork", "ggrepel", "scales", "viridis", "ggsci",
  "ggpubr", "ggridges", "ggbeeswarm",
  # Statistics
  "broom", "MASS",
  # Networking and APIs
  "httr2", "jsonlite", "rvest",
  # Testing and style
  "testthat", "lintr", "styler", "covr",
  # Drug / chemistry helpers
  "webchem"
)

bioc_pkgs <- c(
  # scRNA-seq core
  "SingleCellExperiment", "SummarizedExperiment", "scran", "scater",
  "DropletUtils", "zellkonverter",
  # DE and pathway
  "DESeq2", "edgeR", "limma", "fgsea", "GSEABase", "msigdbr",
  # Signature scoring
  "AUCell", "GSVA", "singscore",
  # GEO / Bioconductor data fetching
  "GEOquery", "rhdf5", "HDF5Array",
  # CMap / LINCS
  "signatureSearch", "signatureSearchData",
  # Visualisation
  "ComplexHeatmap", "EnhancedVolcano",
  # Annotation
  "org.Hs.eg.db", "AnnotationDbi", "biomaRt"
)

github_pkgs <- c(
  # cmapR for LINCS GCTX file parsing
  "cmap/cmapR",
  # Seurat is on CRAN now but Azimuth is GitHub-only
  "satijalab/azimuth"
)

# 4. Install ------------------------------------------------------------------
message("Installing CRAN packages...")
renv::install(cran_pkgs, prompt = FALSE)

message("Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager", prompt = FALSE)
}
renv::install(paste0("bioc::", bioc_pkgs), prompt = FALSE)

message("Installing Seurat (heavy — please be patient)...")
renv::install("Seurat", prompt = FALSE)

message("Installing GitHub packages...")
renv::install(github_pkgs, prompt = FALSE)

# If installation fails
install.packages("pak")
pak::pkg_install(c("bioc::bluster", "bioc::scran"))
renv::snapshot(prompt = FALSE)


# 5. Snapshot -----------------------------------------------------------------
message("Writing renv.lock...")
renv::snapshot(prompt = FALSE)

# 6. Sanity check -------------------------------------------------------------
message("\n--- Setup complete. Sanity checks: ---")
cat("R version:        ", R.version.string, "\n")
cat("renv version:     ", as.character(packageVersion("renv")), "\n")
cat("targets version:  ", as.character(packageVersion("targets")), "\n")
cat("Seurat version:   ", as.character(packageVersion("Seurat")), "\n")
cat("DESeq2 version:   ", as.character(packageVersion("DESeq2")), "\n")
cat("\nNext step: targets::tar_make()\n")
