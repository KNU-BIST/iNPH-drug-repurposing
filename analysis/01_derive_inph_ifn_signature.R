# =============================================================================
# Phase 1 — iNPH IFN-activated monocyte signature derivation
# Script: 01_derive_inph_ifn_signature.R
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
# Derive the 150-UP/150-DOWN iNPH IFN-activated monocyte signature used in
# the cross-disease contextualization and LINCS2 connectivity analyses.
#
# Input
# -----
# data/processed/duy_monocytes_classical.rds
# Canonical processed object containing 22,197 classical monocytes from
# nine iNPH patients, reconstructed by
# 00_prepare_canonical_monocyte_input.R from public GSE292141 raw counts
# and the version-controlled canonical cell-membership manifest.
#
# Outputs
# -------
# data/processed/inph_ifn_signature.rds
# data/processed/de_ifn_subtype_results.rds
# data/processed/subcluster_diagnostic.rds
# data/processed/seu_monocyte_subclustered.rds
# data/processed/isg_reference_panel.rds
#
# Method summary
# --------------
# The canonical monocyte population is normalized, Harmony-corrected by
# patient, subclustered, and scored using predefined monocyte-state panels.
# Pseudobulk differential expression compares IFN-activated subcluster 5
# with classical resting cluster 0 using DESeq2 with a paired
# ~ patient + state design and apeglm shrinkage.
#
# Reproducibility note
# --------------------
# This script consumes the canonical monocyte object produced by Script 00.
# Script 00 reconstructs the exact locked historical cell membership and order
# from public GSE292141 raw count matrices using the version-controlled manifest.
# It does not independently rediscover the original stochastic PBMC cluster
# selection; that historical membership is intentionally frozen by the manifest.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "Seurat",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scuttle",
  "DESeq2",
  "apeglm",
  "harmony",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run renv::restore() before executing this script."
  )
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(DESeq2)
  library(apeglm)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

PROCESSED_DIR <- here::here("data", "processed")

canonical_mono_path <- file.path(
  PROCESSED_DIR,
  "duy_monocytes_classical.rds"
)

if (
  !file.exists(canonical_mono_path) ||
  is.na(file.info(canonical_mono_path)$size) ||
  file.info(canonical_mono_path)$size <= 0
) {
  stop(
    "Canonical monocyte object is missing or empty: ",
    canonical_mono_path,
    "\nRun analysis/00_prepare_canonical_monocyte_input.R first (or use run_full_pipeline.R)."
  )
}

cat(
  "Canonical input MD5: ",
  unname(tools::md5sum(canonical_mono_path)),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Section 1 · Load the reconstructed canonical monocyte population
# -----------------------------------------------------------------------------
# The 22,197-cell population is produced by Script 00 from the public GSE292141
# PBMC raw count matrices and the version-controlled canonical cell-membership
# manifest.
#
# The manifest preserves the exact historical cell membership and order used in
# the original analysis. This is deliberately different from claiming that the
# original stochastic Harmony/clustering-based monocyte-selection step can be
# independently rediscovered exactly. Script 01 therefore begins from the
# reconstructed, locked analysis population and performs all downstream
# normalization, subclustering, differential expression, and signature
# derivation from that object.
#
# -----------------------------------------------------------------------------

# Section 2 · Load monocyte SCE and confirm cohort structure
# -----------------------------------------------------------------------------

mono_sce <- readRDS(canonical_mono_path)

if (!inherits(mono_sce, "SingleCellExperiment")) {
  stop(
    "Canonical monocyte input must be a SingleCellExperiment; observed class: ",
    paste(class(mono_sce), collapse = ", ")
  )
}

required_coldata <- c(
  "patient",
  "tissue",
  "moca_group"
)

missing_coldata <- setdiff(
  required_coldata,
  colnames(SummarizedExperiment::colData(mono_sce))
)

if (length(missing_coldata) > 0L) {
  stop(
    "Canonical monocyte object is missing required colData field(s): ",
    paste(missing_coldata, collapse = ", ")
  )
}

if (!"counts" %in% SummarizedExperiment::assayNames(mono_sce)) {
  stop("Canonical monocyte object does not contain a 'counts' assay.")
}

stopifnot(ncol(mono_sce) == 22197L)
stopifnot(length(unique(mono_sce$patient)) == 9L)
stopifnot(all(as.character(mono_sce$tissue) == "PBMC"))
stopifnot(!anyDuplicated(colnames(mono_sce)))

cat("\nCohort structure:\n")
print(table(mono_sce$patient, mono_sce$moca_group))


# -----------------------------------------------------------------------------
# Section 3 · Subcluster the 22,197 monocytes into states
# -----------------------------------------------------------------------------
# The canonical monocyte population is subclustered to resolve transcriptional
# states that would be obscured by whole-monocyte aggregation. The finalized
# analysis defines the IFN-activated state at this subcluster level and uses a
# within-patient state contrast for downstream pseudobulk differential
# expression.

seu <- Seurat::CreateSeuratObject(
  counts    = SingleCellExperiment::counts(mono_sce),
  meta.data = as.data.frame(SummarizedExperiment::colData(mono_sce))
)

seu <- Seurat::NormalizeData(seu, verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- Seurat::ScaleData(seu, verbose = FALSE)
seu <- Seurat::RunPCA(
  seu,
  npcs = 30,
  seed.use = 42,
  verbose = FALSE
)

# Patient-level batch correction. Without this, monocyte subclusters can
# appear patient-driven (cells from the same patient cluster together) rather
# than reflecting genuine biological subtype structure.
seu <- harmony::RunHarmony(seu, group.by.vars = "patient", verbose = FALSE)

# The finalized analysis uses Louvain clustering (algorithm 1) at resolution
# 0.5. Seeds and the UMAP implementation are specified explicitly to reduce
# software-default drift across future reruns.
seu <- Seurat::FindNeighbors(
  seu,
  reduction = "harmony",
  dims = seq_len(30),
  verbose = FALSE
)

seu <- Seurat::FindClusters(
  seu,
  resolution = 0.5,
  algorithm = 1,
  random.seed = 0,
  cluster.name = "sub_res05",
  verbose = FALSE
)

seu <- Seurat::RunUMAP(
  seu,
  reduction = "harmony",
  dims = seq_len(30),
  umap.method = "uwot",
  metric = "cosine",
  seed.use = 42,
  verbose = FALSE
)

cat("\nSubcluster sizes at resolution 0.5:\n")
print(table(seu$sub_res05))

observed_clusters <- sort(unique(as.character(seu$sub_res05)))

stopifnot(
  identical(
    observed_clusters,
    as.character(0:6)
  )
)


# -----------------------------------------------------------------------------
# Section 4 · Identify the IFN-activated subcluster
# -----------------------------------------------------------------------------
# Each subcluster is summarized using predefined marker panels representing
# major monocyte transcriptional states. In the finalized analysis, subcluster
# 5 shows the strongest type I interferon program (IFN-panel score 1.28 versus
# no more than 0.67 in the remaining subclusters) and is therefore designated
# the IFN-activated state.

panels <- list(
  classical         = c("CD14", "S100A8", "S100A9", "S100A12", "VCAN", "FCN1"),
  non_classical     = c("FCGR3A", "MS4A7", "CDKN1C", "LST1", "CSF1R", "RHOC"),
  intermediate      = c("HLA-DRA", "HLA-DRB1", "HLA-DPB1", "CD74", "CIITA"),
  proinflammatory   = c("TNF", "IL1B", "IL6", "NLRP3", "CCL3", "CCL4",
                        "CXCL8", "PTGS2", "IL1A", "CCL20"),
  interferon        = c("IFITM3", "IFITM2", "ISG15", "MX1", "OAS1",
                        "IRF7", "STAT1"),
  anti_inflammatory = c("CD163", "MRC1", "MAF", "MAFB", "SELENOP", "STAB1")
)

sce  <- Seurat::as.SingleCellExperiment(seu)
logc <- SingleCellExperiment::logcounts(sce)

score_panel <- function(genes, cluster_col) {
  genes <- intersect(genes, rownames(logc))
  if (length(genes) == 0) {
    return(setNames(rep(NA_real_, length(unique(sce[[cluster_col]]))),
                    levels(factor(sce[[cluster_col]]))))
  }
  per_cell <- Matrix::colMeans(logc[genes, , drop = FALSE])
  tapply(per_cell, sce[[cluster_col]], mean)
}

clusters <- sort(unique(sce$sub_res05))
diagnostic <- data.frame(
  cluster = as.character(clusters),
  n_cells = as.integer(table(sce$sub_res05)[as.character(clusters)])
)
for (panel_name in names(panels)) {
  s <- score_panel(panels[[panel_name]], "sub_res05")
  diagnostic[[panel_name]] <- round(as.numeric(s[as.character(clusters)]), 2)
}

cat("\nSubcluster diagnostic at resolution 0.5:\n")
print(diagnostic)

# Numerical sanity check: proinflammatory absent across all subclusters,
# interferon clearly elevated in subcluster 5.
stopifnot(max(diagnostic$proinflammatory) < 0.15)
stopifnot(diagnostic$interferon[diagnostic$cluster == "5"] >
          1.5 * max(diagnostic$interferon[diagnostic$cluster != "5"]))


# -----------------------------------------------------------------------------
# Section 5 · Define the contrast subclusters
# -----------------------------------------------------------------------------
# Subtype-vs-subtype contrast within the same patients. Each patient
# contributes one pseudobulk per state. Patient is included as a covariate
# in the DE design to absorb between-patient variance.

ifn_subcluster          <- "5"   # IFN-activated, n = 940
classical_subcluster    <- "0"   # classical resting reference, n = 4,904

stopifnot(ifn_subcluster %in% as.character(seu$sub_res05))
stopifnot(classical_subcluster    %in% as.character(seu$sub_res05))

sce$subtype_label <- as.character(sce$sub_res05)
keep_cells <- sce$subtype_label %in% c(ifn_subcluster,
                                       classical_subcluster)
sce_two <- sce[, keep_cells]

sce_two$state <- factor(
  ifelse(sce_two$subtype_label == ifn_subcluster,
         "ifn_activated", "classical"),
  levels = c("classical", "ifn_activated")
)

cat("\nCells per (patient × state):\n")
print(table(sce_two$patient, sce_two$state))


# -----------------------------------------------------------------------------
# Section 6 · Pseudobulk aggregation
# -----------------------------------------------------------------------------

sce_two$pb_group <- paste(sce_two$patient, sce_two$state, sep = "_")

pb <- scuttle::aggregateAcrossCells(
  sce_two, ids = sce_two$pb_group, statistics = "sum"
)

# Require at least 30 cells per patient-by-state pseudobulk. This threshold
# retains all nine paired patient contrasts while avoiding extremely small
# pseudobulk libraries.
ncells  <- table(sce_two$pb_group)
keep_pb <- names(ncells)[ncells >= 30]
drop_pb <- names(ncells)[ncells <  30]

if (length(drop_pb) > 0) {
  cat("\nPseudobulks dropped (n_cells < 30):\n")
  print(ncells[drop_pb])
} else {
  cat("\nAll patient × state pseudobulks retained.\n")
}

pb <- pb[, colnames(pb) %in% keep_pb]
pb$state   <- factor(pb$state,   levels = c("classical", "ifn_activated"))
pb$patient <- factor(pb$patient)

# DESeq2 with `~ patient + state` requires every patient to contribute to
# both states for the design to be identifiable. Drop any patient missing
# one state.
state_counts <- table(pb$patient, pb$state)
both_states  <- rownames(state_counts)[apply(state_counts > 0, 1, all)]
pb <- pb[, pb$patient %in% both_states]
pb$patient <- droplevels(pb$patient)

cat("\nFinal paired pseudobulk design:\n")
print(table(pb$patient, pb$state))
stopifnot(length(both_states) == 9)   # all 9 patients should pair


# -----------------------------------------------------------------------------
# Section 7 · DESeq2 with apeglm shrinkage
# -----------------------------------------------------------------------------

dds <- DESeq2::DESeqDataSet(pb, design = ~ patient + state)
dds <- DESeq2::DESeq(dds, parallel = FALSE)

cat("\nDESeq2 model coefficients:\n")
print(DESeq2::resultsNames(dds))

expected_coef <- "state_ifn_activated_vs_classical"

if (!expected_coef %in% DESeq2::resultsNames(dds)) {
  stop(
    "Expected DESeq2 coefficient not found: ",
    expected_coef
  )
}

# apeglm shrinks log fold-changes toward zero based on standard error.
# This is the standard correction for small-n contrasts where raw LFCs
# are unstable.
res <- DESeq2::lfcShrink(
  dds,
  coef = expected_coef,
  type = "apeglm",
  quiet = TRUE
)

de_subtype <- res |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  tibble::as_tibble() |>
  dplyr::arrange(padj)

cat("\n=== Top 30 DE genes by adjusted p-value ===\n")
print(de_subtype |> dplyr::filter(!is.na(padj)) |> head(30))


# -----------------------------------------------------------------------------
# Section 8 · Extract the signature
# -----------------------------------------------------------------------------
# 150 UP / 150 DOWN, both filtered to FDR < 0.05, ranked by log2FC magnitude
# in their respective direction. Disjoint by construction.

passing <- de_subtype |>
  dplyr::filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange))

cat("\nGenes passing FDR < 0.05:", nrow(passing), "\n")

up <- passing |>
  dplyr::filter(log2FoldChange > 0) |>
  dplyr::arrange(dplyr::desc(log2FoldChange)) |>
  dplyr::slice_head(n = 150) |>
  dplyr::pull(gene)

down <- passing |>
  dplyr::filter(log2FoldChange < 0) |>
  dplyr::arrange(log2FoldChange) |>
  dplyr::slice_head(n = 150) |>
  dplyr::pull(gene)

inph_signature <- list(up = up, down = down)

stopifnot(length(up) == 150L)
stopifnot(length(down) == 150L)
stopifnot(!anyDuplicated(up))
stopifnot(!anyDuplicated(down))
stopifnot(length(intersect(up, down)) == 0L)

cat("\nSignature: ", length(up), " UP / ", length(down), " DOWN\n", sep = "")

cat("\nTop 20 upregulated signature genes:\n")
print(head(up, 20))

cat("\nTop 20 downregulated signature genes:\n")
print(head(down, 20))


# -----------------------------------------------------------------------------
# Section 9 · Biological validation
# -----------------------------------------------------------------------------
# The UP signature is evaluated against a predefined 35-gene type I
# interferon-stimulated gene panel. Recovery of at least 25 panel genes in the
# UP component is used as a prespecified biological-coherence check.

expected_up_ifn <- c(
  "IFI6",  "IFI27", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
  "ISG15", "ISG20", "MX1",   "MX2",   "OAS1",  "OAS2",   "OAS3",   "OASL",
  "RSAD2", "USP18", "DDX60", "IRF7",  "IRF9",  "STAT1",  "STAT2",
  "XAF1",  "BST2",  "PLSCR1","EPSTI1","HERC5", "HERC6",  "PARP9",
  "SAMD9", "SAMD9L","GBP1",  "GBP4",  "GBP5"
)

isgs_recovered <- intersect(expected_up_ifn, inph_signature$up)
cat("\n=== Biology validation ===\n")
cat("Expected ISGs recovered in UP signature:",
    length(isgs_recovered), "of", length(expected_up_ifn), "\n")
print(isgs_recovered)

# Prespecified validation criterion: at least 25 of 35 expected ISGs must be
# recovered in the UP signature. Failure triggers a hard stop before downstream
# connectivity analysis.
stopifnot(length(isgs_recovered) >= 25)

cat("\n=== Validation passed. Phase 1 signature is biologically coherent. ===\n")


# -----------------------------------------------------------------------------
# Section 10 · Save final outputs
# -----------------------------------------------------------------------------

dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(inph_signature,
        file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
saveRDS(de_subtype,
        file.path(PROCESSED_DIR, "de_ifn_subtype_results.rds"))
saveRDS(diagnostic,
        file.path(PROCESSED_DIR, "subcluster_diagnostic.rds"))
saveRDS(seu,
        file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"))

# The 35-gene ISG validation panel used in Section 9 is saved so
# 06_audit_manuscript_results.R can read the same canonical list directly.
saveRDS(expected_up_ifn,
        file.path(PROCESSED_DIR, "isg_reference_panel.rds"))

cat("\nFinal outputs saved:\n")
cat("  data/processed/inph_ifn_signature.rds          (Phase 1 signature for LINCS query)\n")
cat("  data/processed/de_ifn_subtype_results.rds      (full DE table)\n")
cat("  data/processed/subcluster_diagnostic.rds       (subcluster annotation table)\n")
cat("  data/processed/seu_monocyte_subclustered.rds   (clustered Seurat object - used by\n")
cat("                                                   scripts 02, 04, 05, 06)\n")


# -----------------------------------------------------------------------------
# Section 11 · Phase 1 closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PHASE 1 COMPLETE\n")
cat("============================================================\n")
cat("  Cells in contrast        :", ncol(sce_two), "\n")
cat("  Patients (paired design) :", length(both_states), "\n")
cat("  Genes tested             :", nrow(de_subtype), "\n")
cat("  Genes pass FDR < 0.05    :", nrow(passing), "\n")
cat("  ISGs recovered           :", length(isgs_recovered), "/",
                                    length(expected_up_ifn), "\n")
cat("  UP / DOWN signature size :", length(up), "/", length(down), "\n")
cat("\n  Next: Phase 2 (cross-disease contextualization)\n")
cat("============================================================\n")

# =============================================================================
# End of Phase 1.
# =============================================================================
