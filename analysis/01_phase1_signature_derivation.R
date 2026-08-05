# =============================================================================
# Phase 1 — Deriving the type I interferon-activated monocyte signature
# =============================================================================
#
# Project       : iNPH Drug Repurposing
# Author        : Masudul Haque, KNUH Department of Neurosurgery
# Final version : 2026-05-04
#
# Purpose
# -------
# This script is the canonical, reproducible derivation of the disease-state
# query signature used in the iNPH drug repurposing project. Running it from
# top to bottom on a fresh clone of the repository (after `source("setup.R")`)
# reproduces the 150-up / 150-down gene signature in
#     data/processed/inph_ifn_signature.rds
#
# Phase 1 deliverable
# -------------------
# A type I interferon-activated monocyte signature, derived by contrasting
# IFN-responsive monocytes (subcluster 5) against classical resting monocytes
# (subcluster 0) within the same iNPH patients. This signature feeds Phase 2
# (specificity validation) and Phase 3 (LINCS connectivity scoring).
#
# Methodological summary
# ----------------------
# Source dataset    : Duy et al. 2025 PNAS, GSE292141
# Cohort            : 9 iNPH patients (Patient 6 dropped — Unknown MOCA)
# Tissue            : Peripheral blood (PBMC) — switched from CSF (see Pivot 2)
# Cells             : 124,047 raw → 122,724 post-QC → 22,197 monocytes (cluster 0)
#                     → 22,197 monocytes resubclustered → 940 IFN-activated +
#                     4,904 classical-resting monocytes used in DE
# DE design         : ~ patient + state, paired across all 9 patients
# Shrinkage         : apeglm
# Result            : 1,231 genes pass FDR < 0.05; 32 of 35 expected ISGs
#                     recovered in the UP signature
#
# Three scientific pivots from original proposal
# ----------------------------------------------
# Pivot 1 (cohort)    : iNPH-vs-healthy → High-MOCA-vs-Low-MOCA within iNPH
#                       Reason: dataset has no healthy controls.
# Pivot 2 (tissue)    : CSF → PBMC monocytes
#                       Reason: CSF cell yield 18-2,503 (fatal imbalance);
#                       PBMC yield 1,343-4,541 (uniform).
# Pivot 3 (signature) : Generic proinflammatory → type I IFN-activated state
#                       Reason: TNF/IL1β biology absent across all monocyte
#                       subclusters (scores ~0.07); IFN-activated subset
#                       distinct (score 1.28 vs ≤0.67 elsewhere).
#
# These pivots are documented in detail in the project follow-up deck.
#
# Dependencies
# ------------
# All packages are pinned in renv.lock. The starting monocyte population
# (22,197 classical monocytes, PBMC, 9 patients) was originally built via
# the targets pipeline (duy_raw -> duy_extracted -> duy_sce -> duy_clustered
# -> duy_annotated -> duy_monocytes, defined in the now-archived _targets.R,
# using R/data_loading.R and R/clustering.R). That pipeline implements the
# abandoned study design and is intentionally excluded from this repository
# (see Archived/). The saved, verified output of that final target -
# data/processed/duy_monocytes_classical.rds - is loaded directly below
# instead, so this script no longer depends on targets or _targets.R at
# runtime. See README.md for details.
#
# FIX THIS SESSION (fresh-clone reproducibility test)
# -----------------------------------------------------------------------
# This script builds `seu`, the fully subclustered Seurat object (sub_res05
# assignments, patient/moca_group metadata, UMAP embedding), and every
# downstream script in this pipeline - 02_lupus_benchmark.R, 05_robustness.R,
# 06_figures_all.R, 08_numbers_audit.R - loads it back from
# data/processed/seu_monocyte_subclustered.rds. But this script never
# actually wrote that file: Section 10 only ever saved inph_ifn_signature.rds,
# de_ifn_subtype_results.rds, and subcluster_diagnostic.rds. This worked
# invisibly on the original machine only because a copy of that file already
# existed there from an earlier session - the same silent-dependency pattern
# already found and fixed in the targets/_targets.R issue. Fixed by adding
# the missing save in Section 10, right after the other three outputs.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 · Environment
# -----------------------------------------------------------------------------

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

# Sanity check — confirm the pre-built monocyte population is present
stopifnot(file.exists(here::here("data", "processed", "duy_monocytes_classical.rds")))


# -----------------------------------------------------------------------------
# Section 1 · Load the pre-built monocyte population
# -----------------------------------------------------------------------------
# This object was produced (once) by the same processing chain described
# above: 124,047 raw PBMC cells -> 122,724 post-QC -> Harmony integration by
# patient -> clustering at resolution 0.5 -> cluster 0 (the only pure
# classical-monocyte cluster) subset out as the 22,197-cell starting
# population used throughout this project.
#
# NOTE ON REPRODUCIBILITY: the clustering step inside this original chain
# was not seeded. Re-deriving this object from raw data (e.g. via
# R/data_loading.R + R/clustering.R, as done for Supplementary Figure S1)
# will NOT reliably reproduce this exact 22,197-cell population, since
# Harmony/Louvain clustering is stochastic. This saved object is therefore
# the canonical starting point for the locked signature and all downstream
# results - do not regenerate it from scratch as a substitute.


# -----------------------------------------------------------------------------
# Section 2 · Load monocyte SCE and confirm cohort structure
# -----------------------------------------------------------------------------

mono_sce <- readRDS(here::here("data", "processed", "duy_monocytes_classical.rds"))
mono_sce

# Confirm: 22,197 monocytes, 9 patients (Patient 6 already dropped by the
# pipeline), tissue is PBMC (not CSF — see Pivot 2 in project documentation).
stopifnot(ncol(mono_sce) == 22197)
stopifnot(unique(mono_sce$tissue) == "PBMC")

cat("\nCohort structure:\n")
print(table(mono_sce$patient, mono_sce$moca_group))


# -----------------------------------------------------------------------------
# Section 3 · Subcluster the 22,197 monocytes into states
# -----------------------------------------------------------------------------
# Why this step
# -------------
# Pseudobulk DE on bulk monocytes (Low MOCA vs High MOCA contrast) produced
# no meaningful signature — the cognition gradient is not detectable at the
# bulk-monocyte level. This is consistent with Duy's published finding being
# subtype-specific rather than a uniform shift across all monocytes.
# Subclustering preserves subtype structure; pseudobulk would average it away.

seu <- Seurat::CreateSeuratObject(
  counts    = SingleCellExperiment::counts(mono_sce),
  meta.data = as.data.frame(SummarizedExperiment::colData(mono_sce))
)

seu <- Seurat::NormalizeData(seu, verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- Seurat::ScaleData(seu, verbose = FALSE)
seu <- Seurat::RunPCA(seu, npcs = 30, verbose = FALSE)

# Patient-level batch correction. Without this, monocyte subclusters can
# appear patient-driven (cells from the same patient cluster together) rather
# than reflecting genuine biological subtype structure.
seu <- harmony::RunHarmony(seu, group.by.vars = "patient", verbose = FALSE)

seu <- Seurat::FindNeighbors(seu, reduction = "harmony",
                             dims = 1:30, verbose = FALSE)
seu <- Seurat::FindClusters(seu, resolution = 0.5,
                            cluster.name = "sub_res05", verbose = FALSE)
seu <- Seurat::RunUMAP(seu, reduction = "harmony",
                       dims = 1:30, verbose = FALSE)

cat("\nSubcluster sizes at resolution 0.5:\n")
print(table(seu$sub_res05))


# -----------------------------------------------------------------------------
# Section 4 · Identify the IFN-activated subcluster
# -----------------------------------------------------------------------------
# We score each subcluster on six biology panels covering canonical monocyte
# states. The original plan was to find a "proinflammatory" subcluster with
# elevated TNF/IL1β/CCL3/CCL4 expression (Pivot 3 outcome).
#
# Empirical finding: no subcluster has elevated proinflammatory scores
# (uniformly ~0.06-0.09 across all 7 subclusters). However, one subcluster
# shows clearly elevated type I interferon activation:
#   subcluster 5 — IFN score 1.28 vs ≤0.67 elsewhere (~2× elevated)
# This becomes our disease state of interest.

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

inflammatory_subcluster <- "5"   # IFN-activated, n = 940
classical_subcluster    <- "0"   # classical resting reference, n = 4,904

stopifnot(inflammatory_subcluster %in% as.character(seu$sub_res05))
stopifnot(classical_subcluster    %in% as.character(seu$sub_res05))

sce$subtype_label <- as.character(sce$sub_res05)
keep_cells <- sce$subtype_label %in% c(inflammatory_subcluster,
                                       classical_subcluster)
sce_two <- sce[, keep_cells]

sce_two$state <- factor(
  ifelse(sce_two$subtype_label == inflammatory_subcluster,
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

# Drop pseudobulks built from too few cells (project standard relaxed from
# 50 to 30 here because the IFN subcluster is small — 940 cells across 9
# patients).
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

# apeglm shrinks log fold-changes toward zero based on standard error.
# This is the standard correction for small-n contrasts where raw LFCs
# are unstable.
res <- DESeq2::lfcShrink(
  dds,
  coef = "state_ifn_activated_vs_classical",
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

stopifnot(length(intersect(up, down)) == 0)   # disjoint
cat("\nSignature: ", length(up), " UP / ", length(down), " DOWN\n", sep = "")
cat("\nTop 20 UP genes (proinflammatory IFN response):\n")
print(head(up, 20))
cat("\nTop 20 DOWN genes (homeostatic / granulocyte-adjacent):\n")
print(head(down, 20))


# -----------------------------------------------------------------------------
# Section 9 · Biological validation
# -----------------------------------------------------------------------------
# The UP signature should be dominated by canonical type I interferon-
# stimulated genes (ISGs). We check against a panel of 35 well-characterized
# ISGs. Expected outcome: at least 25 of 35 land in the UP signature,
# confirming the signature represents real IFN-activated monocyte biology.

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

# Hard validation criterion: at least 25 of 35 expected ISGs must land.
# If this fails, the signature is not the IFN biology we believe it is and
# downstream phases (LINCS) will produce uninterpretable results.
stopifnot(length(isgs_recovered) >= 25)

cat("\n=== Validation passed. Phase 1 signature is biologically coherent. ===\n")


# -----------------------------------------------------------------------------
# Section 10 · Save final outputs
# -----------------------------------------------------------------------------

dir.create(here::here("data/processed"), showWarnings = FALSE, recursive = TRUE)

saveRDS(inph_signature,
        here::here("data/processed/inph_ifn_signature.rds"))
saveRDS(de_subtype,
        here::here("data/processed/de_ifn_subtype_results.rds"))
saveRDS(diagnostic,
        here::here("data/processed/subcluster_diagnostic.rds"))

# FIX THIS SESSION: this save was missing entirely. seu (fully subclustered:
# sub_res05, patient, moca_group, UMAP embedding) is what 02_lupus_benchmark.R,
# 05_robustness.R, 06_figures_all.R, and 08_numbers_audit.R all load back via
# this exact filename. Without this line, the script silently depended on a
# copy of this file already existing from an earlier session - it never
# actually produced its own most heavily-used output. See header note.
saveRDS(seu,
        here::here("data/processed/seu_monocyte_subclustered.rds"))

cat("\nFinal outputs saved:\n")
cat("  data/processed/inph_ifn_signature.rds          (Phase 1 signature for LINCS query)\n")
cat("  data/processed/de_ifn_subtype_results.rds      (full DE table)\n")
cat("  data/processed/subcluster_diagnostic.rds       (subcluster annotation table)\n")
cat("  data/processed/seu_monocyte_subclustered.rds   (clustered Seurat object - used by\n")
cat("                                                   scripts 02, 05, 06, 08)\n")


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
cat("\n  Next: Phase 2 (specificity validation across CSF compendium)\n")
cat("============================================================\n")

# =============================================================================
# End of Phase 1.
# =============================================================================
