# =============================================================================
# Phase 6 — Leave-one-out signature stability + gene-level persistence
# =============================================================================
#
# WHAT CHANGED THIS SESSION (fresh-clone reproducibility test)
# --------------------------------------------------------------------------
# This file previously contained THREE unrelated things pasted together:
#
# 1. An RA "positive control" LINCS query. The exact circular 19+14 gene
#    query traced and retired earlier in this project - all 19 UP genes
#    overlapped the iNPH signature itself, never traceable to any real
#    published RA paper. Dr. Nila's explicit instruction: "should not be
#    mentioned or retained." REMOVED ENTIRELY. Also used tau=FALSE
#    (inconsistent with tau=TRUE everywhere else) and referenced an
#    undefined `lincs_meta` object - the same undefined-object bug found in
#    the retired CLUEio script.
#
# 2. A dead diagnostic fragment inspecting an undefined `pseudobulk` object.
#    Pure scratch code that would error if ever reached. REMOVED.
#
# 3. The actual robustness analysis - kept, but fixed. It depended on
#    targets::tar_read(duy_monocytes) (same bug as script 01) AND re-ran
#    Harmony/Louvain clustering from scratch to recover IFN cluster
#    assignments, unseeded - the same stochasticity caveat documented for
#    Supplementary Figure S1. FIXED by loading the canonical, already-
#    clustered seu_monocyte_subclustered.rds directly instead - removes the
#    targets dependency AND the reclustering stochasticity risk at once.
#
# Also removed: the duplicate LOO figure-generation block. 06_figures_all.R's
# Fig4A already builds this figure, more completely, from the same
# phase6_loo_stability.rds this script saves.
#
# GENE-LEVEL PERSISTENCE (phase6_loo_gene_persistence.csv) - RESOLVED
# --------------------------------------------------------------------------
# The code that originally produced this file was missing entirely from
# this script, despite 06_figures_all.R (Fig4B) and 08_numbers_audit.R both
# reading it. Rebuilt from run_loo_deseq()'s per-iteration top-150 gene
# list, which the function already computed for the Jaccard calculation but
# never returned.
#
# Two candidate definitions of "persistence" were tested and compared
# earlier this session:
#   Option 1 (adopted) - a gene counts as "entering the signature" in an
#     iteration if it lands in that iteration's own top-150-by-log2FC list,
#     the SAME ranking rule 01_phase1_signature_derivation.R uses to define
#     the actual 150-gene signature. No secondary filter.
#   Option 2 (rejected) - same rank-150 selection, plus an additional
#     |log2FC| > 1 requirement borrowed from Figure 1D's volcano-plot
#     coloring threshold. Tested directly: it barely changed the "tracked"
#     count but pushed "stable in 9/9" further from internal consistency,
#     not closer - evidence against this definition, not just a style
#     preference.
# Option 1 was adopted specifically because it is the SAME rule already
# defining the signature everywhere else in the pipeline, not an
# independently invented threshold. Confirmed by manuscript update: Figure
# 4's legend and Results text now state 137 genes tracked / 107 stable in
# 9/9, matching Option 1's output exactly, reproduced identically across
# every run this session.
#
# This version saves directly to the canonical filename
# (phase6_loo_gene_persistence.csv) that 06 and 08 expect - no manual
# rename step required after running this script.
#
# =============================================================================

library(here); library(dplyr); library(tidyr); library(readr)
library(Seurat); library(SingleCellExperiment); library(SummarizedExperiment)
library(Matrix); library(DESeq2)

PROCESSED_DIR <- here::here("data", "processed")

# -----------------------------------------------------------------------------
# Section 1 . Load the canonical, already-clustered monocyte object
# -----------------------------------------------------------------------------
cat("Loading canonical clustered monocyte object...\n")
seu <- readRDS(file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"))
cat("  Cells:", ncol(seu), "| Patients:", length(unique(seu$patient)), "\n")

cat("Loading Phase 1 signature...\n")
inph_signature <- readRDS(file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
cat("  UP:", length(inph_signature$up), "| DOWN:", length(inph_signature$down), "\n")

ifn_cluster <- "5"
stopifnot(sum(seu$sub_res05 == ifn_cluster) == 940)
cat("Cluster 5 (IFN-activated):", sum(seu$sub_res05 == ifn_cluster), "cells - confirmed\n")

# -----------------------------------------------------------------------------
# Section 2 . Build cell-state metadata (cluster 0 + cluster 5 only)
# -----------------------------------------------------------------------------
cell_meta <- data.frame(
  cell    = colnames(seu),
  patient = seu$patient,
  state   = dplyr::case_when(
    seu$sub_res05 == ifn_cluster ~ "IFN",
    seu$sub_res05 == "0"         ~ "classical",
    TRUE                          ~ NA_character_
  ),
  stringsAsFactors = FALSE
)
cell_meta <- cell_meta[!is.na(cell_meta$state), ]
cat("\nCells for LOO (cluster 0 + cluster 5 only):\n")
print(table(cell_meta$state))
# Expect: classical = 4904, IFN = 940

counts_mat <- Seurat::GetAssayData(seu, layer = "counts")

# -----------------------------------------------------------------------------
# Section 3 . Leave-one-out function - returns the top-150 gene vectors it
# computes internally, used for both the Jaccard calculation and the
# gene-level persistence tracking in Section 6
# -----------------------------------------------------------------------------
run_loo_deseq <- function(drop_patient, cell_meta, counts_mat,
                          full_up, full_dn, top_n = 150) {
  meta_loo <- cell_meta |>
    dplyr::filter(patient != drop_patient)

  meta_loo$group <- paste0("P", meta_loo$patient, "_", meta_loo$state)
  groups <- unique(meta_loo$group)

  pb_mat <- sapply(groups, function(g) {
    cells <- meta_loo$cell[meta_loo$group == g]
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })

  pb_col <- data.frame(
    group   = groups,
    patient = factor(gsub("P(\\d+)_.*", "\\1", groups)),
    state   = factor(gsub(".*_(IFN|classical)", "\\1", groups),
                     levels = c("classical", "IFN")),
    row.names = groups,
    stringsAsFactors = FALSE
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(pb_mat),
    colData   = pb_col,
    design    = ~ patient + state
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::lfcShrink(dds, coef = "state_IFN_vs_classical",
                           type = "apeglm", quiet = TRUE)

  res_df <- as.data.frame(res) |>
    dplyr::filter(!is.na(padj), padj < 0.05)

  # Top-150-by-log2FC, same ranking rule as the full-cohort signature
  # derivation in 01_phase1_signature_derivation.R. This IS the definition
  # of "entering the signature" used for gene-level persistence in
  # Section 6 below - no secondary filter.
  loo_up <- rownames(res_df |>
                       dplyr::arrange(dplyr::desc(log2FoldChange)) |> head(top_n))
  loo_dn <- rownames(res_df |>
                       dplyr::arrange(log2FoldChange) |> head(top_n))

  jaccard <- function(a, b)
    length(intersect(a, b)) / length(union(a, b))

  list(
    dropped      = drop_patient,
    n_sig        = nrow(res_df),
    jaccard_up   = round(jaccard(loo_up, full_up), 3),
    jaccard_dn   = round(jaccard(loo_dn, full_dn), 3),
    jaccard_mean = round(mean(c(jaccard(loo_up, full_up),
                                jaccard(loo_dn, full_dn))), 3),
    loo_up       = loo_up
  )
}

# -----------------------------------------------------------------------------
# Section 4 . Run all 9 LOO iterations
# -----------------------------------------------------------------------------
all_patients <- sort(unique(cell_meta$patient))
cat("\nStarting LOO —", length(all_patients),
    "iterations (~4 min each, ~35 min total)\n\n")

loo_results <- list()
for (pt in all_patients) {
  cat(sprintf("Dropping patient %2d ... ", pt))
  t0 <- proc.time()
  loo_results[[as.character(pt)]] <- run_loo_deseq(
    drop_patient = pt,
    cell_meta    = cell_meta,
    counts_mat   = counts_mat,
    full_up      = inph_signature$up,
    full_dn      = inph_signature$down
  )
  elapsed <- round((proc.time() - t0)["elapsed"])
  cat(sprintf("done (%ds) | sig genes: %d | Jaccard mean: %.3f\n",
              elapsed,
              loo_results[[as.character(pt)]]$n_sig,
              loo_results[[as.character(pt)]]$jaccard_mean))
}

# -----------------------------------------------------------------------------
# Section 5 . Summary + save (Jaccard-only table, unchanged shape throughout)
# -----------------------------------------------------------------------------
loo_df <- do.call(rbind, lapply(loo_results, function(x) {
  as.data.frame(x[c("dropped", "n_sig", "jaccard_up", "jaccard_dn", "jaccard_mean")])
}))
rownames(loo_df) <- NULL

cat("\n=== Phase 6 — Leave-One-Out Stability ===\n")
print(loo_df)

cat(sprintf("\nMean Jaccard UP  : %.3f\n", mean(loo_df$jaccard_up)))
cat(sprintf("Mean Jaccard DOWN: %.3f\n", mean(loo_df$jaccard_dn)))
cat(sprintf("Min  Jaccard UP  : %.3f\n", min(loo_df$jaccard_up)))
cat(sprintf("Min  Jaccard DOWN: %.3f\n", min(loo_df$jaccard_dn)))
cat(sprintf("Overall mean     : %.3f\n",
            mean(c(loo_df$jaccard_up, loo_df$jaccard_dn))))
cat(sprintf("\nPass (>=0.70)    : %s\n",
            ifelse(mean(c(loo_df$jaccard_up, loo_df$jaccard_dn)) >= 0.70,
                   "YES", "NO - investigate which patient drives instability")))

saveRDS(loo_df, file.path(PROCESSED_DIR, "phase6_loo_stability.rds"))
cat("\nSaved: phase6_loo_stability.rds\n")

# -----------------------------------------------------------------------------
# Section 6 . Gene-level persistence - saved directly under the canonical
# filename that 06_figures_all.R and 08_numbers_audit.R expect. No manual
# rename needed after this script runs.
# -----------------------------------------------------------------------------
cat("\n=== Gene-level persistence ===\n")

all_loo_up <- unlist(lapply(loo_results, function(x) x$loo_up))
gene_counts_up <- table(all_loo_up)

persist_df <- data.frame(
  gene    = names(gene_counts_up),
  n_iters = as.integer(gene_counts_up),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(in_sig = gene %in% inph_signature$up) |>
  dplyr::arrange(dplyr::desc(n_iters))

n_tracked <- sum(persist_df$in_sig)
n_stable  <- sum(persist_df$n_iters == 9 & persist_df$in_sig)

cat("Genes tracked, in original UP signature, appearing in >=1 LOO iteration\n")
cat("  (manuscript: 137):", n_tracked, "\n")
cat("Of those, stable in all 9 iterations\n")
cat("  (manuscript: 107):", n_stable, "\n")

readr::write_csv(persist_df, file.path(PROCESSED_DIR, "phase6_loo_gene_persistence.csv"))
cat("\nSaved: phase6_loo_gene_persistence.csv (canonical filename, ", nrow(persist_df),
    " rows total, ", n_tracked, " with in_sig=TRUE)\n", sep = "")

cat("\n=== Phase 6 robustness (LOO stability + gene persistence) complete ===\n")
