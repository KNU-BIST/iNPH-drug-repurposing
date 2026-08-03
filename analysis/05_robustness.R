# Phase 6A — Positive control
# RA monocytes have a published IFN/JAK signature (GSE93272 or similar)
# We'll use a manually curated RA ISG panel as a quick proxy positive control

# Published RA-elevated ISGs in monocytes (from literature)
ra_up <- c("IFIT1","IFIT2","IFIT3","ISG15","MX1","OAS1","OAS2",
           "RSAD2","IFI44L","HERC5","CXCL10","IFI6","IFITM3",
           "USP18","SIGLEC1","IFI44","OASL","DDX58","IFIH1","TRIM22")

ra_down <- c("S100A8","S100A9","S100A12","CD14","FCGR3A",
             "VCAN","LYZ","CTSS","AIF1","LST1",
             "LILRB2","MS4A6A","CD36","NCF2","FCER1G")
# Load
lincs_genes <- rhdf5::h5read(lincs2_path, "rownames")[, 1]
cat("LINCS2 genes reloaded:", length(lincs_genes), "\n")


# Convert to Entrez
ra_up_entrez <- symbol_to_entrez(ra_up)
ra_down_entrez <- symbol_to_entrez(ra_down)

cat("RA UP genes in LINCS2:", sum(ra_up_entrez %in% lincs_genes), "\n")
cat("RA DOWN genes in LINCS2:", sum(ra_down_entrez %in% lincs_genes), "\n")

qsig_ra <- qSig(
  query = list(upset = ra_up_entrez, downset = ra_down_entrez),
  gess_method = "LINCS",
  refdb = lincs2_path
)

cat("Running positive control LINCS query (~20 min)...\n")
BiocParallel::register(BiocParallel::SerialParam())
ra_result <- gess_lincs(qsig_ra, sortby = "NCS",
                        tau = FALSE, workers = 1)

ra_table <- result(ra_result)

# Check: do known RA drugs appear in reversal hits?
ra_top_thp1 <- ra_table |>
  dplyr::filter(type == "trt_cp", toupper(cell) == "THP1",
                NCS < 0) |>
  dplyr::arrange(NCS) |>
  dplyr::left_join(lincs_meta |>
                     dplyr::select(pert_id, pert_iname) |>
                     dplyr::rename(pert = pert_id), by = "pert") |>
  head(30)

cat("\nRA positive control — top THP1 reversal hits:\n")
print(ra_top_thp1 |> dplyr::select(pert_iname, cell, NCS, WTCS_FDR))

# Known RA drugs in result
ra_drugs <- c("methotrexate", "tofacitinib", "baricitinib",
              "ruxolitinib", "leflunomide", "hydroxychloroquine",
              "dexamethasone", "prednisolone")

cat("\nKnown RA drugs in top reversal hits:\n")
ra_table |>
  dplyr::left_join(lincs_meta |>
                     dplyr::select(pert_id, pert_iname) |>
                     dplyr::rename(pert = pert_id), by = "pert") |>
  dplyr::filter(tolower(pert_iname) %in% ra_drugs, NCS < 0) |>
  dplyr::arrange(NCS) |>
  dplyr::select(pert_iname, cell, NCS, WTCS_FDR) |>
  print()

saveRDS(ra_result, here::here("data/processed/phase6_positive_control_ra.rds"))






################################################################################

# Inspect the pseudobulk object structure
cat("Class:", class(pseudobulk), "\n")
cat("Names/slots:", paste(names(pseudobulk), collapse = ", "), "\n")

# If it's a DESeqDataSet
if (inherits(pseudobulk, "DESeqDataSet") ||
    inherits(pseudobulk, "SummarizedExperiment")) {
  cat("\nColData columns:\n")
  print(colnames(SummarizedExperiment::colData(pseudobulk)))
  cat("\nColData head:\n")
  print(head(SummarizedExperiment::colData(pseudobulk)))
  cat("\nDimensions:", dim(pseudobulk), "\n")
}

# If it's a list (edgeR DGEList or similar)
if (is.list(pseudobulk) && !is.data.frame(pseudobulk)) {
  cat("\nList elements:", paste(names(pseudobulk), collapse = ", "), "\n")
  if ("samples" %in% names(pseudobulk)) {
    cat("samples columns:", paste(colnames(pseudobulk$samples), collapse=", "), "\n")
    print(head(pseudobulk$samples))
  }
  if ("counts" %in% names(pseudobulk)) {
    cat("counts dim:", dim(pseudobulk$counts), "\n")
    print(head(colnames(pseudobulk$counts)))
  }
}





##
library(targets); library(here); library(SingleCellExperiment)
library(SummarizedExperiment); library(Matrix)
library(Seurat); library(harmony); library(DESeq2); library(dplyr)

# --- Reload from pipeline cache and disk ------------------------------------
cat("Loading iNPH monocytes from targets cache...\n")
inph_mono <- targets::tar_read(duy_monocytes)
cat("  inph_mono:", ncol(inph_mono), "cells x", nrow(inph_mono), "genes\n")

cat("Loading Phase 1 signature...\n")
inph_signature <- readRDS(here::here("data/processed/inph_ifn_signature.rds"))
cat("  UP:", length(inph_signature$up),
    "| DOWN:", length(inph_signature$down), "\n")

# --- Rebuild subclustering to recover IFN cell assignments ------------------
cat("\nRebuilding monocyte subclustering (~5-10 min)...\n")
seu <- Seurat::CreateSeuratObject(
  counts    = SingleCellExperiment::counts(inph_mono),
  meta.data = as.data.frame(SummarizedExperiment::colData(inph_mono))
)
seu <- Seurat::NormalizeData(seu,          verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu,   verbose = FALSE)
seu <- Seurat::ScaleData(seu,              verbose = FALSE)
seu <- Seurat::RunPCA(seu, npcs = 30,      verbose = FALSE)
seu <- harmony::RunHarmony(seu, group.by.vars = "patient", verbose = FALSE)
seu <- Seurat::FindNeighbors(seu, reduction = "harmony",
                             dims = 1:30,  verbose = FALSE)
seu <- Seurat::FindClusters(seu, resolution = 0.5,
                            cluster.name = "sub_res05", verbose = FALSE)

# Verify IFN cluster identity
ifn_panel  <- c("IFITM3","IFITM2","ISG15","MX1","OAS1","IRF7","STAT1")
sce_sub    <- Seurat::as.SingleCellExperiment(seu)
lc         <- SingleCellExperiment::logcounts(sce_sub)
ifn_scores <- tapply(
  Matrix::colMeans(lc[intersect(ifn_panel, rownames(lc)), ]),
  seu$sub_res05, mean
) |> sort(decreasing = TRUE)

cat("\nIFN scores per subcluster:\n"); print(round(ifn_scores, 2))
cat("Cluster sizes:\n"); print(table(seu$sub_res05))

ifn_cluster <- names(ifn_scores)[1]
cat("\nIFN-activated cluster:", ifn_cluster,
    "(", sum(seu$sub_res05 == ifn_cluster), "cells )\n")

# Confirm it's cluster 5 with 940 cells as before
stopifnot(sum(seu$sub_res05 == ifn_cluster) == 940)
cat("✓ Subclustering reproduced exactly\n")

# --- Build cell-state metadata ----------------------------------------------
cell_meta <- data.frame(
  cell    = colnames(seu),
  patient = seu$patient,
  state   = ifelse(seu$sub_res05 == ifn_cluster, "IFN", "classical"),
  stringsAsFactors = FALSE
)
cat("\nCell-state counts:\n"); print(table(cell_meta$state))
cat("Patients:", sort(unique(cell_meta$patient)), "\n")








################################################################################
# Phase 6B — Leave-one-out signature stability
# Load Phase 1 pipeline objects needed for LOO

library(targets); library(DESeq2); library(here)

# Fix cell_meta — Phase 1 compared cluster 5 (IFN) vs cluster 0 only
# not all 21,257 classical cells. Match Phase 1 exactly.
cell_meta$state[seu$sub_res05 != ifn_cluster &
                  seu$sub_res05 != "0"] <- NA

cell_meta <- cell_meta[!is.na(cell_meta$state), ]
cat("Cells for LOO (cluster 0 + cluster 5 only):\n")
print(table(cell_meta$state))
# Should be: classical=4904, IFN=940

# --- LOO function -----------------------------------------------------------
counts_mat <- SingleCellExperiment::counts(inph_mono)

run_loo_deseq <- function(drop_patient, cell_meta, counts_mat,
                          full_up, full_dn, top_n = 150) {
  meta_loo <- cell_meta |>
    dplyr::filter(patient != drop_patient)

  meta_loo$group <- paste0("P", meta_loo$patient, "_", meta_loo$state)
  groups <- unique(meta_loo$group)

  # Pseudobulk: sum counts per patient × state
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
                                jaccard(loo_dn, full_dn))), 3)
  )
}

# --- Run all 9 LOO iterations -----------------------------------------------
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

# --- Summary ----------------------------------------------------------------
loo_df <- do.call(rbind, lapply(loo_results, as.data.frame))
rownames(loo_df) <- NULL

cat("\n=== Phase 6 — Leave-One-Out Stability ===\n")
print(loo_df)

cat(sprintf("\nMean Jaccard UP  : %.3f\n", mean(loo_df$jaccard_up)))
cat(sprintf("Mean Jaccard DOWN: %.3f\n", mean(loo_df$jaccard_dn)))
cat(sprintf("Min  Jaccard UP  : %.3f\n", min(loo_df$jaccard_up)))
cat(sprintf("Min  Jaccard DOWN: %.3f\n", min(loo_df$jaccard_dn)))
cat(sprintf("Overall mean     : %.3f\n",
            mean(c(loo_df$jaccard_up, loo_df$jaccard_dn))))
cat(sprintf("\nPass (≥0.70)     : %s\n",
            ifelse(mean(c(loo_df$jaccard_up, loo_df$jaccard_dn)) >= 0.70,
                   "YES ✓", "NO — investigate which patient drives instability")))

saveRDS(loo_df, here::here("data/processed/phase6_loo_stability.rds"))
cat("Saved.\n")



################################################################################
library(ggplot2)

loo_plot_df <- loo_df |>
  dplyr::mutate(
    patient_label  = paste0("P", dropped),
    moca_group     = ifelse(dropped %in% c(1,2,3,4,5), "High MOCA", "Low MOCA"),
    pass           = jaccard_mean >= 0.70
  ) |>
  tidyr::pivot_longer(cols = c(jaccard_up, jaccard_dn),
                      names_to = "set", values_to = "jaccard") |>
  dplyr::mutate(
    set = ifelse(set == "jaccard_up", "UP genes", "DOWN genes")
  )

p_loo <- ggplot(loo_plot_df,
                aes(x = reorder(patient_label, -jaccard_mean),
                    y = jaccard, color = moca_group, shape = set)) +
  geom_hline(yintercept = 0.70, linetype = "dashed",
             color = "#d4624a", linewidth = 0.7) +
  geom_point(size = 4, alpha = 0.9) +
  annotate("text", x = 0.6, y = 0.715, label = "Threshold = 0.70",
           hjust = 0, size = 3.2, color = "#d4624a", fontface = "italic") +
  scale_color_manual(
    values = c("High MOCA" = "#7a8c6a", "Low MOCA" = "#b04a36"),
    name = "MOCA group"
  ) +
  scale_shape_manual(
    values = c("UP genes" = 16, "DOWN genes" = 17),
    name = "Gene set"
  ) +
  scale_y_continuous(limits = c(0.60, 0.85), breaks = seq(0.60, 0.85, 0.05)) +
  labs(
    x        = "Patient removed",
    y        = "Jaccard similarity to full signature",
    title    = "Phase 6 — Leave-one-out signature stability",
    subtitle = paste0(
      "Each point = Jaccard similarity between the LOO signature ",
      "and the full 9-patient signature\n",
      "Overall mean: UP = 0.714 | DOWN = 0.733 | Combined = 0.724"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

print(p_loo)
ggsave(here::here("figures/phase6_loo_stability.pdf"),
       p_loo, width = 7, height = 4.5)
ggsave(here::here("figures/phase6_loo_stability.png"),
       p_loo, width = 7, height = 4.5, dpi = 300)
cat("LOO figure saved.\n")
