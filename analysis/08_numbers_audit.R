# =============================================================================
# MANUSCRIPT NUMBERS AUDIT — recompute every key statistic from saved objects
# Paste the full output back for cross-checking against the manuscript text.
# =============================================================================
#
# CHANGES FROM Final_audit.R (this session):
#
#   Sections 1-5, 9, 15: UNCHANGED. Untouched by the Figure 2 / Phase 3 rebuild;
#   master file already records these as verified.
#
#   Sections 6-8 (lupus benchmark): FULLY REPLACED. Old objects
#   (phase2_donorlevel_benchmark.rds, phase2_threshold_sensitivity.rds,
#   phase2_unselected_comparison.rds) were built on scaled data with a
#   circular ISG selector - the "comparable to lupus" claim they supported is
#   invalidated. Now reads phase2_benchmark_v4.rds (donor-balanced threshold,
#   selector-excluded, raw counts), approved by Dr. Nila 2026-07-28:
#   population-level comparison (settled, sub-SLE), cluster-5/ISG-high overlap
#   (new - validates subcluster identity), and IFN-high burden (exploratory).
#
#   Section 10 (LINCS screen): BRD-only filter added to the THP1 count, to
#   match 03_lincs_connectivity.R exactly. Previously this check could in
#   principle diverge from the manuscript's stated methodology (which
#   explicitly excludes 45 non-drug genetic-perturbation entries); it happened
#   to still equal 597 without the filter, but that was not guaranteed.
#
#   Sections 11-13 (compound ID, MOA, PDE class): FULLY REPLACED. Old version
#   read supp_table3_thp1_ranked_moa.csv, sourced from a pipeline with real
#   problems: a hardcoded F:-drive path (now disconnected), an undefined
#   `lincs_meta` object, a tau/ranking mismatch with the corrected script, and
#   a webchem/PubChem dependency. Book1.xlsx (the manual annotation file that
#   pipeline used) was recovered and inspected directly: only 5 rows, only 2
#   with real target/MOA content (rosuvastatin, GDC-0068) - both already known
#   independently. Section 12 now checks LINCS2's own embedded MOAss/t_gn_sym
#   columns instead, diagnostically (coverage here has not been verified
#   before this run - do not assume it will match old figures). Section 13's
#   PDE4 class check explicitly flags that oglemilast's classification is
#   expert-asserted, not database-annotated (confirmed via Book1.xlsx: blank).
#
#   Section 14 (JAK2 inhibitor): original check kept as-is, baricitinib check
#   added alongside it (FDA/EMA-approved JAK1/2 inhibitor, found among the
#   THP1 significant reversers this session - a second, stronger convergent
#   JAK-pathway hit).
#
# =============================================================================

suppressMessages({library(dplyr); library(Seurat)})
PROCESSED_DIR <- here::here("data", "processed")
cat("=================== iNPH MANUSCRIPT NUMBERS AUDIT ===================\n\n")

## ---- 1. COHORT & CELL COUNTS ----
cat("---- 1. COHORT & CELL COUNTS ----\n")
seu <- readRDS(file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"))
cat("Total classical monocytes (manuscript: 22,197):", ncol(seu), "\n")
tab <- table(seu$sub_res05)
cat("Subcluster sizes:\n"); print(tab)
cat("IFN subcluster (cluster 5) size (manuscript: 940):", sum(seu$sub_res05=="5"), "\n")
cat("Cluster 0 size (manuscript: 4,904):", sum(seu$sub_res05=="0"), "\n")
cat("IFN as % of monocytes (manuscript: 4.2%):",
    round(100*sum(seu$sub_res05=="5")/ncol(seu),2), "%\n")
cat("Resting majority % (manuscript: 95.8%):",
    round(100*(1-sum(seu$sub_res05=="5")/ncol(seu)),2), "%\n")
cat("N patients (manuscript: 9):", length(unique(seu$patient)), "\n")
cat("Cells per IFN cluster per patient (manuscript range 50-179):\n")
print(range(table(seu$patient[seu$sub_res05=="5"])))
cat("\n")

## ---- 2. IFN PANEL SCORE PER SUBCLUSTER (manuscript: 1.28 vs <=0.67) ----
cat("---- 2. IFN PANEL SCORES ----\n")
isg_panel <- c("IFITM3","IFITM2","ISG15","MX1","OAS1","IRF7","STAT1")
DefaultAssay(seu) <- "RNA"
expr <- GetAssayData(seu, layer="data")
panel_in <- intersect(rownames(expr), isg_panel)
panel_score <- Matrix::colMeans(expr[panel_in,,drop=FALSE])
ps_by_clust <- tapply(panel_score, seu$sub_res05, mean)
cat("Mean ISG panel score per subcluster:\n"); print(round(ps_by_clust,3))
cat("Cluster 5 (manuscript 1.28):", round(ps_by_clust["5"],3), "\n")
cat("Max of others (manuscript <=0.67):", round(max(ps_by_clust[names(ps_by_clust)!="5"]),3), "\n")
cat("Cluster 0 panel score (manuscript 0.47 in methods):", round(ps_by_clust["0"],3), "\n\n")

## ---- 3. MOCA SUBGROUP WILCOXON (manuscript: median 7.2% vs 3.0%, p=0.556) ----
cat("---- 3. MOCA SUBGROUP COMPARISON ----\n")
ifn_pct <- sapply(unique(seu$patient), function(p){
  idx <- seu$patient==p
  100*sum(idx & seu$sub_res05=="5")/sum(idx)
})
names(ifn_pct) <- unique(seu$patient)
high_pat <- c(1,2,3,4,5); low_pat <- c(7,8,9,10)
hi <- ifn_pct[as.character(high_pat)]; lo <- ifn_pct[as.character(low_pat)]
cat("Per-patient IFN%:\n"); print(round(ifn_pct,2))
cat("High MOCA median (manuscript 3.0%):", round(median(hi,na.rm=TRUE),2), "\n")
cat("Low MOCA median (manuscript 7.2%):", round(median(lo,na.rm=TRUE),2), "\n")
wt <- wilcox.test(lo, hi)
cat("Wilcoxon p (manuscript 0.556):", round(wt$p.value,3), "| W =", wt$statistic, "\n")
cat("Highest single value & its group (manuscript: highest is High MOCA):\n")
print(ifn_pct[which.max(ifn_pct)])
cat("\n")

## ---- 4. DIFFERENTIAL EXPRESSION (manuscript: 1,231 genes; PARP14 padj 1.58e-87) ----
cat("---- 4. DIFFERENTIAL EXPRESSION ----\n")
de <- readRDS(file.path(PROCESSED_DIR, "de_ifn_subtype_results.rds"))
de_df <- as.data.frame(de)
if (!"gene" %in% colnames(de_df)) de_df$gene <- rownames(de_df)
sig <- de_df[!is.na(de_df$padj) & de_df$padj < 0.05, ]
cat("Significant genes at FDR<0.05 (manuscript 1,231):", nrow(sig), "\n")
parp14 <- de_df[de_df$gene=="PARP14",]
cat("PARP14 padj (manuscript 1.58e-87):", format(parp14$padj, scientific=TRUE, digits=3), "\n")
cat("PARP14 raw p:", format(parp14$pvalue, scientific=TRUE, digits=3), "\n\n")

## ---- 5. CANONICAL ISG RECOVERY (reconstructed panel - see note) ----
cat("---- 5. ISG RECOVERY ----\n")
sig_obj <- readRDS(file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
cat("Signature UP genes (manuscript 150):", length(sig_obj$up), "\n")
cat("Signature DOWN genes (manuscript 150):", length(sig_obj$down), "\n\n")

cat("NOTE: the original 35-gene panel behind the manuscript's '32 of 35\n")
cat("(91.4%)' claim could not be located anywhere. The panel below is a\n")
cat("FRESH reconstruction from type I IFN literature categories (classic\n")
cat("antiviral effector ISGs, IFN signaling components, genes standard in\n")
cat("lupus/autoimmune IFN-signature panels), built WITHOUT looking at the\n")
cat("signature gene list first, so it cannot be tuned toward a favorable\n")
cat("overlap. This changes the specific gene list underlying a Methods\n")
cat("claim - confirm with Dr. Nila before this goes back into the text.\n\n")

isg_reference_35 <- c(
  # Classic antiviral effector ISGs
  "ISG15","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3",
  "MX1","MX2","OAS1","OAS2","OAS3","OASL","RSAD2","EIF2AK2",
  "DDX58","IFIH1","BST2",
  # IFN signaling components
  "IRF7","IRF9","STAT1","STAT2",
  # Standard in lupus/autoimmune IFN-signature panels + other canonical ISGs
  "IFI6","IFI27","IFI44","IFI44L","LY6E","SIGLEC1","USP18",
  "HERC5","HERC6","CMPK2","XAF1","CXCL10","PARP14"
)
stopifnot(length(isg_reference_35) == 35)

recovered <- intersect(isg_reference_35, sig_obj$up)
cat("Reconstructed 35-gene panel - genes recovered in signature UP set:\n")
cat(sprintf("  %d of 35 (%.1f%%)\n", length(recovered), 100*length(recovered)/35))
cat("Recovered:\n"); print(sort(recovered))
cat("NOT recovered:\n"); print(sort(setdiff(isg_reference_35, sig_obj$up)))
cat("\n")

## ---- 6. LUPUS BENCHMARK - POPULATION LEVEL (v4, donor-balanced) ----
cat("---- 6. LUPUS BENCHMARK: POPULATION-LEVEL (settled, Dr. Nila-approved) ----\n")
bench <- readRDS(file.path(PROCESSED_DIR, "phase2_benchmark_v4.rds"))
A <- bench$engine_A_primary    # AUCell
B <- bench$engine_B_sensitivity # mean(UP)-mean(DOWN)

cat("AUCell (primary engine):\n")
cat(sprintf("  medians  iNPH-c5 / SLE-hi / HC-hi (expect 0.0307 / 0.0702 / 0.0252): %.4f / %.4f / %.4f\n",
            A$med_pop["inph_c5"], A$med_pop["sle_hi"], A$med_pop["hc_hi"]))
cat(sprintf("  iNPH-c5 vs SLE-hi p (expect 6.42e-06): %.3g\n", A$p_c5_sle))
cat(sprintf("  iNPH-c5 vs HC-hi  p (expect 0.213):    %.3g\n", A$p_c5_hc))

cat("\nmean(UP)-mean(DOWN) (sensitivity engine):\n")
cat(sprintf("  medians  iNPH-c5 / SLE-hi / HC-hi (expect 0.0725 / 0.1886 / 0.0599): %.4f / %.4f / %.4f\n",
            B$med_pop["inph_c5"], B$med_pop["sle_hi"], B$med_pop["hc_hi"]))
cat(sprintf("  iNPH-c5 vs SLE-hi p (expect 1.48e-06): %.3g\n", B$p_c5_sle))
cat(sprintf("  iNPH-c5 vs HC-hi  p (expect 0.356):    %.3g\n", B$p_c5_hc))
cat("\nBoth engines agree: iNPH is genuinely IFN-positive but significantly\n")
cat("below interferon-high SLE. Confirmed across 4 independent runs\n")
cat("(150-cap, uncapped, v3, v4) - treat as settled, not re-litigated.\n\n")

## ---- 7. CLUSTER-5 / ISG-HIGH OVERLAP (subcluster identity validation) ----
cat("---- 7. CLUSTER-5 / ISG-HIGH OVERLAP ----\n")
cat(sprintf("AUCell: %.1f%% of cluster 5 clears the shared IFN-high threshold (expect 95.5%%)\n",
            A$overlap["pct_c5_is_ifnhi"]))
cat(sprintf("mean-diff: %.1f%% of cluster 5 clears the shared IFN-high threshold (expect 99.5%%)\n",
            B$overlap["pct_c5_is_ifnhi"]))
cat(sprintf("Donor-balanced shared threshold - AUCell: %.4f | mean-diff: %.4f\n",
            A$thr_shared, B$thr_shared))
cat("Interpretation: unsupervised clustering (cluster 5) and an independent\n")
cat("marker-panel definition land on almost the same cells - strong\n")
cat("orthogonal validation of the subcluster identity itself.\n\n")

## ---- 8. IFN-HIGH BURDEN (exploratory/supportive, per Dr. Nila) ----
cat("---- 8. IFN-HIGH MONOCYTE BURDEN ----\n")
cat("AUCell:\n")
cat(sprintf("  median burden HC / SLE / iNPH (expect 25.0%% / 91.5%% / 13.5%%): %.1f%% / %.1f%% / %.1f%%\n",
            A$burden_med["hc"], A$burden_med["sle"], A$burden_med["inph"]))
cat(sprintf("  iNPH vs SLE burden p (expect 1.29e-05): %.3g\n", A$p_burden_sle))
cat(sprintf("  iNPH vs HC burden p  (expect 0.0128):   %.3g\n", A$p_burden_hc))
cat("mean-diff:\n")
cat(sprintf("  median burden HC / SLE / iNPH (expect 25.0%% / 91.8%% / 30.1%%): %.1f%% / %.1f%% / %.1f%%\n",
            B$burden_med["hc"], B$burden_med["sle"], B$burden_med["inph"]))
cat(sprintf("  iNPH vs SLE burden p (expect 0.00098): %.3g\n", B$p_burden_sle))
cat(sprintf("  iNPH vs HC burden p  (expect 0.491):   %.3g\n", B$p_burden_hc))
cat("NOTE per Dr. Nila: iNPH-vs-SLE burden is robust across both engines -\n")
cat("safe as the main benchmark conclusion. iNPH-vs-HC burden is scoring-\n")
cat("method-dependent (AUCell: significantly below; mean-diff: not\n")
cat("significant) - report as such, do NOT state a true reduction below\n")
cat("healthy baseline.\n\n")

## ---- 9. BULK VS SUBCLUSTER SCORING (manuscript: bulk mean -0.101) ----
cat("---- 9. BULK VS SUBCLUSTER SCORING ----\n")
up_in <- intersect(rownames(expr), sig_obj$up)
down_in <- intersect(rownames(expr), sig_obj$down)
score <- Matrix::colMeans(expr[up_in,,drop=FALSE]) - Matrix::colMeans(expr[down_in,,drop=FALSE])
cat("UP genes present:", length(up_in), "| DOWN present:", length(down_in), "\n")
cat("ALL monocytes mean (manuscript -0.101):", round(mean(score),4),
    "| median:", round(median(score),4), "\n")
cat("IFN-activated mean (manuscript +0.129):", round(mean(score[seu$sub_res05=="5"]),4),
    "| median:", round(median(score[seu$sub_res05=="5"]),4), "\n\n")

## ---- 10. LINCS2 SCREEN (BRD-only filter now matches 03_lincs_connectivity.R exactly) ----
cat("---- 10. LINCS2 SCREEN ----\n")
res <- readRDS(file.path(PROCESSED_DIR, "lincs2_gess_result.rds")) |> signatureSearch::result()
cat("Total signatures (manuscript 136,460):", nrow(res), "\n")
cat("Unique cell lines (manuscript 58):", length(unique(res$cell)), "\n")

all_cpds <- res |> dplyr::filter(type == "trt_cp") |> dplyr::distinct(pert)
brd_cpds <- all_cpds |> dplyr::filter(grepl("^BRD-", pert))
cat("Total trt_cp compounds (manuscript 30,456):", nrow(all_cpds), "\n")
cat("Non-BRD genetic-perturbation entries removed (manuscript 45):",
    nrow(all_cpds) - nrow(brd_cpds), "\n")
cat("BRD-only compounds remaining (manuscript 30,411):", nrow(brd_cpds), "\n")

thp1 <- res |> dplyr::filter(cell=="THP1", type=="trt_cp", grepl("^BRD-", pert))
cat("THP1 non-zero, BRD-only compounds (manuscript 597):",
    sum(thp1$NCS!=0 & !is.na(thp1$NCS)), "\n")
cat("THP1 strong reversers NCS<-0.5 (manuscript 288):", sum(thp1$NCS < -0.5, na.rm=TRUE), "\n")
cat("THP1 strong reversers as % (manuscript 48.2%):",
    round(100*sum(thp1$NCS< -0.5,na.rm=TRUE)/sum(thp1$NCS!=0 & !is.na(thp1$NCS)),1), "%\n\n")

## ---- 11. FINAL CANDIDATE LIST (from script 03's saved objects, Dr. Nila's review) ----
cat("---- 11. FINAL CANDIDATE LIST ----\n")
step2 <- readRDS(file.path(PROCESSED_DIR, "lincs2_step2_significant.rds"))
final_candidates <- readRDS(file.path(PROCESSED_DIR, "lincs2_final_candidates.rds"))
cat("Step-2 significant reversers, NCS<-0.5 AND FDR<0.05 (manuscript 138):", nrow(step2), "\n")
cat("Final candidates - all 3 criteria, FDA-approved (manuscript 8):", nrow(final_candidates), "\n")
print(final_candidates |> dplyr::select(pert, pert_iname, thp1_NCS, thp1_FDR, approval_status))
cat("\n")

## ---- 12. MOA/TARGET ANNOTATION (LINCS2-embedded columns; DIAGNOSTIC) ----
cat("---- 12. MOA/TARGET ANNOTATION COVERAGE ----\n")
cat("NOTE: the prior pipeline (supp_table3_thp1_ranked_moa.csv) depended on\n")
cat("an F:-drive file, an undefined object, and a tau/ranking mismatch with\n")
cat("the corrected script - retired. Book1.xlsx was recovered and inspected:\n")
cat("only 5 rows, only 2 with real target/MOA content (rosuvastatin: HMGCR /\n")
cat("HMGCR inhibitor; GDC-0068: AKT1 / Akt inhibitor), both already known\n")
cat("independently. This section checks LINCS2's OWN embedded columns\n")
cat("instead. Coverage below is a FRESH finding, not verified before this\n")
cat("run - it will very likely NOT match old figures (184/288, 8 of 13 PDE),\n")
cat("since both the annotation source and the underlying scores changed.\n\n")

thp1_strong <- res |>
  dplyr::filter(type == "trt_cp", toupper(cell) == "THP1", grepl("^BRD-", pert), NCS < -0.5)
cat("THP1-specific strong reversers, BRD-only (manuscript 288):", nrow(thp1_strong), "\n")
cat(sprintf("  Non-NA t_gn_sym: %d (%.1f%%)\n", sum(!is.na(thp1_strong$t_gn_sym)),
            100*mean(!is.na(thp1_strong$t_gn_sym))))
cat(sprintf("  Non-NA MOAss:    %d (%.1f%%)\n", sum(!is.na(thp1_strong$MOAss)),
            100*mean(!is.na(thp1_strong$MOAss))))
cat(sprintf("  Non-NA PCIDss:   %d (%.1f%%)\n", sum(!is.na(thp1_strong$PCIDss)),
            100*mean(!is.na(thp1_strong$PCIDss))))

if (sum(!is.na(thp1_strong$MOAss)) == 0) {
  cat("\nRESULT: MOAss is entirely empty for THP1 strong reversers in this\n")
  cat("LINCS2 release. The 'X of 288 (Y%) MOA-annotated' manuscript claim\n")
  cat("cannot be reproduced from embedded columns as written. Options:\n")
  cat("  (a) drop the coverage claim from the manuscript entirely,\n")
  cat("  (b) restrict it explicitly to the reviewed final-candidate set\n")
  cat("      (already fully identified by hand via Dr. Nila's review),\n")
  cat("  (c) source a proper external CLUE.io Touchstone bulk export.\n")
  cat("Decide with Dr. Nila before this goes back into the text.\n")
} else {
  cat("\nMOAss IS populated for a subset - distinct values found:\n")
  print(thp1_strong |> dplyr::filter(!is.na(MOAss)) |>
          dplyr::count(MOAss, sort = TRUE) |> head(20))
}
cat("\n")

## ---- 13. PDE4 INHIBITOR CLASS CHECK ----
cat("---- 13. PDE4 INHIBITOR CLASS (roflumilast, oglemilast) ----\n")
cat("Book1.xlsx confirmed NO target/MOA entry for oglemilast (blank).\n")
cat("Its 'PDE4 inhibitor' classification is expert-asserted (correct,\n")
cat("standard pharmacology) but NOT database-annotated. If the manuscript\n")
cat("states PDE4-class compounds were identified 'via CLUE/MOA annotation,'\n")
cat("that wording needs to change - oglemilast was identified by name.\n")
cat("Embedded LINCS2 annotation for both, for comparison:\n")
pde_check <- res |>
  dplyr::filter(pert %in% c("BRD-K03194791", "BRD-K92303087"), toupper(cell) == "THP1") |>
  dplyr::select(pert, cell, NCS, t_gn_sym, MOAss)
print(pde_check)
cat("\n")

## ---- 14. JAK-PATHWAY EVIDENCE: ORIGINAL HIT + BARICITINIB ----
cat("---- 14. JAK-PATHWAY COMPOUNDS ----\n")
cat("Original JAK2 inhibitor (manuscript: HEPG2 -1.33; A375 -0.946; ASC -0.741; not profiled in THP1)\n")
jak <- res |> dplyr::filter(pert=="BRD-K06817181", type=="trt_cp") |>
  dplyr::arrange(NCS) |> dplyr::select(cell, NCS)
print(head(jak, 6))
cat("THP1 profiled for this compound? (manuscript: not profiled):",
    "THP1" %in% jak$cell, "\n\n")

cat("Baricitinib (BRD-K53581288) - NEW finding this session, FDA/EMA-approved\n")
cat("JAK1/2 inhibitor, found among THP1 significant reversers. Unlike the\n")
cat("compound above, this one IS profiled in THP1 and IS an approved drug -\n")
cat("a stronger, more direct convergent JAK-pathway data point:\n")
baricitinib <- res |> dplyr::filter(pert=="BRD-K53581288", type=="trt_cp") |>
  dplyr::arrange(NCS) |> dplyr::select(cell, NCS, WTCS_FDR)
print(head(baricitinib, 10))
cat("THP1 profiled for baricitinib?:", "THP1" %in% baricitinib$cell, "\n\n")

## ---- 15. LEAVE-ONE-OUT STABILITY (manuscript: UP 0.714, DOWN 0.733, mean 0.724) ----
cat("---- 15. LOO STABILITY ----\n")
loo <- readRDS(file.path(PROCESSED_DIR, "phase6_loo_stability.rds"))
cat("Mean Jaccard UP (manuscript 0.714):", round(mean(loo$jaccard_up),3), "\n")
cat("Mean Jaccard DOWN (manuscript 0.733):", round(mean(loo$jaccard_dn),3), "\n")
cat("Combined mean (manuscript 0.724):", round(mean(c(loo$jaccard_up,loo$jaccard_dn)),3), "\n")
persist <- read.csv(file.path(PROCESSED_DIR, "phase6_loo_gene_persistence.csv"))
cat("Genes tracked (manuscript 103):", nrow(persist[persist$in_sig,]), "\n")
cat("Genes stable in 9/9 (manuscript 96):", sum(persist$n_iters==9 & persist$in_sig), "\n\n")

## ---- 16. FIGURE OUTPUT: EXISTENCE, FRESHNESS, KEY NUMBERS ----
cat("---- 16. FIGURE OUTPUT: EXISTENCE, FRESHNESS, KEY NUMBERS ----\n")
cat("NOTE: figure subtitles in 06_figures_all.R are computed LIVE from the\n")
cat("same data objects this audit reads - there is no separate 'baked-in\n")
cat("number' to fall out of sync. What CAN go stale is the image FILE\n")
cat("itself, if underlying data changes but figures aren't regenerated.\n")
cat("This checks existence + freshness, and reprints the key numbers each\n")
cat("panel's subtitle should show, to eyeball against the actual PNG.\n\n")

FIGURES_DIR <- here::here("figures")
CLEAN_DIR   <- file.path(FIGURES_DIR, "clean")

fig_deps <- list(
  "fig1_panelA_umap.png"                  = "seu_monocyte_subclustered.rds",
  "fig1_panelB_isg_violin.png"            = "seu_monocyte_subclustered.rds",
  "fig1_panelC_moca_box.png"              = "seu_monocyte_subclustered.rds",
  "fig1_panelD_volcano.png"               = "de_ifn_subtype_results.rds",
  "fig2_panelA_thp1_ncs.png"              = "lincs2_gess_result.rds",
  "fig2_panelB_thp1_hits.png"             = c("lincs2_compound_summary.rds","lincs2_final_candidates.rds"),
  "fig3_panelA_rosuvastatin_celllines.png"= "lincs2_gess_result.rds",
  "fig3_panelB_roflumilast_celllines.png" = "lincs2_gess_result.rds",
  "fig4_panelA_loo_stability.png"         = "phase6_loo_stability.rds",
  "fig4_panelB_gene_persistence.png"      = "phase6_loo_gene_persistence.csv",
  "suppfig_S1_panelA_pbmc_umap.png"       = "pbmc_annotated_sce.rds",
  "suppfig_S1_panelB_marker_dotplot.png"  = "pbmc_annotated_sce.rds",
  "suppfig_S1_panelC_subcluster_comp.png" = "seu_monocyte_subclustered.rds",
  "suppfig_S2_heatmap.png"                = "seu_monocyte_subclustered.rds",
  "suppfig_S3_bulk_vs_subcluster.png"     = "seu_monocyte_subclustered.rds",
  "suppfig_S4_panelA_population.png"      = "phase2_benchmark_v4.rds",
  "suppfig_S4_panelB_burden.png"          = "phase2_benchmark_v4.rds",
  "suppfig_S4_panelC_overlap.png"         = "phase2_benchmark_v4.rds",
  "suppfig_S5_baricitinib_celllines.png"  = "lincs2_gess_result.rds"
)

cat(sprintf("%-42s %-8s %-10s %s\n", "FIGURE", "EXISTS", "STATUS", "SOURCE DATA"))
cat(strrep("-", 100), "\n")
n_missing <- 0; n_stale <- 0
for (fig_name in names(fig_deps)) {
  fig_path <- file.path(CLEAN_DIR, fig_name)
  if (!file.exists(fig_path)) {
    n_missing <- n_missing + 1
    cat(sprintf("%-42s %-8s %-10s\n", fig_name, "NO", "MISSING"))
    next
  }
  fig_mtime <- file.info(fig_path)$mtime
  deps <- fig_deps[[fig_name]]
  dep_mtimes <- sapply(deps, function(d) {
    p <- file.path(PROCESSED_DIR, d)
    if (file.exists(p)) file.info(p)$mtime else as.POSIXct(NA)
  })
  newest_dep <- suppressWarnings(max(dep_mtimes, na.rm = TRUE))
  is_stale <- fig_mtime < newest_dep
  if (is_stale) n_stale <- n_stale + 1
  cat(sprintf("%-42s %-8s %-10s %s\n", fig_name, "yes",
              if (is_stale) "STALE" else "current", paste(deps, collapse = ", ")))
}
cat(strrep("-", 100), "\n")
cat(sprintf("Missing: %d | Stale: %d | OK: %d\n\n",
            n_missing, n_stale, length(fig_deps) - n_missing - n_stale))
if (n_missing > 0 || n_stale > 0) {
  cat("ACTION: re-run analysis/06_figures_all.R before treating figures as final.\n\n")
} else {
  cat("All checked figures exist and are newer than every data object they\n")
  cat("depend on.\n\n")
}

cat("---- Key numbers each panel's subtitle should show (cross-check by eye) ----\n\n")

cat("Fig 2A (THP1 NCS distribution):\n")
cat(sprintf("  Expect: '%d compounds with non-zero NCS | %d strong reversers (NCS < -0.5)'\n",
            sum(thp1$NCS != 0 & !is.na(thp1$NCS)), sum(thp1$NCS < -0.5, na.rm = TRUE)))

cat("\nSupp Fig S1 (PBMC UMAP):\n")
cat(sprintf("  Expect monocyte count in red = %d (must equal seu_monocyte_subclustered.rds ncol)\n",
            ncol(seu)))

cat("\nSupp Fig S4 Panel A (population-level score):\n")
cat(sprintf("  Expect medians: iNPH %.4f, SLE-high %.4f, HC-high %.4f\n",
            A$med_pop["inph_c5"], A$med_pop["sle_hi"], A$med_pop["hc_hi"]))
cat(sprintf("  Expect: iNPH vs SLE-high p = %.3g | iNPH vs HC-high p = %.3g\n",
            A$p_c5_sle, A$p_c5_hc))

cat("\nSupp Fig S4 Panel B (burden):\n")
cat(sprintf("  Expect medians: HC %.1f%%, SLE %.1f%%, iNPH %.1f%%\n",
            A$burden_med["hc"], A$burden_med["sle"], A$burden_med["inph"]))

cat("\nSupp Fig S4 Panel C (overlap):\n")
cat(sprintf("  Expect: %.1f%% of cluster-5 cells clear the shared IFN-high threshold\n",
            A$overlap["pct_c5_is_ifnhi"]))

cat("\nSupp Fig S5 (baricitinib):\n")
compound_summary_full <- readRDS(file.path(PROCESSED_DIR, "lincs2_compound_summary.rds"))
bari_rank_check <- compound_summary_full |>
  dplyr::filter(!is.na(thp1_NCS)) |> dplyr::arrange(thp1_NCS) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::filter(pert == "BRD-K53581288") |> dplyr::pull(rank)
cat(sprintf("  Expect: ranked %d by THP1 NCS\n", bari_rank_check))
cat("\n")

cat("=================== END AUDIT ===================\n")
