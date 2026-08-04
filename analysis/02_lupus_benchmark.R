# =============================================================================
# Phase 2 - Cross-disease benchmark, v4 (Dr. Nila's donor-balanced refinement)
# =============================================================================
# Changes from v3 (02_lupus_benchmark_v3_nila.R), per Dr. Nila 2026-07 review:
#
#   1. PER-GROUP TOP-QUARTILE INTENSITY ARM REMOVED ENTIRELY.
#      v3 showed this is structurally broken: a top-25%-of-own-distribution
#      threshold means something different in a ~4%-active population (iNPH)
#      than in a ~70%-active one (SLE) - the iNPH "IFN-high" set was 83% cells
#      that are not cluster 5, dragging the median toward / below zero. Not a
#      coding bug; a property of the specification. Confirmed with Dr. Nila -
#      do not report these numbers anywhere in the manuscript.
#
#   2. BURDEN THRESHOLD IS NOW DONOR-BALANCED. Previously the "shared" IFN-high
#      threshold was the 75th percentile of ISG score across all POOLED HC
#      cells, which implicitly weights each HC donor by how many cells were
#      recovered for them. Now: compute each HC donor's OWN 75th-percentile
#      ISG score, then take the MEDIAN of those donor-level values. One donor,
#      one vote. This single scalar threshold is then applied uniformly to
#      every cell in every group.
#
#   3. ONE THRESHOLD NOW DRIVES BOTH THE OVERLAP CHECK AND THE BURDEN ANALYSIS
#      (v3 used two different threshold definitions for these; now consistent).
#
#   4. BURDEN IS TESTED AT DONOR LEVEL (Wilcoxon, iNPH vs SLE and iNPH vs HC),
#      consistent with donor-level inference used throughout this pipeline.
#      This is an addition beyond Dr. Nila's explicit request - flag back if
#      you'd rather report burden as plain percentages without a test.
#
# Unchanged from v3: population-level comparison (unfiltered cluster 5 vs
# SLE-high / HC-high) is SETTLED and not revisited here - three independent
# runs (150-cap, uncapped v2, v3) agree. Selector-excluded scoring; raw counts
# only; chunked AUCell streaming, no per-donor cap; AUCell primary engine,
# mean(UP)-mean(DOWN) sensitivity engine.
#
# CAVEAT CARRIED INTO THE MANUSCRIPT: burden for iNPH crosses the Duy/Perez
# dataset boundary and is more batch-sensitive than the within-Perez SLE-vs-HC
# comparison. Per Dr. Nila: present as exploratory/supportive, not headline.
#
# Output: phase2_benchmark_v4.rds
# =============================================================================

suppressMessages({
  library(here); library(dplyr); library(tibble); library(Matrix)
  library(Seurat); library(AUCell); library(rhdf5)
})

PROCESSED_DIR <- here::here("data", "processed")
H5 <- here::here("data", "raw", "perez2022_lupus",
                 "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad")
stopifnot(file.exists(H5))
set.seed(1)

READ_BLOCK   <- 10000L
AUC_MAXRANK  <- 1500L
ISG_QUANTILE <- 0.75
MIN_CELLS    <- 10L

INPH_CELLTYPE_COL    <- NA_character_
INPH_CLASSICAL_VALUE <- "classical monocyte"

# -----------------------------------------------------------------------------
# Section 0 . Signature + selector panel  (unchanged from v3)
# -----------------------------------------------------------------------------
sig       <- readRDS(file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
# The canonical 15-gene ISG selector panel, used to define interferon-high
# SLE/HC donors and excluded from the scoring signature to prevent
# selector/scoring circularity (v4 methodology). Previously loaded from
# phase2_donorlevel_benchmark.rds, an object from an earlier script version
# that no longer exists anywhere in the repository - this script has
# likely never run successfully from scratch since that file went missing.
# Recovered directly from phase2_benchmark_v4.rds$excluded_genes, the
# trusted, already-verified object every downstream table and figure in
# this project actually reads from, so this is confirmed to be the real
# panel, not a reconstruction.
ISG_PANEL <- c("ISG15", "IFI6", "IFI44", "IFI44L", "MX1", "MX2", "OAS1",
               "OAS2", "OAS3", "IFIT1", "IFIT3", "RSAD2", "STAT1", "USP18",
               "LY6E")
stopifnot(length(ISG_PANEL) == 15)

up_excl  <- setdiff(sig$up, ISG_PANEL)
down_set <- sig$down
stopifnot(length(up_excl) == length(sig$up) - 15L)
cat(sprintf("Selector-excluded UP: %d | DOWN: %d | ISG selector: %d\n",
            length(up_excl), length(down_set), length(ISG_PANEL)))
gene_sets <- list(isg_selector = ISG_PANEL, up_excl = up_excl, down = down_set)

# -----------------------------------------------------------------------------
# Section 1 . /obs metadata + CSR pointers + gene names  (unchanged from v3)
# -----------------------------------------------------------------------------
read_obs_cat <- function(h5, field) {
  g <- rhdf5::h5read(h5, paste0("/obs/", field))
  if (is.list(g) && !is.null(g$categories))
    as.character(g$categories)[as.integer(g$codes) + 1L] else as.character(g)
}
cat("\nReading cell metadata (/obs) ...\n")
cell_type <- read_obs_cat(H5, "cell_type")
disease   <- read_obs_cat(H5, "disease")
donor     <- read_obs_cat(H5, "ind_cov")
n_cells   <- length(cell_type)
is_mono   <- cell_type == "classical monocyte"
cat(sprintf("Total cells: %d | classical monocytes: %d (all used - no cap)\n",
            n_cells, sum(is_mono)))

indptr <- as.numeric(rhdf5::h5read(H5, "/raw/X/indptr"))
stopifnot(length(indptr) == n_cells + 1L)
feat       <- rhdf5::h5read(H5, "/raw/var/feature_name")
gene_names <- make.unique(as.character(feat$categories)[as.integer(feat$codes) + 1L])
n_genes    <- length(gene_names)

# -----------------------------------------------------------------------------
# Section 2 . Scoring helpers  (unchanged from v3)
# -----------------------------------------------------------------------------
auc_scores <- function(counts_mat, sets, maxrank) {
  rk  <- AUCell::AUCell_buildRankings(counts_mat, plotStats = FALSE, verbose = FALSE)
  auc <- AUCell::AUCell_calcAUC(sets, rk, aucMaxRank = maxrank, verbose = FALSE)
  t(AUCell::getAUC(auc))
}
lognorm_cols <- function(cm) {
  ls <- Matrix::colSums(cm); ls[ls == 0] <- 1
  log1p(Matrix::t(Matrix::t(cm) / ls) * 1e4)
}
meandiff_vec <- function(logm, up, dn) {
  ui <- intersect(rownames(logm), up); di <- intersect(rownames(logm), dn)
  Matrix::colMeans(logm[ui, , drop = FALSE]) - Matrix::colMeans(logm[di, , drop = FALSE])
}
isg_meanlog <- function(logm, panel) {
  pidx <- intersect(rownames(logm), panel); Matrix::colMeans(logm[pidx, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# Section 3 . Stream ALL Perez monocytes in blocks  (unchanged from v3)
# -----------------------------------------------------------------------------
mono_global_id <- integer(0)
A_score <- numeric(0); A_isg <- numeric(0)
B_score <- numeric(0); B_isg <- numeric(0)

block_starts <- seq(1L, n_cells, by = READ_BLOCK)
cat(sprintf("\nStreaming %d blocks ...\n", length(block_starts)))
for (bi in seq_along(block_starts)) {
  c1 <- block_starts[bi]; c2 <- min(c1 + READ_BLOCK - 1L, n_cells)
  mono_cols <- which(is_mono[c1:c2])
  if (length(mono_cols) == 0L) next

  s <- indptr[c1] + 1L; e <- indptr[c2 + 1L]
  d   <- as.numeric(rhdf5::h5read(H5, "/raw/X/data",    index = list(s:e)))
  gi  <- as.integer(rhdf5::h5read(H5, "/raw/X/indices", index = list(s:e))) + 1L
  nnz <- diff(indptr[c1:(c2 + 1L)])
  jj  <- rep.int(seq_len(c2 - c1 + 1L), nnz)
  blk <- Matrix::sparseMatrix(i = gi, j = jj, x = d,
                              dims = c(n_genes, c2 - c1 + 1L),
                              dimnames = list(gene_names, NULL), repr = "C")
  blk <- blk[, mono_cols, drop = FALSE]

  a  <- auc_scores(blk, gene_sets, AUC_MAXRANK)
  lg <- lognorm_cols(blk)
  A_score <- c(A_score, a[, "up_excl"] - a[, "down"])
  A_isg   <- c(A_isg,   a[, "isg_selector"])
  B_score <- c(B_score, meandiff_vec(lg, up_excl, down_set))
  B_isg   <- c(B_isg,   isg_meanlog(lg, ISG_PANEL))
  mono_global_id <- c(mono_global_id, (c1:c2)[mono_cols])

  rm(blk, lg, a, d, gi, jj, nnz)
  if (bi %% 20L == 0L) { gc(FALSE); cat(sprintf("  block %d/%d\n", bi, length(block_starts))) }
}
lupus <- tibble(donor = donor[mono_global_id], disease = disease[mono_global_id],
                A_score = A_score, A_isg = A_isg, B_score = B_score, B_isg = B_isg)
is_sle <- lupus$disease == "systemic lupus erythematosus"
is_hc  <- lupus$disease == "normal"
cat(sprintf("Scored: SLE %d | HC %d monocytes\n", sum(is_sle), sum(is_hc)))

# -----------------------------------------------------------------------------
# Section 4 . iNPH monocytes - full compartment  (unchanged from v3)
# -----------------------------------------------------------------------------
seu <- readRDS(file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"))
keep_cells <- rep(TRUE, ncol(seu))
if (!is.na(INPH_CELLTYPE_COL)) {
  stopifnot(INPH_CELLTYPE_COL %in% colnames(seu@meta.data))
  keep_cells <- as.character(seu@meta.data[[INPH_CELLTYPE_COL]]) == INPH_CLASSICAL_VALUE
  cat(sprintf("Restricted to '%s' == '%s': %d of %d cells\n",
              INPH_CELLTYPE_COL, INPH_CLASSICAL_VALUE, sum(keep_cells), ncol(seu)))
}
inph_counts <- GetAssayData(seu, layer = "counts")[, keep_cells, drop = FALSE]
is_c5  <- as.logical(seu$sub_res05[keep_cells] == "5")
a_inph <- auc_scores(inph_counts, gene_sets, AUC_MAXRANK)
lg_inph <- lognorm_cols(inph_counts)
inph <- tibble(donor = as.character(seu$patient[keep_cells]), is_c5 = is_c5,
               A_score = a_inph[, "up_excl"] - a_inph[, "down"],
               A_isg   = a_inph[, "isg_selector"],
               B_score = meandiff_vec(lg_inph, up_excl, down_set),
               B_isg   = isg_meanlog(lg_inph, ISG_PANEL))
cat(sprintf("iNPH compartment: %d cells | cluster 5: %d\n", nrow(inph), sum(is_c5)))

# -----------------------------------------------------------------------------
# Section 5 . Donor-balanced shared threshold -> overlap + burden (donor-level)
# -----------------------------------------------------------------------------
donor_agg <- function(df, col) {
  df |> dplyr::group_by(donor) |>
    dplyr::summarise(mean_sig = mean(.data[[col]]), n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n >= MIN_CELLS)
}
Wp <- function(a, b) if (nrow(a) > 0 && length(b) > 0)
  suppressWarnings(wilcox.test(a$mean_sig, b))$p.value else NA_real_

donor_q75 <- function(df, isg_col) {
  df |> dplyr::group_by(donor) |>
    dplyr::filter(dplyr::n() >= MIN_CELLS) |>
    dplyr::summarise(q75 = quantile(.data[[isg_col]], ISG_QUANTILE, na.rm = TRUE),
                      .groups = "drop")
}
donor_burden <- function(df, isg_col, thr) {
  df |> dplyr::group_by(donor) |>
    dplyr::summarise(burden_pct = 100 * mean(.data[[isg_col]] >= thr, na.rm = TRUE),
                      n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n >= MIN_CELLS)
}

build_engine <- function(sc, isg, label) {
  # --- POPULATION-LEVEL comparison: SETTLED, unchanged, not revisited ---
  thr_sle_isg <- as.numeric(quantile(lupus[[isg]][is_sle], ISG_QUANTILE, na.rm = TRUE))
  thr_hc_isg  <- as.numeric(quantile(lupus[[isg]][is_hc],  ISG_QUANTILE, na.rm = TRUE))
  sle_hi  <- donor_agg(lupus[is_sle & lupus[[isg]] >= thr_sle_isg, ], sc)
  hc_hi   <- donor_agg(lupus[is_hc  & lupus[[isg]] >= thr_hc_isg,  ], sc)
  inph_c5 <- donor_agg(inph[inph$is_c5, ], sc)

  # --- Dr. Nila's donor-balanced shared threshold ---
  hc_thr_tbl <- donor_q75(lupus[is_hc, ], isg)
  thr_shared <- median(hc_thr_tbl$q75, na.rm = TRUE)

  # --- OVERLAP using the shared threshold (now consistent with burden) ---
  inph_hi_mask <- inph[[isg]] >= thr_shared
  ov <- c(n_ifnhi = sum(inph_hi_mask), n_c5 = sum(inph$is_c5),
          n_both  = sum(inph_hi_mask & inph$is_c5),
          pct_ifnhi_in_c5 = 100 * sum(inph_hi_mask & inph$is_c5) / max(1, sum(inph_hi_mask)),
          pct_c5_is_ifnhi = 100 * sum(inph_hi_mask & inph$is_c5) / max(1, sum(inph$is_c5)))

  # --- BURDEN, donor-level, shared threshold ---
  burden_hc   <- donor_burden(lupus[is_hc, ],  isg, thr_shared)
  burden_sle  <- donor_burden(lupus[is_sle, ], isg, thr_shared)
  burden_inph <- donor_burden(inph,            isg, thr_shared)
  p_bd_sle <- if (nrow(burden_inph) > 0 && nrow(burden_sle) > 0)
    suppressWarnings(wilcox.test(burden_inph$burden_pct, burden_sle$burden_pct))$p.value else NA_real_
  p_bd_hc  <- if (nrow(burden_inph) > 0 && nrow(burden_hc) > 0)
    suppressWarnings(wilcox.test(burden_inph$burden_pct, burden_hc$burden_pct))$p.value else NA_real_

  list(label = label, thr_shared = thr_shared,
       med_pop = c(inph_c5 = median(inph_c5$mean_sig), sle_hi = median(sle_hi$mean_sig),
                   hc_hi = median(hc_hi$mean_sig)),
       n_pop = c(inph_c5 = nrow(inph_c5), sle_hi = nrow(sle_hi), hc_hi = nrow(hc_hi)),
       p_c5_sle = Wp(inph_c5, sle_hi$mean_sig), p_c5_hc = Wp(inph_c5, hc_hi$mean_sig),
       # Donor-level tables, kept so the population comparison can be PLOTTED
       # (previously computed but discarded - only summary stats were saved,
       # which is why S4 initially had no population-level visual panel).
       pop_tables = list(inph_c5 = inph_c5, sle_hi = sle_hi, hc_hi = hc_hi),
       overlap = ov,
       burden_tables = list(hc = burden_hc, sle = burden_sle, inph = burden_inph),
       burden_med = c(hc = median(burden_hc$burden_pct), sle = median(burden_sle$burden_pct),
                      inph = median(burden_inph$burden_pct)),
       burden_n = c(hc = nrow(burden_hc), sle = nrow(burden_sle), inph = nrow(burden_inph)),
       p_burden_sle = p_bd_sle, p_burden_hc = p_bd_hc)
}

A <- build_engine("A_score", "A_isg", "AUCell (raw counts) [PRIMARY]")
B <- build_engine("B_score", "B_isg", "mean(UP)-mean(DOWN) log-norm [SENSITIVITY]")

# -----------------------------------------------------------------------------
# Section 6 . Report
# -----------------------------------------------------------------------------
report <- function(E) {
  cat(sprintf("\n--- %s ---\n", E$label))
  cat(sprintf("  Shared donor-balanced IFN-high threshold (median of HC donors' own 75th pct): %.4f\n",
              E$thr_shared))
  cat("  POPULATION-LEVEL (settled - unchanged from prior runs, not revisited)\n")
  cat(sprintf("    medians  iNPH-c5 / SLE-hi / HC-hi : %.4f / %.4f / %.4f\n",
              E$med_pop["inph_c5"], E$med_pop["sle_hi"], E$med_pop["hc_hi"]))
  cat(sprintf("    n donors iNPH-c5 / SLE-hi / HC-hi : %d / %d / %d\n",
              E$n_pop["inph_c5"], E$n_pop["sle_hi"], E$n_pop["hc_hi"]))
  cat(sprintf("    iNPH-c5 vs SLE-hi p=%.3g | vs HC-hi p=%.3g\n", E$p_c5_sle, E$p_c5_hc))
  cat("  OVERLAP (shared threshold vs cluster 5 - validation of subcluster identity)\n")
  cat(sprintf("    IFN-high cells %d | cluster-5 %d | both %d\n",
              E$overlap["n_ifnhi"], E$overlap["n_c5"], E$overlap["n_both"]))
  cat(sprintf("    %.1f%% of IFN-high cells are in cluster 5 | %.1f%% of cluster 5 is IFN-high\n",
              E$overlap["pct_ifnhi_in_c5"], E$overlap["pct_c5_is_ifnhi"]))
  cat("  BURDEN (donor-level %% of monocytes above shared threshold; exploratory/supportive)\n")
  cat(sprintf("    median burden  HC / SLE / iNPH : %.1f%% / %.1f%% / %.1f%%\n",
              E$burden_med["hc"], E$burden_med["sle"], E$burden_med["inph"]))
  cat(sprintf("    n donors       HC / SLE / iNPH : %d / %d / %d\n",
              E$burden_n["hc"], E$burden_n["sle"], E$burden_n["inph"]))
  cat(sprintf("    iNPH vs SLE burden p=%.3g | iNPH vs HC burden p=%.3g\n",
              E$p_burden_sle, E$p_burden_hc))
  cat("    NOTE: iNPH burden crosses the Duy/Perez dataset boundary - batch-\n")
  cat("          sensitive. Report as exploratory/supportive, not headline.\n")
}
cat("\n================ v4 BENCHMARK (donor-balanced threshold) ================\n")
report(A); report(B)

saveRDS(list(engine_A_primary = A, engine_B_sensitivity = B,
             params = list(read_block = READ_BLOCK, auc_maxrank = AUC_MAXRANK,
                           isg_quantile = ISG_QUANTILE, min_cells = MIN_CELLS,
                           capped = FALSE,
                           threshold_basis = "donor-balanced: median of HC donors' own 75th-pct ISG score",
                           per_group_intensity_arm = "REMOVED - structurally broken, see header notes"),
             excluded_genes = ISG_PANEL),
        file.path(PROCESSED_DIR, "phase2_benchmark_v4.rds"))
cat("\nSaved: phase2_benchmark_v4.rds\n")
