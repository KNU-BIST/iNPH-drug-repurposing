# =============================================================================
# Phase 2 — Cross-disease contextualization of the iNPH IFN state
# Script: 02_cross_disease_contextualization.R
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
# Compare the iNPH IFN-activated monocyte state with interferon-high
# classical monocytes from systemic lupus erythematosus and healthy-control
# donors.
#
# Primary analysis
# ----------------
# Donor-level AUCell scoring using the iNPH signature after exclusion of the
# 15 genes used to define interferon-high comparator cells.
#
# Sensitivity analysis
# --------------------
# Mean log-normalized expression of the remaining UP genes minus mean
# log-normalized expression of the DOWN genes.
#
# Additional outputs
# ------------------
# - donor-level population comparison;
# - shared-threshold IFN-high monocyte burden;
# - overlap between cluster 5 and marker-defined IFN-high status.
#
# Output
# ------
# data/processed/phase2_contextualization.rds
# =============================================================================

# -----------------------------------------------------------------------------
# Section 0 . Environment and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "tibble",
  "Matrix",
  "Seurat",
  "AUCell",
  "rhdf5"
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
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(Seurat)
  library(AUCell)
  library(rhdf5)
})

PROCESSED_DIR <- here::here("data", "processed")
dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)

SIGNATURE_PATH <- file.path(
  PROCESSED_DIR,
  "inph_ifn_signature.rds"
)

INPH_SEURAT_PATH <- file.path(
  PROCESSED_DIR,
  "seu_monocyte_subclustered.rds"
)

H5 <- here::here(
  "data",
  "raw",
  "perez2022_lupus",
  "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad"
)

required_inputs <- c(
  SIGNATURE_PATH,
  INPH_SEURAT_PATH,
  H5
)

input_info <- file.info(required_inputs)

missing_or_empty <- required_inputs[
  !file.exists(required_inputs) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(missing_or_empty) > 0L) {
  stop(
    "Required input file(s) missing or empty: ",
    paste(missing_or_empty, collapse = ", ")
  )
}

if (!rhdf5::H5Fis_hdf5(H5)) {
  stop(
    "Perez comparator file is not a valid HDF5/H5AD file: ",
    H5
  )
}

cat(
  "iNPH signature MD5: ",
  unname(tools::md5sum(SIGNATURE_PATH)),
  "\n",
  sep = ""
)

cat(
  "Subclustered monocyte object MD5: ",
  unname(tools::md5sum(INPH_SEURAT_PATH)),
  "\n",
  sep = ""
)

set.seed(1)

READ_BLOCK   <- 10000L
AUC_MAXRANK  <- 1500L
ISG_QUANTILE <- 0.75
MIN_CELLS    <- 10L


# -----------------------------------------------------------------------------
# Section 1 . Signature + selector panel
# -----------------------------------------------------------------------------
sig <- readRDS(SIGNATURE_PATH)

if (
  !is.list(sig) ||
  !all(c("up", "down") %in% names(sig))
) {
  stop(
    "inph_ifn_signature.rds must contain a list with 'up' and 'down' elements."
  )
}

if (
  length(sig$up) != 150L ||
  length(sig$down) != 150L ||
  anyDuplicated(sig$up) ||
  anyDuplicated(sig$down) ||
  length(intersect(sig$up, sig$down)) != 0L
) {
  stop(
    "The iNPH signature does not contain the expected disjoint ",
    "150-UP/150-DOWN gene sets."
  )
}

# Canonical 15-gene ISG selector used in the finalized donor-level analysis.
# These genes define interferon-high SLE and healthy-control monocytes and
# are excluded from the iNPH scoring signature to prevent direct reuse of
# selector genes during scoring.
ISG_PANEL <- c("ISG15", "IFI6", "IFI44", "IFI44L", "MX1", "MX2", "OAS1",
               "OAS2", "OAS3", "IFIT1", "IFIT3", "RSAD2", "STAT1", "USP18",
               "LY6E")

stopifnot(length(ISG_PANEL) == 15L)
stopifnot(all(ISG_PANEL %in% sig$up))

up_excl  <- setdiff(sig$up, ISG_PANEL)
down_set <- sig$down

stopifnot(length(up_excl) == 135L)
stopifnot(length(down_set) == 150L)
stopifnot(length(intersect(up_excl, ISG_PANEL)) == 0L)

cat(sprintf("Selector-excluded UP: %d | DOWN: %d | ISG selector: %d\n",
            length(up_excl), length(down_set), length(ISG_PANEL)))
gene_sets <- list(isg_selector = ISG_PANEL, up_excl = up_excl, down = down_set)

# -----------------------------------------------------------------------------
# Section 2 . /obs metadata + CSR pointers + gene names
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

n_cells <- length(cell_type)

if (
  length(disease) != n_cells ||
  length(donor) != n_cells
) {
  stop(
    "Perez /obs metadata fields have inconsistent lengths."
  )
}

is_mono <- cell_type == "classical monocyte"
n_mono_expected <- 307429L
n_mono_observed <- sum(is_mono)

cat(
  sprintf(
    "Total cells: %d | classical monocytes: %d (all used - no cap)\n",
    n_cells,
    n_mono_observed
  )
)

if (n_mono_observed != n_mono_expected) {
  stop(
    "Perez classical-monocyte count mismatch: expected ",
    n_mono_expected,
    ", observed ",
    n_mono_observed,
    ". Confirm the H5AD version and cell-type annotation before proceeding."
  )
}

if (anyNA(donor[is_mono])) {
  stop(
    "Missing donor identifiers detected among Perez classical monocytes."
  )
}

indptr <- as.numeric(rhdf5::h5read(H5, "/raw/X/indptr"))
stopifnot(length(indptr) == n_cells + 1L)
feat       <- rhdf5::h5read(H5, "/raw/var/feature_name")
gene_names <- make.unique(as.character(feat$categories)[as.integer(feat$codes) + 1L])
n_genes    <- length(gene_names)

# -----------------------------------------------------------------------------
# Section 3 . Scoring helpers
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
  up_present <- intersect(rownames(logm), up)
  dn_present <- intersect(rownames(logm), dn)

  if (length(up_present) == 0L || length(dn_present) == 0L) {
    stop("UP or DOWN scoring genes are absent from the expression matrix.")
  }

  Matrix::colMeans(logm[up_present, , drop = FALSE]) -
    Matrix::colMeans(logm[dn_present, , drop = FALSE])
}

isg_meanlog <- function(logm, panel) {
  panel_present <- intersect(rownames(logm), panel)

  if (length(panel_present) < 3L) {
    stop("Insufficient ISG selector genes found in the expression matrix.")
  }

  Matrix::colMeans(logm[panel_present, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# Section 4 . Stream ALL Perez monocytes in blocks
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
lupus <- tibble(
  donor = donor[mono_global_id],
  disease = disease[mono_global_id],
  A_score = A_score,
  A_isg = A_isg,
  B_score = B_score,
  B_isg = B_isg
)

if (nrow(lupus) != n_mono_expected) {
  stop(
    "Streaming/scoring mismatch: expected ",
    n_mono_expected,
    " Perez classical monocytes, scored ",
    nrow(lupus),
    "."
  )
}

is_sle <- lupus$disease == "systemic lupus erythematosus"
is_hc  <- lupus$disease == "normal"

if (!all(is_sle | is_hc)) {
  unexpected_disease <- sort(
    unique(
      lupus$disease[!(is_sle | is_hc)]
    )
  )

  stop(
    "Unexpected disease label(s) among Perez classical monocytes: ",
    paste(unexpected_disease, collapse = ", ")
  )
}

cat(
  sprintf(
    "Scored: SLE %d | HC %d monocytes\n",
    sum(is_sle),
    sum(is_hc)
  )
)

# -----------------------------------------------------------------------------
# Section 5 . iNPH monocytes - full compartment
# -----------------------------------------------------------------------------
seu <- readRDS(INPH_SEURAT_PATH)

if (!inherits(seu, "Seurat")) {
  stop(
    "seu_monocyte_subclustered.rds must contain a Seurat object; observed class: ",
    paste(class(seu), collapse = ", ")
  )
}

required_inph_meta <- c(
  "patient",
  "sub_res05"
)

missing_inph_meta <- setdiff(
  required_inph_meta,
  colnames(seu@meta.data)
)

if (length(missing_inph_meta) > 0L) {
  stop(
    "Subclustered iNPH object is missing required metadata field(s): ",
    paste(missing_inph_meta, collapse = ", ")
  )
}

if (ncol(seu) != 22197L) {
  stop(
    "Subclustered iNPH object should contain 22,197 canonical monocytes; observed ",
    ncol(seu),
    "."
  )
}

if (length(unique(as.character(seu$patient))) != 9L) {
  stop(
    "Subclustered iNPH object should contain nine patients."
  )
}

is_c5 <- as.character(seu$sub_res05) == "5"

if (sum(is_c5) != 940L) {
  stop(
    "IFN-activated cluster 5 count mismatch: expected 940 cells, observed ",
    sum(is_c5),
    "."
  )
}

inph_counts <- Seurat::GetAssayData(
  seu,
  layer = "counts"
)

a_inph <- auc_scores(
  inph_counts,
  gene_sets,
  AUC_MAXRANK
)

lg_inph <- lognorm_cols(
  inph_counts
)

inph <- tibble(
  donor = as.character(seu$patient),
  is_c5 = is_c5,
  A_score = a_inph[, "up_excl"] - a_inph[, "down"],
  A_isg = a_inph[, "isg_selector"],
  B_score = meandiff_vec(lg_inph, up_excl, down_set),
  B_isg = isg_meanlog(lg_inph, ISG_PANEL)
)

cat(
  sprintf(
    "iNPH compartment: %d cells | cluster 5: %d\n",
    nrow(inph),
    sum(is_c5)
  )
)


# -----------------------------------------------------------------------------
# Section 6 . Donor-balanced shared threshold -> overlap + burden (donor-level)
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
  # Population-level contextual comparison:
  # SLE-high and HC-high subsets are defined using the top quartile of the
  # corresponding comparator-group ISG-score distribution. The complete,
  # unfiltered iNPH cluster 5 population is retained.
  thr_sle_isg <- as.numeric(quantile(lupus[[isg]][is_sle], ISG_QUANTILE, na.rm = TRUE))
  thr_hc_isg  <- as.numeric(quantile(lupus[[isg]][is_hc],  ISG_QUANTILE, na.rm = TRUE))
  sle_hi  <- donor_agg(lupus[is_sle & lupus[[isg]] >= thr_sle_isg, ], sc)
  hc_hi   <- donor_agg(lupus[is_hc  & lupus[[isg]] >= thr_hc_isg,  ], sc)
  inph_c5 <- donor_agg(inph[inph$is_c5, ], sc)

  # Shared burden threshold:
  # Each healthy-control donor contributes its own 75th-percentile ISG score.
  # The median of those donor-specific thresholds is applied uniformly to all
  # cells in all three groups.
  hc_thr_tbl <- donor_q75(lupus[is_hc, ], isg)
  thr_shared <- median(hc_thr_tbl$q75, na.rm = TRUE)

  # --- OVERLAP using the shared donor-balanced threshold ---
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

  list(
    label = label,  thr_shared = thr_shared,
    med_pop = c(inph_c5 = median(inph_c5$mean_sig, na.rm = TRUE),
      sle_hi  = median(sle_hi$mean_sig, na.rm = TRUE),
      hc_hi   = median(hc_hi$mean_sig, na.rm = TRUE)),
    n_pop = c(inph_c5 = nrow(inph_c5), sle_hi  = nrow(sle_hi), hc_hi   = nrow(hc_hi)),
    p_c5_sle = Wp(inph_c5, sle_hi$mean_sig),
    p_c5_hc  = Wp(inph_c5, hc_hi$mean_sig),
    pop_tables = list(inph_c5 = inph_c5,  sle_hi  = sle_hi,  hc_hi   = hc_hi),
    overlap = ov,
    burden_tables = list(hc   = burden_hc, sle  = burden_sle, inph = burden_inph),
    burden_med = c(
      hc   = median(burden_hc$burden_pct, na.rm = TRUE),
      sle  = median(burden_sle$burden_pct, na.rm = TRUE),
      inph = median(burden_inph$burden_pct, na.rm = TRUE)),
    burden_n = c(
      hc   = nrow(burden_hc),
      sle  = nrow(burden_sle),
      inph = nrow(burden_inph)),
    p_burden_sle = p_bd_sle,
    p_burden_hc  = p_bd_hc)
}

A <- build_engine("A_score", "A_isg", "AUCell (raw counts) [PRIMARY]")
B <- build_engine("B_score", "B_isg", "mean(UP)-mean(DOWN) log-norm [SENSITIVITY]")

# -----------------------------------------------------------------------------
# Section 7 . Report
# -----------------------------------------------------------------------------
report <- function(E) {
  cat(sprintf("\n--- %s ---\n", E$label))
  cat(sprintf("  Shared donor-balanced IFN-high threshold (median of HC donors' own 75th pct): %.4f\n",
              E$thr_shared))
  cat("  POPULATION-LEVEL DONOR COMPARISON\n")
  cat(sprintf("    medians  iNPH-c5 / SLE-hi / HC-hi : %.4f / %.4f / %.4f\n",
              E$med_pop["inph_c5"], E$med_pop["sle_hi"], E$med_pop["hc_hi"]))
  cat(sprintf("    n donors iNPH-c5 / SLE-hi / HC-hi : %d / %d / %d\n",
              E$n_pop["inph_c5"], E$n_pop["sle_hi"], E$n_pop["hc_hi"]))
  cat(sprintf("    iNPH-c5 vs SLE-hi p=%.3g | vs HC-hi p=%.3g\n", E$p_c5_sle, E$p_c5_hc))
  cat("  CONCORDANCE (shared threshold vs cluster 5 marker-defined IFN-high status)\n")
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
cat("\n================ CROSS-DISEASE CONTEXTUALIZATION ================\n")
report(A); report(B)

contextualization_result <- list(
  engine_A_primary = A,
  engine_B_sensitivity = B,
  params = list(
    read_block = READ_BLOCK,
    auc_maxrank = AUC_MAXRANK,
    isg_quantile = ISG_QUANTILE,
    min_cells = MIN_CELLS,
    capped = FALSE,
    threshold_basis =
      "donor-balanced: median of HC donors' own 75th-pct ISG score",
    per_group_intensity_arm = "not_estimated",
    per_group_intensity_note = paste(
      "Not estimated because group-specific quantile thresholds select",
      "biologically different populations when IFN-high prevalence differs",
      "substantially across groups."
    )
  ),
  excluded_genes = ISG_PANEL,
  input_provenance = list(
    inph_signature_md5 =
      unname(tools::md5sum(SIGNATURE_PATH)),
    subclustered_monocyte_md5 =
      unname(tools::md5sum(INPH_SEURAT_PATH)),
    perez_h5ad_basename =
      basename(H5),
    perez_classical_monocytes =
      n_mono_expected
  )
)

output_path <- file.path(
  PROCESSED_DIR,
  "phase2_contextualization.rds"
)

saveRDS(
  contextualization_result,
  output_path
)

if (
  !file.exists(output_path) ||
  is.na(file.info(output_path)$size) ||
  file.info(output_path)$size <= 0
) {
  stop(
    "Failed to create the Phase 2 output file: ",
    output_path
  )
}

cat(
  "\nSaved and verified: ",
  output_path,
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 8 . Closing summary
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PHASE 2 COMPLETE\n")
cat("============================================================\n")
cat("  Perez classical monocytes :", nrow(lupus), "\n")
cat("  iNPH canonical monocytes  :", nrow(inph), "\n")
cat("  iNPH cluster 5 cells      :", sum(inph$is_c5), "\n")
cat("  Selector-excluded UP genes:", length(up_excl), "\n")
cat("  DOWN genes                :", length(down_set), "\n")
cat("  Primary scoring           : AUCell\n")
cat("  Sensitivity scoring       : mean(UP)-mean(DOWN)\n")
cat("============================================================\n")
