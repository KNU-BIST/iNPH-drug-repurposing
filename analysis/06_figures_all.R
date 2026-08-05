# =============================================================================
# iNPH Drug Repurposing — ALL FIGURES (main + supplementary), single script
# =============================================================================
#
# Replaces: All_figures_v4.R  +  05_supplementary_figures_V4.R
#
# FIGURE STRUCTURE (agreed with Dr. Nila, 2026-07-28)
# ---------------------------------------------------
# MAIN
#   Fig 1  IFN subcluster identification + DE     A UMAP  B ISG violin
#                                                 C MoCA box  D volcano
#   Fig 2  LINCS2 drug screen                     A NCS histogram
#                                                 B top-30 THP1 reversal hits
#   Fig 3  Lead compound cell-line profiles       A rosuvastatin  B roflumilast
#   Fig 4  Robustness (promoted from supp)        A LOO stability
#                                                 B per-gene persistence
# SUPPLEMENTARY
#   S1  PBMC annotation           A UMAP  B lineage dotplot  C subcluster comp
#   S2  Pseudobulk heatmap (18 samples)
#   S3  Bulk vs subcluster scoring
#   S4  Lupus cross-disease benchmark (v4)  A burden  B cluster-5 overlap
#   S5  Baricitinib cell-line profile
#   S6  Proposed mechanism schematic
#
# WHAT CHANGED vs the two source scripts
# --------------------------------------
# 1. NO targets DEPENDENCY. Both source scripts called targets::tar_read().
#    The targets pipeline implements the abandoned study design and is not
#    published, so figures depending on it could not be reproduced from the
#    repository. Fixed as follows:
#      - main script's `inph_mono` load was dead code (never used) -> removed
#      - S1 now builds its own PBMC object from the public GEO data via
#        R/data_loading.R + R/clustering.R (SECTION 2), cached to disk
#      - S2 pseudobulk now uses counts from the published Seurat object
#
# 2. LUPUS BENCHMARK MOVED main Fig 2 -> supplementary S4, and rebuilt on
#    phase2_benchmark_v4.rds. The old panels read
#    phase2_donorlevel_benchmark.rds / phase2_step1_celllevel_comparison.rds,
#    both computed on scaled data with a circular ISG selector, and carried
#    the subtitle "Comparable to interferon-high SLE monocytes (n.s.)" - the
#    claim that was subsequently disproven. Do not resurrect those objects.
#
# 3. STALE PATHS REMOVED. The main script read lincs_meta from
#    Sys.getenv("INPH_EXTERNAL_DATA") (now deliberately unset, so it resolved
#    to a broken path) and loaded lincs2_compound_named.rds, written by the
#    retired CLUE.io annotation pipeline. Both replaced with project-relative
#    reads of objects written by 03_lincs_connectivity.R.
#
# 4. MOA / DRUG-CLASS COLOURING REMOVED from the THP1 hits panel. Per
#    Dr. Nila 2026-07-28: the aggregate MOA coverage and class-split analyses
#    are dropped, because LINCS2's embedded MOA columns are empty (0 of 288
#    strong reversers annotated) and the original 184/288 figure could not be
#    traced to any source. Compounds are now coloured by the manuscript tier
#    (primary lead / FDA-approved / not approved), which is reproducible from
#    the verified approval table in 03_lincs_connectivity.R.
#
# 5. JAK2 INHIBITOR REPLACED BY BARICITINIB (S5). BRD-K06817181 was not
#    profiled in THP1 and is not an approved drug. Baricitinib is FDA/EMA
#    approved, profiled in THP1 (NCS -1.23, FDR 8.17e-04), and ranks 39th of
#    all compounds with a THP1 score. Framed as a post hoc, mechanistically
#    aligned secondary candidate - NOT as JAK-class convergence, and the
#    outside-the-top-30 selection is stated in the subtitle.
#
# 6. JAK2 BOX REMOVED from the mechanism schematic (S6), which now shows only
#    the two primary leads.
#
# HOW TO RUN
# ----------
#   Section 0  config              (once)
#   Section 1  data loading        (once per session)
#   Section 2  PBMC build for S1   (once ever - caches to disk; slow, see note)
#   Sections 3-4  figures
#
# OUTPUT: every figure saved twice -
#   figures/with_titles/  full version (titles, subtitles, captions)
#   figures/clean/        publication version (titles stripped, tags kept)
# =============================================================================

library(here); library(dplyr); library(ggplot2); library(tidyr)
library(stringr); library(readr); library(ggrepel); library(patchwork)

# =============================================================================
# SECTION 0 · CONFIGURATION
# =============================================================================

PROCESSED_DIR <- here::here("data", "processed")
RAW_DIR       <- here::here("data", "raw")
FIGURES_DIR   <- here::here("figures")
dir.create(FIGURES_DIR, showWarnings = FALSE)

PANEL_DPI    <- 600
COMBINED_DPI <- 600
BASE_SIZE    <- 11

PAL <- list(
  disease    = "#d4624a", diseaseSoft = "#fbe9e7",
  drug       = "#5e9ea0", neutral     = "#94a3b8",
  accent     = "#B8860B", accentSoft  = "#faf5e6",
  mocaHigh   = "#7a8c6a", mocaLow     = "#b04a36",
  good       = "#2d6a4f", goodSoft    = "#d8e8dd",
  navy       = "#0D1B2A", purple      = "#6d28d9",
  purpleLight= "#ede9fe",
  bg         = "white",   text        = "#1a2332",
  textSoft   = "#475569", textMute    = "#94a3b8",
  rule       = "#e2e8f0", intFill     = "#f1f5f9",
  intBorder  = "#cbd5e1",
  inhibit    = "#991b1b", inhibitFill = "#fff1f2",
  resting    = "#9DB4C8", ifn         = "#D85A30"
)

THEME_BASE <- theme_minimal(base_size = BASE_SIZE) +
  theme(
    plot.title       = element_text(face = "bold", size = BASE_SIZE + 1,
                                    color = PAL$text),
    plot.subtitle    = element_text(size = BASE_SIZE - 1.5,
                                    color = PAL$textSoft),
    plot.caption     = element_text(size = BASE_SIZE - 2.5,
                                    color = PAL$textMute),
    panel.grid.major = element_line(color = PAL$rule, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(color = "black", face = "bold",
                                    size  = BASE_SIZE - 1.5),
    axis.title       = element_text(color = "black", face = "bold",
                                    size  = BASE_SIZE - 1),
    legend.text      = element_text(size = BASE_SIZE - 1.5),
    legend.title     = element_text(size = BASE_SIZE - 1, face = "bold"),
    plot.background  = element_rect(fill = PAL$bg, color = NA),
    panel.background = element_rect(fill = PAL$bg, color = NA)
  )

# Subcluster palette, shared by Fig 1A/1B and S1C
CLUSTER_COLORS <- c(
  "0" = PAL$textMute, "1" = "#7a8c6a", "2" = PAL$drug,
  "3" = "#9b59b6", "4" = PAL$accent, "5" = PAL$disease, "6" = "#8B5A3C"
)

# Manuscript tier palette (replaces the retired MOA/drug-class palette)
TIER_COLORS <- c(
  "Primary lead"     = PAL$accent,
  "FDA-approved"     = PAL$drug,
  "Not FDA-approved" = PAL$textMute
)

# ── Dual output: titled + clean ──────────────────────────────────────────────
DIR_TITLED <- file.path(FIGURES_DIR, "with_titles")
DIR_CLEAN  <- file.path(FIGURES_DIR, "clean")
dir.create(DIR_TITLED, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CLEAN,  showWarnings = FALSE, recursive = TRUE)

strip_titles <- function(p) {
  blank <- ggplot2::theme(plot.title    = ggplot2::element_blank(),
                          plot.subtitle = ggplot2::element_blank(),
                          plot.caption  = ggplot2::element_blank())
  if (inherits(p, "patchwork")) {
    p & blank
  } else {
    p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) + blank
  }
}

save_panel <- function(p, name, w, h, dpi = PANEL_DPI) {
  # Defensive: clear any stray open graphics device before writing, since a
  # device left open by an earlier failed S2 pheatmap png()/dev.off() pair
  # in this session could interfere with subsequent writes. Harmless no-op
  # if nothing is open.
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()

  paths <- c(
    file.path(DIR_TITLED, paste0(name, ".pdf")),
    file.path(DIR_TITLED, paste0(name, ".png")),
    file.path(DIR_CLEAN,  paste0(name, ".pdf")),
    file.path(DIR_CLEAN,  paste0(name, ".png"))
  )
  pc <- strip_titles(p)
  plots <- list(p, p, pc, pc)

  before <- sapply(paths, function(f)
    if (file.exists(f)) file.info(f)$mtime else as.POSIXct(NA))

  for (i in seq_along(paths)) {
    if (grepl("\\.png$", paths[i]))
      ggsave(paths[i], plots[[i]], width = w, height = h, dpi = dpi)
    else
      ggsave(paths[i], plots[[i]], width = w, height = h)
  }

  # Verify: file exists AND its timestamp moved forward (or is brand new).
  # This is what actually tells us a fresh file landed, not just that
  # ggsave() returned without erroring.
  after <- sapply(paths, function(f)
    if (file.exists(f)) file.info(f)$mtime else as.POSIXct(NA))
  ok <- file.exists(paths) & (is.na(before) | after > before)

  for (i in seq_along(paths))
    cat(if (ok[i]) "  OK:   " else "  FAIL: ", paths[i], "\n")

  if (!all(ok))
    warning("save_panel('", name, "'): not all 4 files verified fresh - see FAIL lines above.")
  else
    message("  Saved + verified (titled + clean): ", name)
}

save_combined <- function(p, name, w, h, dpi = COMBINED_DPI) {
  save_panel(p, name, w, h, dpi)
}

# =============================================================================
# SECTION 1 · DATA LOADING  (all project-relative; no targets, no env vars)
# =============================================================================
message("Loading data...")

inph_signature <- readRDS(file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
de_res         <- readRDS(file.path(PROCESSED_DIR, "de_ifn_subtype_results.rds"))
seu            <- readRDS(file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"))

# Phase 2 - v4 benchmark ONLY (see header note 2)
bench_v4 <- readRDS(file.path(PROCESSED_DIR, "phase2_benchmark_v4.rds"))
BENCH_A  <- bench_v4$engine_A_primary       # AUCell  (primary engine)
BENCH_B  <- bench_v4$engine_B_sensitivity   # mean-diff (sensitivity)

# Phase 3 - objects written by 03_lincs_connectivity.R
result_table   <- readRDS(file.path(PROCESSED_DIR, "lincs2_gess_result.rds")) |>
  signatureSearch::result()
compound_named <- readRDS(file.path(PROCESSED_DIR, "lincs2_compound_summary.rds"))
final_cands    <- readRDS(file.path(PROCESSED_DIR, "lincs2_final_candidates.rds"))

# compound_named was saved BEFORE the name-backfill fix applied to
# step2_significant in 03_lincs_connectivity.R, so 6 of the top-30 compounds
# still carry the literal placeholder "NAME NOT IN LOCAL DB" - identical
# strings, which break any plot needing unique per-row labels (Fig 2B).
# Backfill real names from step2_significant where available (covers
# methyl-fasudil, ZM-39923, CHEMBL-374350, C-646, identified by hand this
# session), then make any STILL-placeholder name unique via its BRD ID, since
# the two genuinely unidentified compounds both show as "UNIDENTIFIED" even
# after backfill and would otherwise collide with each other.
PLACEHOLDER_NAMES <- c("NAME NOT IN LOCAL DB", "UNIDENTIFIED")
step2_names <- readRDS(file.path(PROCESSED_DIR, "lincs2_step2_significant.rds")) |>
  dplyr::select(pert, pert_iname_fixed = pert_iname)

compound_named <- compound_named |>
  dplyr::left_join(step2_names, by = "pert") |>
  dplyr::mutate(
    pert_iname = dplyr::if_else(
      pert_iname %in% PLACEHOLDER_NAMES & !is.na(pert_iname_fixed) &
        !(pert_iname_fixed %in% PLACEHOLDER_NAMES),
      pert_iname_fixed, pert_iname),
    pert_iname = dplyr::if_else(
      pert_iname %in% PLACEHOLDER_NAMES | is.na(pert_iname),
      paste0("Unidentified (", pert, ")"), pert_iname)) |>
  dplyr::select(-pert_iname_fixed)

# Phase 6 - robustness
loo_df     <- readRDS(file.path(PROCESSED_DIR, "phase6_loo_stability.rds"))
persist_df <- readr::read_csv(
  file.path(PROCESSED_DIR, "phase6_loo_gene_persistence.csv"),
  show_col_types = FALSE)

cell_meta <- data.frame(
  cell       = colnames(seu),
  patient    = seu$patient,
  moca_group = as.character(seu$moca_group),
  subcluster = as.character(seu$sub_res05),
  stringsAsFactors = FALSE
)

# Compound BRD IDs used across panels
BRD_ROSUVASTATIN <- "BRD-K82941592"
BRD_ROFLUMILAST  <- "BRD-K03194791"
BRD_BARICITINIB  <- "BRD-K53581288"

message("Data loaded.")

# =============================================================================
# SECTION 2 · PBMC OBJECT FOR S1  (standalone; replaces the targets dependency)
# =============================================================================
# Builds the annotated PBMC SingleCellExperiment from the public GEO data using
# the project's own helper functions, then caches it. Only needs to run once.
#
# IMPORTANT - THIS IS NOT BIT-REPRODUCIBLE. run_qc_and_clustering() does not
# set a seed, and both UMAP and Louvain clustering are stochastic. A rebuild
# gives a visually similar but non-identical embedding, and cluster numbering
# may differ. S1A/S1B colour by cell_type (not cluster ID), so the panels
# remain valid, but the figure will not byte-match a previous run. State this
# in the Methods rather than implying exact reproducibility.
#
# COMPUTE - this is the heaviest step in the pipeline (all PBMC cells, all
# cell types, through PCA + Harmony + UMAP). Run it with memory free.
# =============================================================================

PBMC_SCE_PATH <- file.path(PROCESSED_DIR, "pbmc_annotated_sce.rds")
DUY_EXTRACTED <- file.path(RAW_DIR, "duy", "extracted")   # per-sample mtx triplets

if (file.exists(PBMC_SCE_PATH)) {
  message("Loading cached PBMC object for S1...")
  pbmc_sce <- readRDS(PBMC_SCE_PATH)
} else {
  message("Building PBMC object for S1 (slow, runs once)...")
  source(here::here("R", "data_loading.R"))
  source(here::here("R", "clustering.R"))

  set.seed(1)   # does not make the original reproducible, but fixes future runs

  pbmc_sce <- load_duy_samples(
    extracted_dir     = DUY_EXTRACTED,
    tissue_filter     = "PBMC",
    drop_unknown_moca = TRUE
  ) |>
    fix_gene_symbols() |>
    run_qc_and_clustering(min_genes_per_cell = 200, max_pct_mt = 20,
                          resolution = 0.5, n_pcs = 30) |>
    identify_monocytes()

  saveRDS(pbmc_sce, PBMC_SCE_PATH)
  message("  Cached to: ", PBMC_SCE_PATH)
}

# =============================================================================
# SECTION 3 · MAIN FIGURES
# =============================================================================

W_UMAP    <- 4.2;  H_UMAP    <- 4.0
W_ISG     <- 3.5;  H_ISG     <- 4.0
W_MoCA    <- 2.8;  H_MoCA    <- 4.0
W_VOLCANO <- 10.5; H_VOLCANO <- 5.5
W_NCS     <- 5.5;  H_NCS     <- 4.5
W_THP1    <- 7.5;  H_THP1    <- 9.5
W_CELLLINE<- 6.0;  H_CELLLINE<- 6.0
W_LOO     <- 8.0;  H_LOO     <- 4.5
W_PERSIST <- 7.0;  H_PERSIST <- 6.0

# -----------------------------------------------------------------------------
# FIG 1A — Monocyte UMAP
# -----------------------------------------------------------------------------
{
  P1A_PT_SIZE <- 0.4

  umap_df <- as.data.frame(Seurat::Embeddings(seu, "umap")) |>
    setNames(c("UMAP1", "UMAP2")) |>
    dplyr::mutate(subcluster = as.character(seu$sub_res05))

  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = subcluster)) +
    geom_point(size = P1A_PT_SIZE, alpha = 0.5) +
    scale_color_manual(values = CLUSTER_COLORS, name = "Subcluster") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 1)) +
    labs(title = "Monocyte subclusters \u2014 UMAP", x = "UMAP 1", y = "UMAP 2") +
    THEME_BASE +
    theme(panel.grid = element_blank(),
          axis.text  = element_blank(), axis.ticks = element_blank())
}
save_panel(p_umap, "fig1_panelA_umap", W_UMAP, H_UMAP)

# -----------------------------------------------------------------------------
# FIG 1B — ISG panel score per subcluster
# -----------------------------------------------------------------------------
{
  P1B_IFN_PANEL <- c("IFITM3","IFITM2","ISG15","MX1","OAS1","IRF7","STAT1")
  P1B_REF_LINE  <- 1.0

  lc <- Seurat::GetAssayData(seu, layer = "data")
  common_genes <- intersect(P1B_IFN_PANEL, rownames(lc))
  score_df <- data.frame(
    ifn_score  = Matrix::colMeans(lc[common_genes, ]),
    subcluster = as.character(seu$sub_res05)
  )

  p_isg_violin <- ggplot(score_df,
                         aes(x = subcluster, y = ifn_score, fill = subcluster)) +
    geom_violin(alpha = 0.7, color = NA, scale = "width") +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.85) +
    scale_fill_manual(values = CLUSTER_COLORS, guide = "none") +
    geom_hline(yintercept = P1B_REF_LINE, linetype = "dashed",
               color = PAL$disease, linewidth = 0.5) +
    annotate("text", x = Inf, y = P1B_REF_LINE,
             label = " ref = 1.0", hjust = 0, vjust = -0.4,
             size = 2.8, color = PAL$disease, fontface = "italic") +
    labs(title = "ISG panel score per subcluster",
         subtitle = paste0("Genes: ", paste(common_genes, collapse = ", ")),
         x = "Subcluster", y = "Mean log-normalized ISG expression") +
    THEME_BASE
}
save_panel(p_isg_violin, "fig1_panelB_isg_violin", W_ISG, H_ISG)

# -----------------------------------------------------------------------------
# FIG 1C — IFN monocyte % by MoCA group
# -----------------------------------------------------------------------------
{
  P1C_HIGH_PATIENTS <- c(1,2,3,4,5)
  P1C_LOW_PATIENTS  <- c(7,8,9,10)
  P1C_MoCA_COLORS   <- c("High" = PAL$mocaHigh, "Low" = PAL$mocaLow)

  ifn_prop_df <- seu@meta.data |>
    dplyr::filter(!is.na(moca_group)) |>
    dplyr::group_by(patient) |>
    dplyr::summarise(n_total = dplyr::n(),
                     n_ifn   = sum(sub_res05 == "5"),
                     ifn_pct = n_ifn / n_total * 100,
                     .groups = "drop") |>
    dplyr::mutate(moca_grp = dplyr::case_when(
      patient %in% P1C_HIGH_PATIENTS ~ "High",
      patient %in% P1C_LOW_PATIENTS  ~ "Low",
      TRUE ~ "Unknown")) |>
    dplyr::filter(moca_grp != "Unknown")

  w_result <- wilcox.test(ifn_pct ~ moca_grp, data = ifn_prop_df, exact = TRUE)
  pval_str <- paste0("p = ", signif(w_result$p.value, 3))
  cat("Figure 1C (computed live): W =", w_result$statistic,
      " p =", signif(w_result$p.value, 4), "\n")

  p_moca_box <- ggplot(ifn_prop_df,
                       aes(x = moca_grp, y = ifn_pct, fill = moca_grp)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(aes(color = moca_grp), width = 0.08, size = 2.8, alpha = 0.9) +
    geom_text(aes(label = paste0("P", patient)),
              vjust = -0.8, size = 2.3, color = PAL$textMute) +
    annotate("text", x = 1.5, y = max(ifn_prop_df$ifn_pct) * 1.15,
             label = pval_str, size = 3.0, color = PAL$text, fontface = "italic") +
    scale_fill_manual(values  = P1C_MoCA_COLORS, guide = "none") +
    scale_color_manual(values = P1C_MoCA_COLORS, guide = "none") +
    labs(title    = "IFN monocyte frequency by cognitive group",
         subtitle = "% CD14\u207a monocytes in subcluster 5 (n.s., Wilcoxon)",
         x = "MoCA group", y = "IFN-activated monocytes (%)") +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank())
}
save_panel(p_moca_box, "fig1_panelC_moca_box", W_MoCA, H_MoCA)

# -----------------------------------------------------------------------------
# FIG 1D — Volcano
# -----------------------------------------------------------------------------
{
  P1D_FDR_THRESH <- 0.05; P1D_LFC_THRESH <- 1.0
  P1D_N_LABEL_UP <- 20;   P1D_N_LABEL_DN <- 10
  P1D_CAP_Y      <- 100

  volc_df <- as.data.frame(de_res) |>
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
    dplyr::mutate(
      neg_log10_fdr = pmin(-log10(padj), P1D_CAP_Y),
      is_capped     = -log10(padj) > P1D_CAP_Y,
      direction = dplyr::case_when(
        padj < P1D_FDR_THRESH & log2FoldChange >  P1D_LFC_THRESH ~ "UP",
        padj < P1D_FDR_THRESH & log2FoldChange < -P1D_LFC_THRESH ~ "DOWN",
        TRUE ~ "NS"),
      direction = factor(direction, levels = c("UP","DOWN","NS")))

  if (!"gene" %in% colnames(volc_df))
    volc_df <- tibble::rownames_to_column(volc_df, "gene")

  top_up  <- volc_df |> dplyr::filter(direction == "UP") |>
    dplyr::arrange(dplyr::desc(log2FoldChange)) |> head(P1D_N_LABEL_UP)
  top_dn  <- volc_df |> dplyr::filter(direction == "DOWN") |>
    dplyr::arrange(log2FoldChange) |> head(P1D_N_LABEL_DN)
  top_sig <- volc_df |> dplyr::arrange(dplyr::desc(neg_log10_fdr)) |> dplyr::slice(1)
  label_df <- dplyr::bind_rows(top_up, top_dn, top_sig) |>
    dplyr::distinct(gene, .keep_all = TRUE)

  n_up <- sum(volc_df$direction == "UP")
  n_dn <- sum(volc_df$direction == "DOWN")
  max_abs_x <- max(abs(volc_df$log2FoldChange), na.rm = TRUE) * 1.05

  p_volcano <- ggplot(volc_df,
                      aes(x = log2FoldChange, y = neg_log10_fdr, color = direction)) +
    geom_point(data = dplyr::filter(volc_df, direction == "NS"),
               size = 0.7, alpha = 0.25, shape = 16) +
    geom_point(data = dplyr::filter(volc_df, direction != "NS", !is_capped),
               size = 1.5, alpha = 0.85, shape = 16) +
    geom_point(data = dplyr::filter(volc_df, is_capped),
               aes(y = P1D_CAP_Y), size = 2.0, shape = 24,
               fill = PAL$disease, color = PAL$disease, alpha = 0.9) +
    ggrepel::geom_text_repel(data = label_df, aes(label = gene), size = 2.7,
                             max.overlaps = 35, segment.color = PAL$textMute,
                             segment.size = 0.25, box.padding = 0.30) +
    geom_hline(yintercept = -log10(P1D_FDR_THRESH), linetype = "dashed",
               color = PAL$textMute, linewidth = 0.4) +
    geom_vline(xintercept = c(-P1D_LFC_THRESH, P1D_LFC_THRESH),
               linetype = "dashed", color = PAL$textMute, linewidth = 0.35) +
    geom_vline(xintercept = 0, color = PAL$rule, linewidth = 0.35) +
    scale_color_manual(
      values = c("UP" = PAL$disease, "DOWN" = PAL$drug, "NS" = PAL$textMute),
      labels = c("UP" = paste0("UP (n=", n_up, ")"),
                 "DOWN" = paste0("DOWN (n=", n_dn, ")"),
                 "NS" = "Not significant"),
      name = NULL,
      guide = guide_legend(override.aes = list(shape = 16, size = 3))) +
    coord_cartesian(xlim = c(-max_abs_x, max_abs_x), ylim = c(0, P1D_CAP_Y + 5)) +
    annotate("text", x = max_abs_x * 0.65, y = P1D_CAP_Y + 3,
             label = "IFN-activated (UP)", hjust = 0.5, size = 3.2,
             color = PAL$disease, fontface = "bold") +
    annotate("text", x = -max_abs_x * 0.65, y = P1D_CAP_Y + 3,
             label = "Classical resting (DOWN)", hjust = 0.5, size = 3.2,
             color = PAL$drug, fontface = "bold") +
    annotate("text", x = Inf, y = P1D_CAP_Y + 1, label = "\u25b2 capped at 100",
             hjust = 1.05, size = 2.3, color = PAL$textMute, fontface = "italic") +
    labs(x = "Shrunken log\u2082 fold change (apeglm)",
         y = "\u2212log\u2081\u2080 (FDR)",
         title = "Differential expression: IFN-activated vs classical resting",
         subtitle = paste0("DESeq2 | ~ patient + state | 1,231 genes at FDR < 0.05 | ",
                           "n=940 IFN-activated vs n=4,904 classical | 9 patients")) +
    THEME_BASE +
    theme(legend.position = "right")
}
save_panel(p_volcano, "fig1_panelD_volcano", W_VOLCANO, H_VOLCANO)

# -----------------------------------------------------------------------------
# FIG 2A — THP1 NCS distribution
# -----------------------------------------------------------------------------
{
  P2A_BINS      <- 50
  P2A_THRESHOLD <- -0.5

  thp1_all <- result_table |>
    dplyr::filter(cell == "THP1", type == "trt_cp",
                  grepl("^BRD-", pert), !is.na(NCS))
  P2A_RSVS_NCS <- thp1_all$NCS[thp1_all$pert == BRD_ROSUVASTATIN][1]

  fig2a_df <- thp1_all |>
    dplyr::filter(NCS != 0) |>
    dplyr::select(pert, NCS) |>
    dplyr::mutate(connectivity = factor(
      ifelse(NCS < 0, "Negative connectivity", "Positive connectivity"),
      levels = c("Negative connectivity", "Positive connectivity")))

  n_nonzero  <- nrow(fig2a_df)
  n_strong   <- sum(fig2a_df$NCS < P2A_THRESHOLD)
  cat("Figure 2A: non-zero THP1 compounds =", n_nonzero,
      "| strong reversers =", n_strong,
      sprintf("(%.1f%%)\n", 100 * n_strong / n_nonzero))

  x_min <- floor(min(fig2a_df$NCS) * 10) / 10
  x_max <- ceiling(max(fig2a_df$NCS) * 10) / 10
  ymax  <- max(hist(fig2a_df$NCS, breaks = P2A_BINS, plot = FALSE)$counts)
  y_top <- ceiling(ymax / 10) * 10

  p_ncs_hist <- ggplot(fig2a_df, aes(x = NCS, fill = connectivity)) +
    geom_histogram(bins = P2A_BINS, color = "white",
                   linewidth = 0.15, alpha = 0.95) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.40) +
    geom_vline(xintercept = P2A_THRESHOLD, linetype = "dashed",
               color = "grey45", linewidth = 0.45) +
    annotate("point", x = P2A_RSVS_NCS, y = y_top * 0.04,
             shape = 21, size = 3.0, stroke = 0.9,
             color = PAL$accent, fill = "white") +
    scale_fill_manual(values = c("Negative connectivity" = PAL$disease,
                                 "Positive connectivity" = PAL$drug),
                      name = NULL) +
    scale_x_continuous(breaks = seq(floor(x_min), ceiling(x_max), by = 0.5),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = seq(0, y_top, by = 10),
                       expand = expansion(mult = c(0, 0.02)),
                       limits = c(0, y_top)) +
    labs(x = "Normalized connectivity score in THP1 cells",
         y = "Number of compounds",
         title = "LINCS2 connectivity score distribution \u2014 THP1",
         subtitle = paste0(n_nonzero, " compounds with non-zero NCS | ",
                           n_strong, " strong reversers (NCS < -0.5) | ",
                           "open circle = rosuvastatin")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    THEME_BASE +
    theme(legend.position = "top", legend.justification = "center",
          legend.box.just = "center", panel.grid.minor = element_blank()) +
    coord_cartesian(xlim = c(x_min - 0.05, x_max + 0.05), clip = "off")
}
save_panel(p_ncs_hist, "fig2_panelA_thp1_ncs", W_NCS, H_NCS)

# -----------------------------------------------------------------------------
# FIG 2B — Top-30 THP1 reversal hits, coloured by MANUSCRIPT TIER
# (MOA/drug-class colouring removed - see header note 4)
# -----------------------------------------------------------------------------
{
  P2B_N_HITS      <- 30
  P2B_LABEL_SIZE  <- 13
  P2B_SIZE_APPR   <- 5.5
  P2B_SIZE_OTHER  <- 3.2
  # X-axis text sized to match Panel A's APPARENT size once both panels are
  # scaled to the same displayed height in PowerPoint. Panel B's canvas
  # (H_THP1) is taller than Panel A's (H_NCS), so it gets shrunk more when
  # matched to A's height - text here is pre-inflated by that exact ratio.
  P2B_XTITLE_SIZE <- 10  * (H_THP1 / H_NCS)   # ~21
  P2B_XTEXT_SIZE  <- 9.5 * (H_THP1 / H_NCS)   # ~20

  fig2b_df <- compound_named |>
    dplyr::filter(!is.na(thp1_NCS)) |>
    dplyr::arrange(thp1_NCS) |>
    head(P2B_N_HITS) |>
    dplyr::left_join(final_cands |> dplyr::select(pert, tier), by = "pert") |>
    dplyr::mutate(
      drug_name = dplyr::coalesce(pert_iname, pert),
      category = dplyr::case_when(
        !is.na(tier) & tier == "Primary lead" ~ "Primary lead",
        !is.na(tier)                          ~ "FDA-approved",
        TRUE                                  ~ "Not FDA-approved"),
      category = factor(category, levels = names(TIER_COLORS)),
      is_approved = category != "Not FDA-approved") |>
    dplyr::mutate(drug_name = factor(drug_name, levels = rev(drug_name)))

  p_thp1 <- ggplot(fig2b_df, aes(x = thp1_NCS, y = drug_name, color = category)) +
    geom_segment(aes(x = 0, xend = thp1_NCS, y = drug_name, yend = drug_name),
                 linewidth = 0.55, color = PAL$rule) +
    geom_point(aes(size = is_approved), alpha = 0.95) +
    scale_size_manual(values = c("TRUE" = P2B_SIZE_APPR, "FALSE" = P2B_SIZE_OTHER),
                      guide = "none") +
    scale_color_manual(values = TIER_COLORS, name = NULL, drop = FALSE) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = "Normalized connectivity score in THP1 cells", y = NULL,
         title = "Top 30 reversal compounds \u2014 THP1",
         subtitle = paste0("Ranked by THP1-specific NCS | larger points = ",
                           "FDA-approved | approval status individually verified"),
         caption = paste0("Approval status verified compound-by-compound from ",
                          "regulatory sources; no aggregate MOA-class analysis ",
                          "performed (see Methods).")) +
    THEME_BASE +
    theme(legend.position = "right",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", size = P2B_LABEL_SIZE),
          axis.text.x = element_text(size = P2B_XTEXT_SIZE),
          axis.title.x = element_text(size = P2B_XTITLE_SIZE),
          legend.text = element_text(size = P2B_LABEL_SIZE - 2),
          legend.title = element_text(size = P2B_LABEL_SIZE - 1, face = "bold"))
}
save_panel(p_thp1, "fig2_panelB_thp1_hits", W_THP1, H_THP1)

# -----------------------------------------------------------------------------
# Helper: cell-line NCS profile for one compound (Fig 3A/3B, S5)
# -----------------------------------------------------------------------------
cellline_panel <- function(brd_id, drug_label, subtitle_extra = NULL,
                           fdr_thresh = 0.05, highlight = PAL$accent) {
  df <- result_table |>
    dplyr::filter(pert == brd_id, type == "trt_cp") |>
    dplyr::mutate(
      sig       = WTCS_FDR < fdr_thresh,
      direction = dplyr::case_when(NCS < 0 ~ "Reversal", NCS > 0 ~ "Activation",
                                   TRUE ~ "Neutral"),
      cell_label = factor(cell, levels = cell[order(NCS)]))

  thp1_idx <- which(levels(df$cell_label) == "THP1")
  n_lines  <- nrow(df)
  n_rev    <- sum(df$NCS < 0 & df$sig)
  thp1_ncs <- df$NCS[df$cell == "THP1"]

  sub <- paste0(n_lines, " LINCS2 cell lines | significant reversal in ",
                n_rev, " | THP1 NCS = ",
                if (length(thp1_ncs)) round(thp1_ncs[1], 2) else "not profiled")
  if (!is.null(subtitle_extra)) sub <- paste0(sub, "\n", subtitle_extra)

  p <- ggplot(df, aes(x = NCS, y = cell_label, fill = direction)) +
    geom_col(width = 0.62, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig, formatC(WTCS_FDR, digits = 2, format = "e"), ""),
                  hjust = ifelse(NCS < 0, 1.1, -0.1)),
              size = 2.4, color = PAL$textMute) +
    scale_fill_manual(values = c("Reversal" = PAL$disease,
                                 "Activation" = PAL$drug,
                                 "Neutral" = PAL$textMute), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.22, 0.25))) +
    coord_cartesian(clip = "off") +
    labs(x = "NCS (Normalized Connectivity Score)", y = NULL,
         title = paste0(drug_label, " \u2014 NCS across LINCS2 cell lines"),
         subtitle = sub,
         caption = paste0(brd_id, " = ", tolower(drug_label))) +
    THEME_BASE +
    theme(legend.position = "top", panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 30, b = 5.5, l = 15))

  if (length(thp1_idx) == 1)
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = thp1_idx - 0.5, ymax = thp1_idx + 0.5,
                      fill = highlight, alpha = 0.10)
  p
}

# -----------------------------------------------------------------------------
# FIG 3A / 3B — Lead compound cell-line profiles
# -----------------------------------------------------------------------------
p_rosuva <- cellline_panel(BRD_ROSUVASTATIN, "Rosuvastatin",
                           "Gold highlight = THP1 (only monocytic line)")
p_roflu  <- cellline_panel(BRD_ROFLUMILAST,  "Roflumilast",
                           "Gold highlight = THP1 (only monocytic line)")

save_panel(p_rosuva, "fig3_panelA_rosuvastatin_celllines", W_CELLLINE, H_CELLLINE)
save_panel(p_roflu,  "fig3_panelB_roflumilast_celllines",  W_CELLLINE, H_CELLLINE)

# -----------------------------------------------------------------------------
# FIG 4A — Leave-one-out signature stability  (promoted from supplementary)
# -----------------------------------------------------------------------------
{
  P4A_HIGH <- c(1,2,3,4,5); P4A_LOW <- c(7,8,9,10)
  P4A_THRESHOLD <- 0.70
  P4A_SHAPES <- c("UP genes" = 16, "DOWN genes" = 17)
  P4A_YLIM   <- c(0.58, 0.86)

  patient_order <- loo_df |>
    dplyr::arrange(dplyr::desc(jaccard_mean)) |>
    dplyr::mutate(patient_label = paste0("P", dropped)) |>
    dplyr::pull(patient_label)

  p4a_df <- loo_df |>
    dplyr::mutate(
      patient_label = factor(paste0("P", dropped), levels = patient_order),
      moca_group = dplyr::case_when(dropped %in% P4A_HIGH ~ "High MOCA",
                                    dropped %in% P4A_LOW  ~ "Low MOCA",
                                    TRUE ~ "Unknown")) |>
    tidyr::pivot_longer(c(jaccard_up, jaccard_dn),
                        names_to = "set", values_to = "jaccard") |>
    dplyr::mutate(set = dplyr::recode(set, jaccard_up = "UP genes",
                                      jaccard_dn = "DOWN genes"))

  p4a_text <- loo_df |>
    dplyr::mutate(patient_label = factor(paste0("P", dropped), levels = patient_order),
                  label = formatC(jaccard_mean, digits = 3, format = "f"))

  p_loo <- ggplot(p4a_df, aes(x = patient_label, y = jaccard,
                              color = moca_group, shape = set)) +
    geom_hline(yintercept = P4A_THRESHOLD, linetype = "dashed",
               color = PAL$disease, linewidth = 0.75) +
    annotate("text", x = Inf, y = P4A_THRESHOLD + 0.013,
             label = paste0("Threshold = ", P4A_THRESHOLD),
             hjust = 1.05, size = 2.8, color = PAL$disease, fontface = "italic") +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text(data = p4a_text,
              aes(x = patient_label, y = P4A_YLIM[1] + 0.004, label = label),
              color = PAL$textMute, size = 2.6, inherit.aes = FALSE) +
    scale_color_manual(values = c("High MOCA" = PAL$mocaHigh,
                                  "Low MOCA"  = PAL$mocaLow),
                       labels = c("High MOCA" = "High MoCA (23\u201329)",
                                  "Low MOCA"  = "Low MoCA (19\u201320)"),
                       name = "MoCA group (excluded patient)") +
    scale_shape_manual(values = P4A_SHAPES, name = "Gene set") +
    scale_y_continuous(limits = P4A_YLIM,
                       breaks = seq(P4A_YLIM[1], P4A_YLIM[2], 0.05)) +
    labs(x = "Patient removed from analysis",
         y = "Jaccard similarity to full 9-patient signature",
         title = "Leave-one-out signature stability",
         subtitle = paste0("Overall mean: UP = 0.714 | DOWN = 0.733 | ",
                           "Combined = 0.724 | Threshold = 0.70")) +
    THEME_BASE +
    theme(legend.position = "right",
          axis.text.x = element_text(size = BASE_SIZE - 1))
}
save_panel(p_loo, "fig4_panelA_loo_stability", W_LOO, H_LOO)

# -----------------------------------------------------------------------------
# FIG 4B — Per-gene persistence vs DE significance
# -----------------------------------------------------------------------------
{
  P4B_TIER_COLORS <- c("Core stable (9/9)" = PAL$good,
                       "High (7\u20138/9)"  = "#8CB89A",
                       "Variable (\u22646/9)" = PAL$disease)
  P4B_ISG_LABEL <- c("PARP14","HERC5","OAS3","IFIT3","ISG15","IFIT1",
                     "OAS1","RSAD2","MX1","IFIT2","SIGLEC1","CXCL10")

  de_for_fdr <- as.data.frame(de_res)
  if (!"gene" %in% colnames(de_for_fdr)) de_for_fdr$gene <- rownames(de_for_fdr)

  p4b_df <- persist_df |>
    dplyr::filter(in_sig) |>
    dplyr::left_join(dplyr::select(de_for_fdr, gene, padj), by = "gene") |>
    dplyr::mutate(
      neglog10fdr = -log10(padj),
      tier = dplyr::case_when(
        n_iters == 9        ~ "Core stable (9/9)",
        n_iters %in% c(7,8) ~ "High (7\u20138/9)",
        TRUE                ~ "Variable (\u22646/9)"),
      tier = factor(tier, levels = names(P4B_TIER_COLORS)),
      label = ifelse(gene %in% P4B_ISG_LABEL, gene, NA_character_))

  n_core <- sum(p4b_df$tier == "Core stable (9/9)")
  n_high <- sum(p4b_df$tier == "High (7\u20138/9)")
  n_var  <- sum(p4b_df$tier == "Variable (\u22646/9)")

  p_persist <- ggplot(p4b_df, aes(x = tier, y = neglog10fdr, fill = tier)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = PAL$textMute, linewidth = 0.5) +
    annotate("text", x = Inf, y = -log10(0.05), label = "FDR = 0.05",
             hjust = 1.05, vjust = -0.5, size = 2.8,
             color = PAL$textMute, fontface = "italic") +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.35) +
    geom_jitter(aes(color = tier), width = 0.12, height = 0,
                size = 2.0, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.6, color = PAL$text,
                             max.overlaps = 20, min.segment.length = 0,
                             segment.size = 0.25, box.padding = 0.3, na.rm = TRUE) +
    scale_fill_manual(values = P4B_TIER_COLORS, guide = "none") +
    scale_color_manual(values = P4B_TIER_COLORS, guide = "none") +
    labs(x = "Leave-one-out stability tier",
         y = expression(bold(-log[10]~"(FDR), full DE analysis")),
         title = "LOO per-gene persistence versus DE significance",
         subtitle = paste0("137 genes entering the signature in ≥1 iteration | ",
                           "core 9/9 n = ", n_core, "; high 7–8/9 n = ", n_high,
                           "; variable ≤6/9 n = ", n_var)) +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(color = "black", face = "bold"))
}
save_panel(p_persist, "fig4_panelB_gene_persistence", W_PERSIST, H_PERSIST)

# =============================================================================
# SECTION 4 · SUPPLEMENTARY FIGURES
# =============================================================================

# -----------------------------------------------------------------------------
# S1 — PBMC annotation (A UMAP, B lineage dotplot, C subcluster composition)
# -----------------------------------------------------------------------------
{
  S1_CELLTYPE_COL <- "cell_type"
  S1_CELL_COLORS  <- c("monocyte" = PAL$disease, "other" = PAL$neutral)
  S1_MARKERS <- c("CD14","LYZ","S100A8","S100A9","VCAN","FCN1",
                  "CD3D","CD3E","CD19","MS4A1",
                  "GNLY","NKG7","FCER1A","CST3","FCGR3A","MS4A7")
  S1_PT_SIZE <- 0.25
  W_S1 <- 13; H_S1 <- 9

  # ── A: PBMC UMAP by cell type ──
  umap_coords <- SingleCellExperiment::reducedDim(pbmc_sce, "UMAP")
  cell_types  <- as.character(
    SummarizedExperiment::colData(pbmc_sce)[[S1_CELLTYPE_COL]])

  umap_all <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2],
                         cell_type = cell_types, stringsAsFactors = FALSE)
  rm(umap_coords, cell_types); gc()

  pS1a <- ggplot(umap_all, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
    geom_point(size = S1_PT_SIZE, alpha = 0.4) +
    scale_color_manual(values = S1_CELL_COLORS, name = "Cell type",
                       na.value = PAL$neutral) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "PBMC clusters \u2014 UMAP", x = "UMAP 1", y = "UMAP 2") +
    THEME_BASE +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank())

  # ── B: lineage marker dotplot ──
  markers_present <- intersect(S1_MARKERS, rownames(pbmc_sce))
  sce_subset <- pbmc_sce[markers_present, ]
  seu_all <- Seurat::CreateSeuratObject(
    counts    = SingleCellExperiment::counts(sce_subset),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce_subset))
  ) |> Seurat::NormalizeData(verbose = FALSE)
  rm(sce_subset); gc()

  pS1b <- Seurat::DotPlot(seu_all, features = markers_present,
                          group.by = S1_CELLTYPE_COL) +
    scale_color_gradient2(low = PAL$drug, mid = "white", high = PAL$disease,
                          midpoint = 0, name = "Avg.\nexpression") +
    labs(title = "Canonical lineage marker expression", x = NULL, y = NULL) +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2))
  rm(seu_all); gc()

  # ── C: subcluster composition per patient ──
  S1C_HIGH <- c(1,2,3,4,5); S1C_LOW <- c(7,8,9,10)
  comp_df <- cell_meta |>
    dplyr::group_by(patient, subcluster) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      moca = dplyr::case_when(patient %in% S1C_HIGH ~ "High",
                              patient %in% S1C_LOW  ~ "Low", TRUE ~ "NA"),
      patient_label = factor(
        paste0("P", patient, " (", moca, ")"),
        levels = paste0("P", c(S1C_HIGH, S1C_LOW), " (",
                        c(rep("High", length(S1C_HIGH)),
                          rep("Low", length(S1C_LOW))), ")")),
      subcluster = factor(subcluster, levels = as.character(0:6)))

  pS1c <- ggplot(comp_df, aes(x = patient_label, y = prop, fill = subcluster)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = CLUSTER_COLORS, name = "Subcluster") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(title = "Monocyte subcluster composition per patient",
         x = NULL, y = "Proportion of monocytes") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2),
          panel.grid.major.x = element_blank())

  save_panel(pS1a, "suppfig_S1_panelA_pbmc_umap",       5.5, 5)
  save_panel(pS1b, "suppfig_S1_panelB_marker_dotplot",  7.5, 5)
  save_panel(pS1c, "suppfig_S1_panelC_subcluster_comp", 9,   4.5)

  figS1 <- (pS1a + pS1b + patchwork::plot_layout(widths = c(1, 1.35))) / pS1c +
    patchwork::plot_layout(heights = c(1, 0.9)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_combined(figS1, "suppfig_S1_combined", W_S1, H_S1)
  message("Supplementary Figure S1 saved.")
}

# -----------------------------------------------------------------------------
# S2 — Pseudobulk heatmap (counts now from the PUBLISHED Seurat object)
# -----------------------------------------------------------------------------
{
  S2_N_UP <- 30; S2_N_DOWN <- 30
  W_S2 <- 11; H_S2 <- 10

  HEAT_COLORS <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
  ANN_COLORS <- list(
    State      = c("IFN" = "#5c0099", "classical" = "#d9d9d9"),
    MoCA_group = c("High (23\u201329)" = "#7a8c6a", "Low (19\u201320)" = "#b04a36"),
    Direction  = c("UP" = "#e07b3c", "DOWN" = "#1f78b4"))

  counts_mat <- Seurat::GetAssayData(seu, layer = "counts")

  cm <- cell_meta |>
    dplyr::filter(subcluster %in% c("0", "5")) |>
    dplyr::mutate(state = ifelse(subcluster == "5", "IFN", "classical"),
                  col_name = paste0("P", patient, "_", state))

  pb_both <- sapply(unique(cm$col_name), function(g) {
    Matrix::rowSums(counts_mat[, cm$cell[cm$col_name == g], drop = FALSE])
  })

  pb_meta <- cm |>
    dplyr::select(col_name, patient, state, moca_group) |>
    dplyr::distinct() |>
    dplyr::arrange(state, patient) |>
    dplyr::mutate(moca_group = dplyr::case_when(
      moca_group %in% c("High", "high") ~ "High (23\u201329)",
      moca_group %in% c("Low",  "low")  ~ "Low (19\u201320)",
      TRUE ~ as.character(moca_group)))

  col_order  <- pb_meta$col_name
  pb_ordered <- pb_both[, col_order]
  pb_logcpm  <- log1p(sweep(pb_ordered, 2, Matrix::colSums(pb_ordered) / 1e6, "/"))

  sig_up      <- head(inph_signature$up,   S2_N_UP)
  sig_down    <- head(inph_signature$down, S2_N_DOWN)
  sig_present <- intersect(c(sig_up, sig_down), rownames(pb_logcpm))

  heat_scaled <- t(scale(t(as.matrix(pb_logcpm[sig_present, ]))))

  ann_col <- data.frame(State = pb_meta$state, MoCA_group = pb_meta$moca_group,
                        row.names = col_order)
  ann_row <- data.frame(Direction = ifelse(sig_present %in% sig_up, "UP", "DOWN"),
                        row.names = sig_present)
  stopifnot(all(unique(ann_col$MoCA_group) %in% names(ANN_COLORS$MoCA_group)))

  gap_col <- sum(pb_meta$state == "classical")

  for (dest in c(DIR_TITLED, DIR_CLEAN)) {
    for (ext in c("pdf", "png")) {
      out_path <- file.path(dest, paste0("suppfig_S2_heatmap.", ext))
      if (ext == "pdf") pdf(out_path, width = W_S2, height = H_S2)
      else png(out_path, width = W_S2 * COMBINED_DPI, height = H_S2 * COMBINED_DPI,
               res = COMBINED_DPI)
      pheatmap::pheatmap(
        heat_scaled, color = HEAT_COLORS,
        annotation_col = ann_col, annotation_row = ann_row,
        annotation_colors = ANN_COLORS,
        cluster_rows = TRUE, cluster_cols = TRUE,
        cutree_cols = 2, cutree_rows = 2,
        show_rownames = TRUE, show_colnames = FALSE,
        fontsize_row = 6.5, fontsize = 9, border_color = NA, gaps_col = gap_col,
        main = if (dest == DIR_TITLED) paste0(
          "Top ", S2_N_UP, " UP + ", S2_N_DOWN,
          " DOWN signature genes across 18 pseudobulk samples\n",
          "Row-scaled log-CPM | 9 iNPH patients \u00d7 2 cell states") else "",
        legend_breaks = c(-2, -1, 0, 1, 2),
        legend_labels = c("-2", "-1", "0", "1", "2\n(scaled\nexpression)"))
      dev.off()
    }
  }
  message("Supplementary Figure S2 saved (titled + clean).")
}

# -----------------------------------------------------------------------------
# S3 — Bulk vs subcluster scoring
# -----------------------------------------------------------------------------
{
  Seurat::DefaultAssay(seu) <- "RNA"
  expr    <- Seurat::GetAssayData(seu, layer = "data")
  up_in   <- intersect(rownames(expr), inph_signature$up)
  down_in <- intersect(rownames(expr), inph_signature$down)
  score   <- Matrix::colMeans(expr[up_in, , drop = FALSE]) -
    Matrix::colMeans(expr[down_in, , drop = FALSE])

  s3_df <- rbind(
    data.frame(score = score, group = "All iNPH monocytes\n(22,197 cells)"),
    data.frame(score = score[seu$sub_res05 == "5"], group = "IFN-activated\n(940 cells)"))
  s3_df$group <- factor(s3_df$group,
    levels = c("All iNPH monocytes\n(22,197 cells)", "IFN-activated\n(940 cells)"))

  bulk_mean <- round(mean(score), 3)
  ifn_mean  <- round(mean(score[seu$sub_res05 == "5"]), 3)
  cat("S3: bulk mean =", bulk_mean, "| IFN-activated mean =", ifn_mean, "\n")

  pS3 <- ggplot(s3_df, aes(x = group, y = score, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, color = "grey40", linewidth = 0.3) +
    geom_boxplot(width = 0.16, outlier.shape = NA, color = "grey25",
                 fill = "white", alpha = 0.85, linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    scale_fill_manual(values = c(PAL$resting, PAL$ifn), guide = "none") +
    labs(x = NULL, y = "iNPH IFN signature score\n(mean UP \u2212 mean DOWN)",
         title = "Bulk vs IFN-subcluster scoring",
         subtitle = paste0("Bulk mean = ", bulk_mean,
                           " (resting majority masks the signal); ",
                           "IFN-activated mean = ", ifn_mean)) +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank())

  save_panel(pS3, "suppfig_S3_bulk_vs_subcluster", 5.5, 4.5)
  message("Supplementary Figure S3 saved.")
}

{
  # ============================================================
  # TUNABLE KNOBS
  # ============================================================
  S4_TEXT_SCALE  <- 1.15
  S4_ROW_HEIGHTS <- c(1, 0.65)
  S4_BOX_WIDTH   <- 0.6
  W_S4_COMBINED  <- 12
  H_S4_COMBINED  <- 9

  # REMOVED: S4_LEGEND_C_X / S4_LEGEND_C_Y and legend.position.inside.
  # "Inside" legends are computed as a fraction of panel space, so they
  # break the moment the row gets squeezed by patchwork::plot_layout()
  # in the combined figure — that's exactly why labels kept escaping.
  # A top-outside legend has a fixed, predictable slot that can never
  # overlap the bar or run outside the panel border.

  s4_text_theme <- theme(
    plot.title    = element_text(size = (BASE_SIZE + 1)   * S4_TEXT_SCALE),
    plot.subtitle = element_text(size = (BASE_SIZE - 1.5) * S4_TEXT_SCALE),
    plot.caption  = element_text(size = (BASE_SIZE - 2.5) * S4_TEXT_SCALE),
    axis.text     = element_text(size = (BASE_SIZE - 1.5) * S4_TEXT_SCALE),
    axis.title    = element_text(size = (BASE_SIZE - 1)   * S4_TEXT_SCALE),
    legend.text   = element_text(size = (BASE_SIZE - 2.5) * S4_TEXT_SCALE),
    legend.title  = element_text(size = (BASE_SIZE - 1)   * S4_TEXT_SCALE)
  )

  S4_COLORS <- c("Healthy control" = PAL$neutral,
                 "SLE"             = PAL$drug,
                 "iNPH"            = PAL$disease)

  # ── A & B: unchanged from last version ──
  pop_df <- dplyr::bind_rows(
    BENCH_A$pop_tables$hc_hi   |> dplyr::mutate(group = "Healthy control"),
    BENCH_A$pop_tables$sle_hi  |> dplyr::mutate(group = "SLE"),
    BENCH_A$pop_tables$inph_c5 |> dplyr::mutate(group = "iNPH")) |>
    dplyr::mutate(group = factor(group, levels = names(S4_COLORS)))

  pS4a <- ggplot(pop_df, aes(x = group, y = mean_sig, fill = group)) +
    geom_boxplot(width = S4_BOX_WIDTH, outlier.shape = NA, alpha = 0.55, color = "grey25") +
    geom_jitter(width = 0.10, size = 1.4, alpha = 0.5, color = PAL$textSoft) +
    scale_fill_manual(values = S4_COLORS, guide = "none") +
    labs(x = NULL, y = "AUCell signature score (donor mean)",
         title = "Population-level IFN signature score",
         subtitle = paste0(
           "Donor-level, unfiltered cluster 5 | medians: iNPH ",
           round(BENCH_A$med_pop["inph_c5"], 4), ", SLE-high ",
           round(BENCH_A$med_pop["sle_hi"], 4), ", HC-high ",
           round(BENCH_A$med_pop["hc_hi"], 4), "\n",
           "iNPH vs SLE-high p = ", signif(BENCH_A$p_c5_sle, 3),
           " | iNPH vs HC-high p = ", signif(BENCH_A$p_c5_hc, 3)),
         caption = "Confirmed identical across 4 independent runs (150-cap, uncapped, v3, v4).") +
    THEME_BASE + s4_text_theme +
    theme(panel.grid.major.x = element_blank())

  burden_df <- dplyr::bind_rows(
    BENCH_A$burden_tables$hc   |> dplyr::mutate(group = "Healthy control"),
    BENCH_A$burden_tables$sle  |> dplyr::mutate(group = "SLE"),
    BENCH_A$burden_tables$inph |> dplyr::mutate(group = "iNPH")) |>
    dplyr::mutate(group = factor(group, levels = names(S4_COLORS)))

  pS4b <- ggplot(burden_df, aes(x = group, y = burden_pct, fill = group)) +
    geom_boxplot(width = S4_BOX_WIDTH, outlier.shape = NA, alpha = 0.55, color = "grey25") +
    geom_jitter(width = 0.10, size = 1.4, alpha = 0.5, color = PAL$textSoft) +
    scale_fill_manual(values = S4_COLORS, guide = "none") +
    labs(x = NULL, y = "IFN-high monocytes per donor (%)",
         title = "IFN-high monocyte burden",
         subtitle = paste0(
           "Donor-level, shared healthy-anchored threshold | medians: HC ",
           round(BENCH_A$burden_med["hc"], 1), "%, SLE ",
           round(BENCH_A$burden_med["sle"], 1), "%, iNPH ",
           round(BENCH_A$burden_med["inph"], 1), "%\n",
           "iNPH vs SLE p = ", signif(BENCH_A$p_burden_sle, 3),
           " | iNPH vs HC p = ", signif(BENCH_A$p_burden_hc, 3),
           " (HC comparison is scoring-method-dependent — see Methods)"),
         caption = paste0(
           "Exploratory: the iNPH comparison crosses the Duy/Perez dataset ",
           "boundary and is batch-sensitive.")) +
    THEME_BASE + s4_text_theme +
    theme(panel.grid.major.x = element_blank())

  # ── C: cluster-5 / ISG-high overlap ──
  ov <- BENCH_A$overlap
  overlap_df <- data.frame(
    category = factor(c("Above shared IFN-high threshold", "Below shared IFN-high threshold"),
                      levels = c("Above shared IFN-high threshold", "Below shared IFN-high threshold")),
    n = c(ov[["n_both"]], ov[["n_c5"]] - ov[["n_both"]]))

  pS4c <- ggplot(overlap_df, aes(x = "", y = n, fill = category)) +
    # FIX for the gap: bumped width from 0.55 to 0.92 so the bar fills
    # almost the entire row height. On a single-category flipped axis,
    # geom_col width IS the thing that creates the whitespace band you
    # were seeing above/below the bar.
    geom_col(width = 0.92) +
    # Trim the residual discrete-axis padding further (default reserves
    # ~0.6 units of category space on each side even at width < 1).
    scale_x_discrete(expand = expansion(add = 0.05)) +
    # Keeps the bar flush against x = 0, small buffer on the far end.
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_fill_manual(values = c("Above shared IFN-high threshold" = PAL$disease,
                                 "Below shared IFN-high threshold" = PAL$neutral), name = NULL) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Cluster 5 cells (n = 940)",
         title = "Cluster 5 corresponds to the IFN-high population",
         subtitle = paste0(
           sprintf("%.1f%%", ov[["pct_c5_is_ifnhi"]]),
           " of cluster-5 cells clear the shared IFN-high threshold ",
           "(AUCell; ", sprintf("%.1f%%", BENCH_B$overlap[["pct_c5_is_ifnhi"]]),
           " under mean-difference)\n",
           "Population-level score: iNPH ",
           round(BENCH_A$med_pop["inph_c5"], 4), " vs SLE-high ",
           round(BENCH_A$med_pop["sle_hi"], 4),
           " (p = ", signif(BENCH_A$p_c5_sle, 3), ")"),
         caption = paste0(
           "Unsupervised clustering and marker-panel definition identify ",
           "essentially the same cells.")) +
    THEME_BASE + s4_text_theme +
    theme(
      legend.position  = "top",       # or "right"/"left" if you prefer
      legend.direction = "vertical",
      legend.justification = "right",
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # <- ncol = 1, not nrow = 1

  save_panel(pS4a, "suppfig_S4_panelA_population", 4, 4)
  save_panel(pS4b, "suppfig_S4_panelB_burden",  4, 4)
  save_panel(pS4c, "suppfig_S4_panelC_overlap", 9, 2)

  figS4 <- (pS4a + pS4b + patchwork::plot_layout(widths = c(1, 1))) / pS4c +
    patchwork::plot_layout(heights = S4_ROW_HEIGHTS) +
    patchwork::plot_annotation(tag_levels = "A")
  save_combined(figS4, "suppfig_S4_combined", W_S4_COMBINED, H_S4_COMBINED)
  message("Supplementary Figure S4 saved.")
}

# -----------------------------------------------------------------------------
# S5 — Baricitinib cell-line profile  (replaces the retired JAK2 panel)
# -----------------------------------------------------------------------------
{
  bari_rank <- compound_named |>
    dplyr::filter(!is.na(thp1_NCS)) |>
    dplyr::arrange(thp1_NCS) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::filter(pert == BRD_BARICITINIB) |>
    dplyr::pull(rank)

  pS5 <- cellline_panel(
    BRD_BARICITINIB, "Baricitinib",
    subtitle_extra = paste0(
      "Post hoc, mechanistically aligned secondary candidate. Ranked ",
      bari_rank, " by THP1 NCS \u2014 identified outside the top-30 review ",
      "by JAK/IFN-pathway name matching, not by rank."),
    highlight = PAL$purple)

  save_panel(pS5, "suppfig_S5_baricitinib_celllines", W_CELLLINE, H_CELLLINE)
  cat("S5: baricitinib rank =", bari_rank, "\n")
  message("Supplementary Figure S5 saved.")
}

# -----------------------------------------------------------------------------
# S6 — Proposed mechanism schematic  (JAK2 inhibitor box REMOVED)
# -----------------------------------------------------------------------------
{
  FS_DRUG  <- 3.8;  FS_DRUG_S <- 2.7
  FS_TGT   <- 3.2;  FS_INTER  <- 3.0
  FS_CONV  <- 3.8;  FS_CONV_S <- 2.5
  FS_OUT   <- 4.0;  FS_OUT_S  <- 2.6
  FS_NOTE  <- 3
  W_MECH   <- 7;    H_MECH    <- 9

  AR <- arrow(length = unit(0.18, "cm"), type = "closed")

  pS6 <- ggplot() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(4, 4, 4, 4),
          aspect.ratio    = 11.2 / 10) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 11.2), expand = FALSE) +

    annotate("text", x = 5, y = 10.9,
             label = paste0("Proposed mechanism: rosuvastatin and roflumilast\n",
                            "may oppose IFN-associated signalling"),
             size = FS_DRUG, fontface = "bold.italic",
             color = PAL$text, hjust = 0.5) +

    # ── ROSUVASTATIN (LEFT) ──
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 9.0, ymax = 10,
             fill = PAL$accent, color = NA) +
    annotate("text", x = 2.2, y = 9.72, label = "Rosuvastatin",
             size = FS_DRUG, fontface = "bold", color = "white", hjust = 0.5) +
    annotate("text", x = 2.2, y = 9.3,
             label = "HMGCR inhibitor  \u00b7  FDA 2003",
             size = FS_DRUG_S, color = "#fef3c7", hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 2.2, y = 9.0, yend = 8.77,
             color = PAL$inhibit, linewidth = 0.9) +
    annotate("segment", x = 1.65, xend = 2.75, y = 8.77, yend = 8.77,
             color = PAL$inhibit, linewidth = 2.0) +
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 7.9, ymax = 8.65,
             fill = PAL$inhibitFill, color = PAL$inhibit, linewidth = 0.6) +
    annotate("text", x = 2.2, y = 8.275, label = "HMG-CoA reductase (HMGCR)",
             size = FS_TGT, fontface = "bold", color = PAL$inhibit, hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 2.2, y = 7.9, yend = 7.42,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 6.6, ymax = 7.35,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 2.2, y = 6.975, label = "Isoprenoid synthesis \u2193",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 2.2, y = 6.6, yend = 6.12,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 5.2, ymax = 6.05,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 2.2, y = 5.77,
             label = "RAC1/RHOA membrane anchoring \u2193",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("text", x = 2.2, y = 5.43, label = "JAK activation \u2193",
             size = FS_INTER - 0.2, color = PAL$textMute, hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 4.05, y = 5.2, yend = 4.55,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +

    # ── ROFLUMILAST (RIGHT) ──
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 9.0, ymax = 10.0,
             fill = PAL$drug, color = NA) +
    annotate("text", x = 7.8, y = 9.72, label = "Roflumilast",
             size = FS_DRUG, fontface = "bold", color = "white", hjust = 0.5) +
    annotate("text", x = 7.8, y = 9.3,
             label = "PDE4 inhibitor  \u00b7  FDA 2011",
             size = FS_DRUG_S, color = "#cffafe", hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 7.8, y = 9.0, yend = 8.77,
             color = PAL$inhibit, linewidth = 0.9) +
    annotate("segment", x = 7.25, xend = 8.35, y = 8.77, yend = 8.77,
             color = PAL$inhibit, linewidth = 2.0) +
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 7.9, ymax = 8.65,
             fill = PAL$inhibitFill, color = PAL$inhibit, linewidth = 0.6) +
    annotate("text", x = 7.8, y = 8.275, label = "Phosphodiesterase 4 (PDE4)",
             size = FS_TGT, fontface = "bold", color = PAL$inhibit, hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 7.8, y = 7.9, yend = 7.42,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 6.6, ymax = 7.35,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 7.8, y = 6.975, label = "Intracellular cAMP \u2191",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 7.8, y = 6.6, yend = 6.12,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 5.2, ymax = 6.05,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 7.8, y = 5.77, label = "PKA activation \u2191",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("text", x = 7.8, y = 5.43,
             label = "NF-\u03baB \u2193  \u00b7  STAT1 phosphorylation \u2193",
             size = FS_INTER - 0.2, color = PAL$textMute, hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 5.95, y = 5.2, yend = 4.55,
             color = PAL$textSoft, linewidth = 0.7, arrow = AR) +

    # ── CONVERGENCE ──
    annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 3.1, ymax = 4.5,
             fill = PAL$navy, color = NA) +
    annotate("text", x = 5.0, y = 4.1, label = "JAK / STAT signalling",
             size = FS_CONV, fontface = "bold", color = "white", hjust = 0.5) +
    annotate("text", x = 5.0, y = 3.6,
             label = "STAT1 \u2193  \u00b7  NF-\u03baB \u2193  \u00b7  IRF \u2193",
             size = FS_CONV_S, color = PAL$textMute, hjust = 0.5) +
    annotate("text", x = 5.0, y = 3.22, label = "(predicted shared endpoint)",
             size = FS_CONV_S - 0.3, color = PAL$textMute,
             fontface = "italic", hjust = 0.5) +
    annotate("segment", x = 5.0, xend = 5.0, y = 3.1, yend = 2.62,
             color = PAL$textSoft, linewidth = 1.0, arrow = AR) +

    # ── OUTCOME ──
    annotate("rect", xmin = 2.2, xmax = 7.8, ymin = 1.5, ymax = 2.6,
             fill = PAL$good, color = NA) +
    annotate("text", x = 5.0, y = 2.3, label = "ISG transcription \u2193",
             size = FS_OUT, fontface = "bold", color = "white", hjust = 0.5) +
    annotate("text", x = 5.0, y = 1.9,
             label = "IFIT1  \u00b7  ISG15  \u00b7  MX1  \u00b7  OAS1  \u00b7  CXCL10  \u00b7  SIGLEC1",
             size = FS_OUT_S, color = "#86efac", hjust = 0.5) +

    annotate("text", x = 5.0, y = 1.05,
             label = paste0("Bars (\u22a3) = enzyme inhibition  \u00b7  ",
                            "Arrows (\u2192) = downstream activation  \u00b7  ",
                            "\u2191/\u2193 = relative change"),
             size = FS_NOTE, color = PAL$textMute,
             fontface = "italic", hjust = 0.5)

  save_panel(pS6, "suppfig_S6_mechanism", W_MECH, H_MECH)
  message("Supplementary Figure S6 saved.")
}

# =============================================================================
# SECTION 5 · COMBINED MAIN FIGURES
# =============================================================================
message("Assembling combined main figures...")

fig1 <- (p_umap + p_isg_violin + p_moca_box +
           patchwork::plot_layout(widths = c(1.4, 1.2, 1))) / p_volcano +
  patchwork::plot_layout(heights = c(1, 1.3)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig1, "fig1_combined", 12, 9.5)

fig2 <- p_ncs_hist / p_thp1 +
  patchwork::plot_layout(heights = c(1, 1.8)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig2, "fig2_combined", 8, 12)

fig3 <- p_rosuva + p_roflu +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig3, "fig3_combined", 12, 6)

fig4 <- p_loo / p_persist +
  patchwork::plot_layout(heights = c(1, 1.4)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig4, "fig4_combined", 11, 10)

message("=== All figures complete (main 1-4, supplementary S1-S6) ===")
message("Titled versions : ", DIR_TITLED)
message("Clean versions  : ", DIR_CLEAN)
