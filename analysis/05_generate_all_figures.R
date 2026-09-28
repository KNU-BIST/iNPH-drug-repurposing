# =============================================================================
# Final manuscript / supplementary figure exports
# Build: FIGURES_FINAL_20260912
# Install: analysis/05_generate_all_figures_FINAL.R
# Based on the user's 05_generate_all_figures_revised(7).R, with 05F merged inline.
# Run from the project root. No other figure-generation script is sourced.
#
# MAIN: 1 state; 2 HC/SLE context; 3 THP1 screen; 4 candidate robustness/nulls;
#       5 candidate cell-line profiles.
# SUPPLEMENT: S1 PBMC annotation; S2 state/P6 sensitivity; S3 pseudobulk heatmap;
#       S4 population-vs-subcluster scoring; S5 conditional signature LOO;
#       S6 baricitinib; S7 CSF context; S8 hypothetical pathway (LAST).
#
# Each category has separate clean/ and with_titles/ folders under
# figures/revision/. Final order/filename manifests are generated.
# Main Figure 4: C=roflumilast, D=rosuvastatin, matching the latest user assembly.
# The user's final Figure 1E typography, plain headers and math axis titles remain.
# S1 uses the existing validated PBMC cache; no raw-data reconstruction, clustering,
# DE fitting, or LINCS querying is performed. Saved-data plot summaries are retained.
# =============================================================================

# =============================================================================
# SECTION 0 · CONFIGURATION AND PREFLIGHT
# =============================================================================

required_packages <- c(
  "here",
  "dplyr",
  "ggplot2",
  "tidyr",
  "readr",
  "ggrepel",
  "patchwork",
  "Seurat",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "Matrix",
  "signatureSearch",
  "pheatmap",
  "scales",
  "tibble",
  "grid",
  "digest"
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
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

PROCESSED_DIR <- here::here("data", "processed")
RAW_DIR       <- here::here("data", "raw")
# Final exports use separate main and supplementary folders.
FIGURES_DIR <- here::here("figures", "revision")
OUTPUT_ROOTS <- c(
  main = file.path(FIGURES_DIR, "main"),
  supplementary = file.path(FIGURES_DIR, "supplementary")
)
EXPORT_RECORDS <- list()

dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)


PANEL_DPI <- getOption("inph.figures.dpi", 600)
COMBINED_DPI <- PANEL_DPI
BASE_SIZE    <- 11
SAVE_TITLED  <- isTRUE(getOption("inph.figures.save_titled", TRUE))
SIZE_OVERRIDES <- getOption("inph.figures.sizes", list())
# Example before sourcing: options(inph.figures.sizes = list(
#   fig1_panelA_umap = c(3.2, 3.0), fig1_panelE_volcano = c(8.5, 4.5)))

if (!is.list(SIZE_OVERRIDES)) {
  stop("Option 'inph.figures.sizes' must be a named list of c(width, height) values.")
}

# Prefer Cairo for PDF export so Unicode labels render consistently on Windows.
# The fallback keeps the script functional on systems without Cairo support.
open_pdf_device <- function(filename, width, height, ...) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      filename = filename,
      width = width,
      height = height,
      ...
    )
  } else {
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      ...
    )
  }
}

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
    text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.65),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = grid::unit(2, "mm"),
    plot.title       = element_text(face = "bold", size = BASE_SIZE + 1,
                                    color = PAL$text),
    plot.subtitle    = element_text(size = BASE_SIZE - 1.5,
                                    color = "black"),
    plot.caption     = element_text(size = BASE_SIZE - 2.5,
                                    color = "black"),
    panel.grid.major = element_line(color = PAL$rule, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(color = "black", face = "bold",
                                    size  = BASE_SIZE - 1.5),
    axis.title       = element_text(color = "black", face = "bold",
                                    size  = BASE_SIZE - 1),
    legend.text      = element_text(size = BASE_SIZE - 1.5, color = "black"),
    legend.title     = element_text(size = BASE_SIZE - 1, face = "bold", color = "black"),
    plot.background  = element_rect(fill = PAL$bg, color = NA),
    panel.background = element_rect(fill = PAL$bg, color = NA)
  )

# Subcluster palette, shared by Fig 1A/1B and S1C
CLUSTER_COLORS <- c(
  "0" = PAL$textMute, "1" = "#7a8c6a", "2" = PAL$drug,
  "3" = "#9b59b6", "4" = PAL$accent, "5" = PAL$disease, "6" = "#8B5A3C"
)

# Candidate-status palette for the top-30 THP1 panel
TIER_COLORS <- c(
  "Highlighted post hoc candidate"     = PAL$accent,
  "FDA-approved"     = "#0072B2",
  "Not FDA-approved" = "#808080"
)

# Route each figure by its FINAL stem; never infer destinations from execution order.
figure_group <- function(stem) {
  if (grepl("^fig[1-5]_", stem)) return("main")
  if (grepl("^suppfig_S[1-8]_", stem)) return("supplementary")
  stop("Unrecognized final figure filename: ", stem)
}
figure_dirs <- function(stem) {
  root <- OUTPUT_ROOTS[[figure_group(stem)]]
  c(with_titles = file.path(root, "with_titles"), clean = file.path(root, "clean"))
}
for (root in OUTPUT_ROOTS) {
  for (variant in c("clean", "with_titles")) {
    dir.create(file.path(root, variant), showWarnings = FALSE, recursive = TRUE)
  }
}
register_exports <- function(stem, paths) {
  EXPORT_RECORDS[[length(EXPORT_RECORDS) + 1L]] <<- data.frame(
    group = figure_group(stem), stem = stem, variant = basename(dirname(paths)),
    format = tools::file_ext(paths), path = paths, stringsAsFactors = FALSE)
}

# Text is independent of the data palette: gray marks/gridlines stay gray.
# Apply to every leaf plot before saving, including text/repel annotations.
black_plot_text <- function(p) {
  if (inherits(p, "patchwork")) {
    for (j in seq_len(length(p))) p[[j]] <- black_plot_text(p[[j]])
    return(p & ggplot2::theme(text = ggplot2::element_text(colour = "black")))
  }
  for (j in seq_along(p$layers)) {
    layer <- p$layers[[j]]
    if (inherits(layer$geom, "GeomText") || inherits(layer$geom, "GeomLabel") ||
        inherits(layer$geom, "GeomTextRepel") || inherits(layer$geom, "GeomLabelRepel")) {
      # White lettering on dark schematic boxes needs its existing contrast.
      col <- layer$aes_params$colour
      white <- is.character(col) && length(col) == 1L &&
        tolower(col) %in% c("white", "#ffffff")
      if (!white) {
        layer$aes_params$colour <- "black"
        layer$aes_params$alpha <- 1
      }
      p$layers[[j]] <- layer
    }
  }
  # Respect deliberately blank title elements in clean variants.
  tt <- p$theme
  for (nm in names(tt)) {
    if (inherits(tt[[nm]], "element_text")) tt[[nm]]$colour <- "black"
  }
  p$theme <- tt
  p + ggplot2::theme(text = ggplot2::element_text(colour = "black"))
}

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
  dirs <- figure_dirs(name)
  p <- black_plot_text(p)
  override <- SIZE_OVERRIDES[[name]]
  if (!is.null(override)) {
    if (is.list(override)) override <- unlist(override, use.names = FALSE)
    if (length(override) != 2L || any(!is.finite(as.numeric(override))) ||
        any(as.numeric(override) <= 0)) {
      stop("Size override for '", name, "' must be c(width, height) > 0.")
    }
    w <- as.numeric(override[[1]])
    h <- as.numeric(override[[2]])
  }
  paths <- c(
    file.path(dirs[["with_titles"]], paste0(name, ".pdf")),
    file.path(dirs[["with_titles"]], paste0(name, ".png")),
    file.path(dirs[["clean"]], paste0(name, ".pdf")),
    file.path(dirs[["clean"]], paste0(name, ".png"))
  )
  pc <- strip_titles(p)
  plots <- list(p, p, pc, pc)

  export_indices <- if (SAVE_TITLED) seq_along(paths) else c(3L, 4L)
  for (i in export_indices) {
    message("Exporting: ", paths[i])
    if (grepl("\\.png$", paths[i])) {
      ggplot2::ggsave(
        filename = paths[i],
        plot = plots[[i]],
        width = w,
        height = h,
        dpi = dpi
      )
    } else {
      ggplot2::ggsave(
        filename = paths[i],
        plot = plots[[i]],
        width = w,
        height = h,
        device = open_pdf_device
      )
    }
  }

  file_ok <- file.exists(paths) &
    !is.na(file.info(paths)$size) &
    file.info(paths)$size > 0

  if (!SAVE_TITLED) file_ok[c(1L, 2L)] <- TRUE

  for (i in export_indices) {
    cat(if (file_ok[i]) "  OK:   " else "  FAIL: ", paths[i], "\n")
  }

  if (!all(file_ok)) {
    stop("save_panel('", name, "') failed to create one or more non-empty files.")
  }

  register_exports(name, paths[export_indices])
  message("Saved and verified: ", name)
}

save_combined <- function(p, name, w, h, dpi = COMBINED_DPI) {
  save_panel(p, name, w, h, dpi)
}

# =============================================================================
# SECTION 1 · LOAD AND VALIDATE FINALIZED ANALYSIS OBJECTS
# =============================================================================

read_figure_rds <- function(path) {
  message("Reading: ", path)
  started <- Sys.time()
  value <- readRDS(path)
  message("Read complete (", round(as.numeric(difftime(Sys.time(), started,
                                                       units = "secs")), 1), " seconds): ", basename(path))
  value
}

message("Loading finalized analysis objects...")

analysis_paths <- c(
  inph_signature = file.path(PROCESSED_DIR, "inph_ifn_signature.rds"),
  de_results = file.path(PROCESSED_DIR, "de_ifn_subtype_results.rds"),
  subclustered_monocytes = file.path(PROCESSED_DIR, "seu_monocyte_subclustered.rds"),
  contextualization = file.path(PROCESSED_DIR, "phase2_contextualization.rds"),
  lincs_result = file.path(PROCESSED_DIR, "lincs2_gess_result.rds"),
  compound_summary = file.path(PROCESSED_DIR, "lincs2_compound_summary.rds"),
  final_candidates = file.path(PROCESSED_DIR, "lincs2_final_candidates.rds"),
  loo_stability = file.path(PROCESSED_DIR, "loo_signature_stability.rds"),
  loo_persistence = file.path(PROCESSED_DIR, "loo_up_gene_persistence.csv")
)

analysis_info <- file.info(unname(analysis_paths))

bad_analysis_inputs <- names(analysis_paths)[
  !file.exists(unname(analysis_paths)) |
    is.na(analysis_info$size) |
    analysis_info$size <= 0
]

if (length(bad_analysis_inputs) > 0L) {
  stop(
    "Required finalized analysis input(s) missing or empty: ",
    paste(bad_analysis_inputs, collapse = ", ")
  )
}

inph_signature <- read_figure_rds(analysis_paths[["inph_signature"]])
de_res <- read_figure_rds(analysis_paths[["de_results"]])
seu <- read_figure_rds(analysis_paths[["subclustered_monocytes"]])

contextualization_results <- read_figure_rds(
  analysis_paths[["contextualization"]]
)

CTX_A <- contextualization_results$engine_A_primary
CTX_B <- contextualization_results$engine_B_sensitivity

result_table <- read_figure_rds(
  analysis_paths[["lincs_result"]]
) |>
  signatureSearch::result()

# Script 03 now saves compound names after the finalized two-source metadata
# resolution. No downstream name backfill is needed here.
compound_named <- read_figure_rds(
  analysis_paths[["compound_summary"]]
)

final_cands <- read_figure_rds(
  analysis_paths[["final_candidates"]]
)

loo_df <- read_figure_rds(
  analysis_paths[["loo_stability"]]
)

persist_df <- readr::read_csv(
  analysis_paths[["loo_persistence"]],
  show_col_types = FALSE
)

if (
  !is.list(inph_signature) ||
  length(inph_signature$up) != 150L ||
  length(inph_signature$down) != 150L ||
  anyDuplicated(inph_signature$up) ||
  anyDuplicated(inph_signature$down)
) {
  stop("Phase 1 signature failed the expected 150-UP/150-DOWN validation.")
}

if (
  !inherits(seu, "Seurat") ||
  ncol(seu) != 22197L ||
  length(unique(as.character(seu$patient))) != 9L ||
  sum(as.character(seu$sub_res05) == "0") != 4904L ||
  sum(as.character(seu$sub_res05) == "5") != 940L
) {
  stop("Subclustered monocyte object does not match the finalized analysis.")
}

if (
  is.null(CTX_A) ||
  is.null(CTX_B)
) {
  stop(
    "phase2_contextualization.rds must contain both primary and sensitivity engines."
  )
}

required_lincs_columns <- c(
  "pert", "type", "cell", "NCS", "WTCS", "WTCS_FDR"
)

if (
  nrow(result_table) != 136460L ||
  length(setdiff(required_lincs_columns, colnames(result_table))) > 0L
) {
  stop("LINCS2 result object does not match the finalized Script 03 output.")
}

if (
  nrow(compound_named) != 30411L ||
  anyDuplicated(compound_named$pert) ||
  sum(!is.na(compound_named$thp1_NCS) & compound_named$thp1_NCS != 0) != 597L
) {
  stop("Compound summary does not match the finalized 30,411/597 accounting.")
}

expected_final_candidate_ids <- c(
  "BRD-K82941592",
  "BRD-K03194791",
  "BRD-K53581288",
  "BRD-A02367930",
  "BRD-K75641298",
  "BRD-K60038276",
  "BRD-K77987382",
  "BRD-K29950728"
)

if (
  nrow(final_cands) != 8L ||
  !setequal(final_cands$pert, expected_final_candidate_ids)
) {
  stop("Final candidate object does not match the finalized eight-compound set.")
}

required_loo_columns <- c(
  "dropped", "n_sig", "jaccard_up", "jaccard_dn", "jaccard_mean"
)

if (
  nrow(loo_df) != 9L ||
  length(setdiff(required_loo_columns, colnames(loo_df))) > 0L ||
  abs(mean(loo_df$jaccard_up) - 0.714) > 0.002 ||
  abs(mean(loo_df$jaccard_dn) - 0.733) > 0.002
) {
  stop("LOO stability object does not match the finalized Script 04 output.")
}

if (
  !all(c("gene", "n_iters", "in_sig") %in% colnames(persist_df)) ||
  sum(persist_df$in_sig) != 137L ||
  sum(persist_df$in_sig & persist_df$n_iters == 9L) != 107L
) {
  stop("LOO UP-gene persistence table does not match the finalized 137/107 result.")
}

cell_meta <- data.frame(
  cell = colnames(seu),
  patient = seu$patient,
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
# SECTION 2 - READ THE EXISTING S1 PBMC CACHE (NO RECONSTRUCTION)
# =============================================================================
PBMC_SCE_PATH <- file.path(PROCESSED_DIR, "pbmc_reconstructed_historical_qc.rds")
if (!file.exists(PBMC_SCE_PATH)) {
  stop("Missing saved S1 PBMC cache: ", PBMC_SCE_PATH,
       ". Restore the validated cache before running this export script.")
}
pbmc_sce <- tryCatch(read_figure_rds(PBMC_SCE_PATH), error = function(e) {
  stop("The saved S1 cache could not be read: ", conditionMessage(e),
       ". No cache was renamed or rebuilt. Restore a readable cache and run again.",
       call. = FALSE)
})
canonical_ids <- colnames(seu)

stopifnot(
  inherits(pbmc_sce, "SingleCellExperiment"),
  ncol(pbmc_sce) == 122724L,
  !anyDuplicated(colnames(pbmc_sce)),
  "UMAP" %in%
    SingleCellExperiment::reducedDimNames(pbmc_sce),
  length(canonical_ids) == 22197L,
  !anyDuplicated(canonical_ids),
  all(canonical_ids %in% colnames(pbmc_sce))
)

pbmc_sce$cell_type <- ifelse(
  colnames(pbmc_sce) %in% canonical_ids,
  "monocyte",
  "other"
)

stopifnot(
  sum(pbmc_sce$cell_type == "monocyte") == 22197L,
  sum(pbmc_sce$cell_type == "other") == 100527L
)

S1_CACHE_READY <- TRUE
S1_CACHE_NOTE <- "Validated existing PBMC cache loaded; no clustering run."
S1_REUSE_RECORDS <- list()
message("S1 cache validated: 122,724 PBMCs; 22,197 canonical monocytes.")

# -----------------------------------------------------------------------------
# Locked revision inputs used by reviewer-facing replacements
# -----------------------------------------------------------------------------
REVISION_RESULTS_DIR <- here::here("results", "revision")
REVISION_FIGURES_DIR <- here::here("figures", "revision")
REVISION_PROCESSED_DIR <- here::here("data", "processed")
EXPECTED_EVIDENCE_SHA <-
  "21e6225b9e86bda29896f11a55e61d1618e9ce5a57dc6d271a26e99e358abf95"

find_revision_file <- function(filename, roots = REVISION_RESULTS_DIR) {
  roots <- roots[file.exists(roots)]
  hits <- unique(unlist(lapply(roots, function(root) {
    direct <- file.path(root, filename)
    recursive <- list.files(
      root,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    c(direct, recursive)
  }), use.names = FALSE))
  hits <- hits[file.exists(hits) & !dir.exists(hits)]
  if (!length(hits)) {
    stop("Required locked revision file not found: ", filename)
  }
  info <- file.info(hits)
  hits[[which.max(info$mtime)]]
}

find_revision_asset <- function(filename, roots = REVISION_FIGURES_DIR) {
  roots <- roots[file.exists(roots)]
  hits <- unique(unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
  }), use.names = FALSE))
  hits <- hits[file.exists(hits) & !dir.exists(hits)]
  if (!length(hits)) return(NA_character_)
  info <- file.info(hits)
  hits[[which.max(info$mtime)]]
}

evidence_lock_path <- find_revision_file("revision_evidence_lock.txt")
if (!identical(
  digest::digest(file = evidence_lock_path, algo = "sha256"),
  EXPECTED_EVIDENCE_SHA
)) {
  stop("The available evidence lock does not match the approved SHA-256.")
}

marker_state_path <- file.path(
  REVISION_PROCESSED_DIR,
  "revision_marker_ifn_states.rds"
)
if (!file.exists(marker_state_path)) {
  marker_state_path <- find_revision_file(
    "revision_marker_ifn_states.rds",
    roots = c(REVISION_PROCESSED_DIR, REVISION_RESULTS_DIR)
  )
}
marker_state <- read_figure_rds(marker_state_path)
marker_donor_counts_all <- as.data.frame(marker_state$donor_summary)
marker_donor_counts <- marker_donor_counts_all |>
  dplyr::filter(cohort == "canonical_nine_patient")
marker_cluster_concordance <- as.data.frame(
  marker_state$canonical_cluster_concordance
)

required_marker_columns <- c(
  "patient", "total_cells", "marker_ifn_high_n", "marker_ifn_high_pct"
)
if (length(setdiff(required_marker_columns, colnames(marker_donor_counts))) > 0L) {
  stop("Marker-defined donor table is missing required columns.")
}
if (
  nrow(marker_donor_counts) != 9L ||
  sum(marker_donor_counts$marker_ifn_high_n, na.rm = TRUE) != 3195L
) {
  stop("Marker-defined donor table does not match the locked 3,195-cell total.")
}
if (
  nrow(marker_cluster_concordance) != 1L ||
  marker_cluster_concordance$marker_high_and_cluster5_n[[1]] != 898L
) {
  stop("Marker/cluster concordance table does not match the locked result.")
}

# Validate the final average-tie object used for main Figure 4 and reviewer figures.
average_tie_path <- file.path(REVISION_RESULTS_DIR, "average_tie_finalization",
  "run_20260907_104848_4a48653337e1", "revision_candidate_robustness_average_ties.rds")
if (!file.exists(average_tie_path)) stop("Missing locked average-tie object: ", average_tie_path)
average_tie_result <- read_figure_rds(average_tie_path)
if (!identical(as.character(average_tie_result$ranking_method), "average")) {
  stop("The finalized robustness object is not the approved average-tie result.")
}

message("Locked revision evidence validated.")


W_UMAP    <- 4.2;  H_UMAP    <- 4.0
W_ISG     <- 3.5;  H_ISG     <- 4.0
# Per-panel dimensions in inches; named SIZE_OVERRIDES use the final filenames below.
W_MoCA    <- 3.4;  H_MoCA    <- 3.8
W_MARKER  <- 4.0;  H_MARKER  <- 3.8
W_VOLCANO <- 9.5; H_VOLCANO <- 5.5
W_NCS     <- 5.5;  H_NCS     <- 4.5
W_THP1    <- 7.5;  H_THP1    <- 9.5
W_CELLLINE<- 6.0;  H_CELLLINE<- 6.0
W_LOO     <- 8.0;  H_LOO     <- 4.5
W_PERSIST <- 7.0;  H_PERSIST <- 6.0

# -----------------------------------------------------------------------------
# Helper: cell-line NCS profile for one compound (Fig 5A/5B, S6)
# -----------------------------------------------------------------------------
cellline_panel <- function(brd_id, drug_label, subtitle_extra = NULL,
                           fdr_thresh = 0.05, highlight = PAL$accent) {
  df <- result_table |>
    dplyr::filter(
      pert == brd_id,
      type == "trt_cp",
      !is.na(NCS)
    ) |>
    dplyr::mutate(
      sig = !is.na(WTCS_FDR) & WTCS_FDR < fdr_thresh,
      direction = dplyr::case_when(
        NCS < 0 ~ "Negative connectivity",
        NCS > 0 ~ "Positive connectivity",
        TRUE ~ "NCS = 0"
      ),
      cell_label = factor(
        cell,
        levels = cell[order(NCS)]
      )
    )

  if (
    nrow(df) == 0L ||
    anyDuplicated(df$cell)
  ) {
    stop(
      "Cell-line profile is missing or duplicated for compound ",
      brd_id,
      "."
    )
  }

  thp1_idx <- which(toupper(levels(df$cell_label)) == "THP1")
  n_lines  <- nrow(df)
  n_rev    <- sum(df$NCS < 0 & df$sig)
  thp1_ncs <- df$NCS[toupper(df$cell) == "THP1"]

  sub <- paste0(n_lines, " LINCS2 cell lines | significant reversal in ",
                n_rev, " | THP1 NCS = ",
                if (length(thp1_ncs)) round(thp1_ncs[1], 2) else "not profiled")
  if (!is.null(subtitle_extra)) sub <- paste0(sub, "\n", subtitle_extra)

  # Display only the saved FDR value; retain NA when unavailable.
  df$fdr_label <- ifelse(
    is.na(df$WTCS_FDR),
    "NA",
    formatC(df$WTCS_FDR, digits = 2, format = "e")
  )
  # Zero-score rows have no bar. Missing FDRs are labeled, never imputed.
  label_df <- df[df$NCS != 0, , drop = FALSE]

  p <- ggplot(df, aes(x = NCS, y = cell_label, fill = direction))
  if (length(thp1_idx) == 1L) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = thp1_idx - 0.5, ymax = thp1_idx + 0.5,
                      fill = highlight, alpha = 0.10)
  }
  p <- p +
    geom_col(width = 0.62, alpha = 0.9) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_text(data = label_df,
              aes(label = fdr_label, hjust = ifelse(NCS < 0, 1.05, -0.05)),
              size = 3.0, color = "black") +
    scale_fill_manual(
      values = c("Negative connectivity" = PAL$disease,
                 "NCS = 0" = PAL$textMute,
                 "Positive connectivity" = PAL$drug),
      limits = c("Negative connectivity", "NCS = 0", "Positive connectivity"),
      breaks = c("Negative connectivity", "NCS = 0", "Positive connectivity"),
      drop = FALSE, name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.34, 0.34))) +
    coord_cartesian(clip = "off") +
    labs(x = "Normalized Connectivity Score", y = NULL,
         title = paste0(drug_label, " \u2014 NCS across LINCS2 cell lines"),
         subtitle = sub,
         caption = paste0(brd_id, " = ", tolower(drug_label),
                          "; labels show saved FDR for every nonzero bar, including FDR >= 0.05; NA = unavailable.")) +
    THEME_BASE +
    theme(legend.position = "top", panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 30, b = 5.5, l = 15))

  p
}

# MAIN FIGURES: 1 through 5

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
    theme(
      panel.grid = element_blank(),

      axis.text = element_text(
        color = "black",
        face = "bold",
        size = BASE_SIZE - 1.5
      ),

      axis.title = element_text(
        color = "black",
        face = "bold",
        size = BASE_SIZE - 1
      ),

      axis.ticks = element_line(
        color = "black",
        linewidth = 0.65
      ),

      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),

      axis.ticks.length = grid::unit(2, "mm")
    )
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
# FIG 1C — Descriptive IFN-response-high frequency by MoCA group
#
# Retained for continuity with the original manuscript and reviewer discussion.
# This is descriptive only: exact patient-linked MoCA scores are unavailable,
# the group sizes are small, and the locked comparison was not significant
# (P = 0.556). No cognition association is inferred.
# -----------------------------------------------------------------------------
{
  P1C_HIGH_PATIENTS <- c(1, 2, 3, 4, 5)
  P1C_LOW_PATIENTS  <- c(7, 8, 9, 10)
  P1C_MOCA_COLORS   <- c("High" = PAL$mocaHigh, "Low" = PAL$mocaLow)

  moca_prop_df <- cell_meta |>
    dplyr::filter(moca_group %in% c("High", "Low")) |>
    dplyr::group_by(patient) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_ifn = sum(subcluster == "5"),
      ifn_pct = 100 * n_ifn / n_total,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      moca_grp = dplyr::case_when(
        patient %in% P1C_HIGH_PATIENTS ~ "High",
        patient %in% P1C_LOW_PATIENTS ~ "Low",
        TRUE ~ "Unknown"
      )
    ) |>
    dplyr::filter(moca_grp != "Unknown") |>
    dplyr::mutate(
      moca_grp = factor(moca_grp, levels = c("High", "Low")),
      patient_label = paste0("P", patient)
    )

  if (nrow(moca_prop_df) != 9L || anyNA(moca_prop_df$ifn_pct)) {
    stop("Figure 1C MoCA grouping failed locked nine-donor validation.")
  }

  moca_test <- stats::wilcox.test(
    ifn_pct ~ moca_grp, data = moca_prop_df, exact = TRUE
  )
  if (abs(as.numeric(moca_test$p.value) - 0.556) > 0.002) {
    stop("Figure 1C MoCA comparison differs from locked P = 0.556.")
  }

  p_moca_box <- ggplot(
    moca_prop_df, aes(x = moca_grp, y = ifn_pct, fill = moca_grp)
  ) +
    geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(aes(color = moca_grp), width = 0.08, size = 2.8, alpha = 0.9) +
    geom_text(aes(label = patient_label), vjust = -0.8, size = 2.3,
              color = "black") +
    annotate(
      "text", x = 1.5, y = max(moca_prop_df$ifn_pct) * 1.16,
      label = "P = 0.556",
      size = 2.7, color = "black", fontface = "italic"
    ) +
    scale_fill_manual(values = P1C_MOCA_COLORS, guide = "none") +
    scale_color_manual(values = P1C_MOCA_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.22))) +
    labs(
      title = "Descriptive IFN-response-high frequency by MoCA group",
      subtitle = "Group-level MoCA labels; exact patient-linked scores unavailable",
      x = "MoCA group", y = "Cluster-5 cells (%)"
    ) +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank())
}
save_panel(p_moca_box, "fig1_panelC_moca_box", W_MoCA, H_MoCA)

# -----------------------------------------------------------------------------
# FIG 1D — Marker-defined IFN-response-high frequency by donor
# -----------------------------------------------------------------------------
{
  marker_freq_df <- marker_donor_counts |>
    dplyr::mutate(
      patient = as.integer(patient),
      patient_label = factor(
        paste0("P", patient),
        levels = paste0("P", sort(unique(patient)))
      ),
      low_cell_yield = total_cells < 30L
    ) |>
    dplyr::arrange(patient)

  if (
    nrow(marker_freq_df) != 9L ||
    anyNA(marker_freq_df$marker_ifn_high_pct) ||
    sum(marker_freq_df$marker_ifn_high_n) != 3195L
  ) {
    stop("Figure 1D marker-defined donor frequencies failed locked validation.")
  }

  stopifnot(all(marker_freq_df$total_cells >= 30L),
            sum(marker_freq_df$total_cells) == 22197L)
  # All canonical donors meet this yield threshold; it needs no legend.
  p_marker_freq <- ggplot(
    marker_freq_df,
    aes(x = patient_label, y = marker_ifn_high_pct)
  ) +
    geom_col(width = 0.72, color = NA, fill = PAL$ifn) +
    geom_text(
      aes(label = paste0(round(marker_ifn_high_pct, 1), "%")),
      vjust = -0.35,
      size = 2.7,
      color = PAL$text
    ) +
    scale_y_continuous(
      labels = scales::label_number(scale = 1),
      limits = c(0, max(marker_freq_df$marker_ifn_high_pct) * 1.18),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Marker-defined IFN-response-high monocytes by donor",
      subtitle = paste0(
        "Descriptive within-iNPH distribution; total marker-high cells = ",
        sum(marker_freq_df$marker_ifn_high_n)
      ),
      x = "Donor",
      y = "Marker-defined IFN-response-high cells (%)"
    ) +
    THEME_BASE +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

}
save_panel(p_marker_freq, "fig1_panelD_marker_frequency", W_MARKER, H_MARKER)

# -----------------------------------------------------------------------------
# FIG 1E — Volcano
# -----------------------------------------------------------------------------
{
  P1E_FDR_THRESH <- 0.05
  P1E_LFC_THRESH <- 1.0

  # Reduce the number of displayed labels to avoid crowding.
  P1E_N_LABEL_UP <- 14
  P1E_N_LABEL_DN <- 6
  P1E_CAP_Y      <- 100

  # Typography tuned for placement as Figure 1E.
  P1E_GENE_LABEL_SIZE  <- 3.0
  P1E_HEADER_SIZE      <- 4.0
  P1E_CAP_NOTE_SIZE    <- 1.8
  P1E_AXIS_TEXT_SIZE   <- 15
  P1E_AXIS_TITLE_SIZE  <- 15.0
  P1E_LEGEND_TEXT_SIZE <- 10
  P1E_TITLE_SIZE       <- 11.0
  P1E_SUBTITLE_SIZE    <- 9


  volc_df <- as.data.frame(de_res) |>
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
    dplyr::mutate(
      neg_log10_fdr = pmin(-log10(padj), P1E_CAP_Y),
      is_capped     = -log10(padj) > P1E_CAP_Y,
      direction = dplyr::case_when(
        padj < P1E_FDR_THRESH & log2FoldChange >  P1E_LFC_THRESH ~ "UP",
        padj < P1E_FDR_THRESH & log2FoldChange < -P1E_LFC_THRESH ~ "DOWN",
        TRUE ~ "NS"),
      direction = factor(direction, levels = c("UP","DOWN","NS")))

  if (!"gene" %in% colnames(volc_df))
    volc_df <- tibble::rownames_to_column(volc_df, "gene")

  top_up  <- volc_df |> dplyr::filter(direction == "UP") |>
    dplyr::arrange(dplyr::desc(log2FoldChange)) |> head(P1E_N_LABEL_UP)
  top_dn  <- volc_df |> dplyr::filter(direction == "DOWN") |>
    dplyr::arrange(log2FoldChange) |> head(P1E_N_LABEL_DN)
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
               aes(y = P1E_CAP_Y), size = 2.0, shape = 24,
               fill = PAL$disease, color = PAL$disease, alpha = 0.9) +
    ggrepel::geom_text_repel(data = label_df, aes(label = gene),
                             size = P1E_GENE_LABEL_SIZE,
                             max.overlaps = 35, segment.color = PAL$textMute,
                             segment.size = 0.25, box.padding = 0.30) +
    geom_hline(yintercept = -log10(P1E_FDR_THRESH), linetype = "dashed",
               color = PAL$textMute, linewidth = 0.4) +
    geom_vline(xintercept = c(-P1E_LFC_THRESH, P1E_LFC_THRESH),
               linetype = "dashed", color = PAL$textMute, linewidth = 0.35) +
    geom_vline(xintercept = 0, color = PAL$rule, linewidth = 0.35) +
    scale_color_manual(
      values = c("UP" = PAL$disease, "DOWN" = PAL$drug, "NS" = PAL$textMute),
      labels = c("UP" = paste0("UP (n=", n_up, ")"),
                 "DOWN" = paste0("DOWN (n=", n_dn, ")"),
                 "NS" = "Not significant"),
      name = NULL,
      guide = guide_legend(override.aes = list(shape = 16, size = 3))) +
    coord_cartesian(
      xlim = c(-max_abs_x, max_abs_x),
      ylim = c(0, P1E_CAP_Y + 5)
    ) +

    annotate(
      "text",
      x = max_abs_x * 0.65,
      y = P1E_CAP_Y + 3,
      label = "IFN-response-high (UP)",
      hjust = 0.5,
      size = P1E_HEADER_SIZE,
      color = "black",
      fontface = "plain"
    ) +

    annotate(
      "text",
      x = -max_abs_x * 0.65,
      y = P1E_CAP_Y + 3,
      label = "Classical resting (DOWN)",
      hjust = 0.5,
      size = P1E_HEADER_SIZE,
      color = "black",
      fontface = "plain"
    ) +

    labs(
      x = expression(
        bold("Shrunken " * log[2] * " fold change (apeglm)")
      ),
      y = expression(
        bold(-log[10] * " (FDR)")
      ),
      title =
        "Differential expression: submitted cluster 5 vs classical resting",
      subtitle = paste0(
        "DESeq2 | ~ patient + state | cluster 5 vs cluster 0 | ",
        "1,231 genes at FDR < 0.05 | 9 patients"
      )
    ) +

    THEME_BASE +

    theme(
      legend.position = "right",

      axis.line.x = element_line(
        color = "black",
        linewidth = 0.8
      ),

      axis.line.y = element_line(
        color = "black",
        linewidth = 0.8
      ),

      axis.ticks.x = element_line(
        color = "black",
        linewidth = 0.65
      ),

      axis.ticks.y = element_line(
        color = "black",
        linewidth = 0.65
      ),

      axis.text.x = element_text(
        family = "sans",
        color = "black",
        face = "bold",
        size = P1E_AXIS_TEXT_SIZE
      ),

      axis.text.y = element_text(
        family = "sans",
        color = "black",
        face = "bold",
        size = P1E_AXIS_TEXT_SIZE
      ),

      axis.title.x = element_text(
        family = "sans",
        color = "black",
        face = "bold",
        size = P1E_AXIS_TITLE_SIZE,
        margin = margin(t = 8)
      ),

      axis.title.y = element_text(
        family = "sans",
        color = "black",
        face = "bold",
        size = P1E_AXIS_TITLE_SIZE,
        margin = margin(r = 8)
      ),

      legend.text = element_text(
        family = "sans",
        color = "black",
        face = "plain",
        size = P1E_LEGEND_TEXT_SIZE
      ),

      legend.key.size = grid::unit(4, "mm")
    )
}
save_panel(p_volcano, "fig1_panelE_volcano", W_VOLCANO, H_VOLCANO)

fig1 <- (p_umap + p_isg_violin + p_moca_box +
           patchwork::plot_layout(widths = c(1.25, 1.05, 1))) /
  (p_marker_freq + p_volcano + patchwork::plot_layout(widths = c(1, 2.55))) +
  patchwork::plot_layout(heights = c(1, 1.15)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig1, "fig1_combined", 12, 11)

# -----------------------------------------------------------------------------
# FIG 2 — Cross-disease contextualization (main Figure 2)
# -----------------------------------------------------------------------------
{
  # ============================================================
  # TUNABLE KNOBS
  # ============================================================
  F2_TEXT_SCALE  <- 1.15
  F2_ROW_HEIGHTS <- c(1, 0.65)
  F2_BOX_WIDTH   <- 0.6
  W_F2_COMBINED  <- 12
  H_F2_COMBINED  <- 9

  # Keep the F2C legend outside the plotting region so the combined
  # patchwork layout remains stable across output sizes.

  f2_text_theme <- theme(
    plot.title    = element_text(size = (BASE_SIZE + 1)   * F2_TEXT_SCALE),
    plot.subtitle = element_text(size = (BASE_SIZE - 1.5) * F2_TEXT_SCALE),
    plot.caption  = element_text(size = (BASE_SIZE - 2.5) * F2_TEXT_SCALE),
    axis.text     = element_text(size = (BASE_SIZE - 1.5) * F2_TEXT_SCALE),
    axis.title    = element_text(size = (BASE_SIZE - 1)   * F2_TEXT_SCALE),
    legend.text   = element_text(size = (BASE_SIZE - 2.5) * F2_TEXT_SCALE),
    legend.title  = element_text(size = (BASE_SIZE - 1)   * F2_TEXT_SCALE)
  )

  F2_COLORS <- c("Healthy control" = PAL$neutral,
                 "SLE"             = PAL$drug,
                 "iNPH"            = PAL$disease)

  # ── A & B: population-level activity and donor-level burden ──
  pop_df <- dplyr::bind_rows(
    CTX_A$pop_tables$hc_hi   |> dplyr::mutate(group = "Healthy control"),
    CTX_A$pop_tables$sle_hi  |> dplyr::mutate(group = "SLE"),
    CTX_A$pop_tables$inph_c5 |> dplyr::mutate(group = "iNPH")) |>
    dplyr::mutate(group = factor(group, levels = names(F2_COLORS)))

  p_context_a <- ggplot(pop_df, aes(x = group, y = mean_sig, fill = group)) +
    geom_boxplot(width = F2_BOX_WIDTH, outlier.shape = NA, alpha = 0.55, color = "grey25") +
    geom_jitter(width = 0.10, size = 1.4, alpha = 0.5, color = PAL$textSoft) +
    scale_fill_manual(values = F2_COLORS, guide = "none") +
    labs(x = NULL, y = "AUCell signature score (donor mean)",
         title = "Population-level IFN signature score",
         subtitle = paste0(
           "Donor-level, unfiltered cluster 5 | medians: iNPH ",
           round(CTX_A$med_pop["inph_c5"], 4), ", SLE-high ",
           round(CTX_A$med_pop["sle_hi"], 4), ", HC-high ",
           round(CTX_A$med_pop["hc_hi"], 4), "\n",
           "iNPH vs SLE-high p = ", signif(CTX_A$p_c5_sle, 3),
           " | iNPH vs HC-high p = ", signif(CTX_A$p_c5_hc, 3)),
         caption = "Primary AUCell analysis with the 15 selector genes excluded from signature scoring.") +
    THEME_BASE + f2_text_theme +
    theme(panel.grid.major.x = element_blank())

  burden_df <- dplyr::bind_rows(
    CTX_A$burden_tables$hc   |> dplyr::mutate(group = "Healthy control"),
    CTX_A$burden_tables$sle  |> dplyr::mutate(group = "SLE"),
    CTX_A$burden_tables$inph |> dplyr::mutate(group = "iNPH")) |>
    dplyr::mutate(group = factor(group, levels = names(F2_COLORS)))

  p_context_b <- ggplot(burden_df, aes(x = group, y = burden_pct, fill = group)) +
    geom_boxplot(width = F2_BOX_WIDTH, outlier.shape = NA, alpha = 0.55, color = "grey25") +
    geom_jitter(width = 0.10, size = 1.4, alpha = 0.5, color = PAL$textSoft) +
    scale_fill_manual(values = F2_COLORS, guide = "none") +
    labs(x = NULL, y = "IFN-high monocytes per donor (%)",
         title = "IFN-high monocyte burden",
         subtitle = paste0(
           "Donor-level, shared healthy-anchored threshold | medians: HC ",
           round(CTX_A$burden_med["hc"], 1), "%, SLE ",
           round(CTX_A$burden_med["sle"], 1), "%, iNPH ",
           round(CTX_A$burden_med["inph"], 1), "%\n",
           "iNPH vs SLE p = ", signif(CTX_A$p_burden_sle, 3),
           " | iNPH vs HC p = ", signif(CTX_A$p_burden_hc, 3),
           " (HC comparison is scoring-method-dependent — see Methods)"),
         caption = paste0(
           "Exploratory: the iNPH comparison crosses the Duy/Perez dataset ",
           "boundary and is batch-sensitive.")) +
    THEME_BASE + f2_text_theme +
    theme(panel.grid.major.x = element_blank())

  # ── C: cluster-5 / ISG-high overlap ──
  ov <- CTX_A$overlap
  overlap_df <- data.frame(
    category = factor(c("Above shared IFN-high threshold", "Below shared IFN-high threshold"),
                      levels = c("Above shared IFN-high threshold", "Below shared IFN-high threshold")),
    n = c(ov[["n_both"]], ov[["n_c5"]] - ov[["n_both"]]))

  p_context_c <- ggplot(overlap_df, aes(x = "", y = n, fill = category)) +
    geom_col(width = 0.92) +
    scale_x_discrete(expand = expansion(add = 0.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_fill_manual(values = c("Above shared IFN-high threshold" = PAL$disease,
                                 "Below shared IFN-high threshold" = PAL$neutral), name = NULL) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Cluster 5 cells (n = 940)",
         title = "Concordance of cluster 5 with the shared IFN-high threshold",
         subtitle = paste0(
           sprintf("%.1f%%", ov[["pct_c5_is_ifnhi"]]),
           " of cluster-5 cells clear the shared IFN-high threshold ",
           "(AUCell; ", sprintf("%.1f%%", CTX_B$overlap[["pct_c5_is_ifnhi"]]),
           " under mean-difference)\n",
           "Population-level score: iNPH ",
           round(CTX_A$med_pop["inph_c5"], 4), " vs SLE-high ",
           round(CTX_A$med_pop["sle_hi"], 4),
           " (p = ", signif(CTX_A$p_c5_sle, 3), ")"),
         caption = paste0(
           "Descriptive concordance; selector genes were excluded from ",
           "the iNPH signature-scoring gene set.")) +
    THEME_BASE + f2_text_theme +
    theme(
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = "right"
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))

  save_panel(p_context_a, "fig2_panelA_cross_disease_population", 4, 4)
  save_panel(p_context_b, "fig2_panelB_cross_disease_burden",  4, 4)
  save_panel(p_context_c, "fig2_panelC_cross_disease_concordance", 9, 2)

  fig2 <- (p_context_a + p_context_b + patchwork::plot_layout(widths = c(1, 1))) / p_context_c +
    patchwork::plot_layout(heights = F2_ROW_HEIGHTS) +
    patchwork::plot_annotation(tag_levels = "A")
  save_combined(fig2, "fig2_combined", W_F2_COMBINED, H_F2_COMBINED)
  message("Main Figure 2 saved.")
}

# -----------------------------------------------------------------------------
# FIG 3A — THP1 NCS distribution (zero-aware)
# -----------------------------------------------------------------------------
{
  P2A_BINS      <- 50
  P2A_THRESHOLD <- -0.5

  thp1_all <- result_table |>
    dplyr::filter(toupper(cell) == "THP1", type == "trt_cp",
                  grepl("^BRD-", pert), !is.na(NCS))
  P2A_RSVS_NCS <- thp1_all$NCS[thp1_all$pert == BRD_ROSUVASTATIN][1]

  fig3a_df <- thp1_all |>
    dplyr::filter(NCS != 0) |>
    dplyr::select(pert, NCS) |>
    dplyr::mutate(connectivity = factor(
      ifelse(NCS < 0, "Negative connectivity", "Positive connectivity"),
      levels = c("Negative connectivity", "Positive connectivity")))

  n_zero     <- sum(thp1_all$NCS == 0)
  n_nonzero  <- nrow(fig3a_df)
  n_strong   <- sum(fig3a_df$NCS < P2A_THRESHOLD)
  cat("Figure 3A: non-zero THP1 compounds =", n_nonzero,
      "| strong reversers =", n_strong,
      sprintf("(%.1f%%)\n", 100 * n_strong / n_nonzero))

  x_min <- floor(min(fig3a_df$NCS) * 10) / 10
  x_max <- ceiling(max(fig3a_df$NCS) * 10) / 10
  ymax  <- max(hist(fig3a_df$NCS, breaks = P2A_BINS, plot = FALSE)$counts)
  y_top <- ceiling(ymax / 10) * 10

  p_ncs_hist <- ggplot(fig3a_df, aes(x = NCS, fill = connectivity)) +
    geom_histogram(bins = P2A_BINS, color = "white",
                   linewidth = 0.15, alpha = 0.95) +
    # The exact-zero mass is retained as an explicit rug/bar annotation rather
    # than silently dropping 779 compounds from the accounting.
    geom_rug(
      data = data.frame(NCS = rep(0, n_zero)),
      aes(x = NCS),
      sides = "b",
      inherit.aes = FALSE,
      color = PAL$textMute,
      alpha = 0.22,
      length = grid::unit(0.015, "npc")
    ) +
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
         subtitle = paste0(n_nonzero, " non-zero compounds; ", n_zero,
                           " exact-zero compounds retained in accounting | ",
                           n_strong, " strong reversers (NCS < -0.5) | ",
                           "open circle = rosuvastatin")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    THEME_BASE +
    theme(legend.position = "top", legend.justification = "center",
          legend.box.just = "center", panel.grid.minor = element_blank()) +
    coord_cartesian(xlim = c(x_min - 0.05, x_max + 0.05), clip = "off")
}
save_panel(p_ncs_hist, "fig3_panelA_thp1_ncs_zero_aware", 6, 6)

# -----------------------------------------------------------------------------
# FIG 3B — Top-30 THP1 reversal hits, coloured by candidate status
# -----------------------------------------------------------------------------
{
  P3B_N_HITS      <- 30
  P3B_LABEL_SIZE  <- 10
  P3B_SIZE_APPR   <- 5.5
  P3B_SIZE_OTHER  <- 3.2
  P3B_XTITLE_SIZE <- BASE_SIZE - 1
  P3B_XTEXT_SIZE <- BASE_SIZE - 1.5

  fig3b_df <- compound_named |>
    dplyr::filter(!is.na(thp1_NCS), thp1_NCS != 0) |>
    dplyr::arrange(thp1_NCS) |>
    dplyr::slice_head(n = P3B_N_HITS) |>
    dplyr::left_join(
      final_cands |>
        dplyr::select(pert, tier),
      by = "pert"
    ) |>
    dplyr::mutate(
      drug_name = dplyr::case_when(
        is.na(pert_iname) |
          pert_iname == "" |
          pert_iname == "UNIDENTIFIED" ~
          paste0("Unidentified (", pert, ")"),
        TRUE ~ pert_iname
      ),
      category = dplyr::case_when(
        pert %in% c(BRD_ROSUVASTATIN, BRD_ROFLUMILAST) ~ "Highlighted post hoc candidate",
        !is.na(tier) ~ "FDA-approved",
        TRUE ~ "Not FDA-approved"
      ),
      category = factor(category, levels = names(TIER_COLORS)),
      is_approved = category != "Not FDA-approved",
      drug_name = factor(
        drug_name,
        levels = rev(unique(drug_name))
      )
    )

  if (
    nrow(fig3b_df) != 30L ||
    anyDuplicated(as.character(fig3b_df$drug_name))
  ) {
    stop("Figure 2B top-30 compound labels are not uniquely defined.")
  }

  p_thp1 <- ggplot(fig3b_df, aes(x = thp1_NCS, y = drug_name, color = category)) +
    geom_segment(aes(x = 0, xend = thp1_NCS, y = drug_name, yend = drug_name),
                 linewidth = 0.55, color = PAL$rule) +
    geom_point(aes(size = is_approved), alpha = 0.95, show.legend = TRUE) +
    scale_size_manual(values = c("TRUE" = P3B_SIZE_APPR, "FALSE" = P3B_SIZE_OTHER),
                      guide = "none") +
    scale_color_manual(values = TIER_COLORS, name = NULL,
                       limits = names(TIER_COLORS), breaks = names(TIER_COLORS),
                       drop = FALSE) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 3.2, alpha = 1)
      )
    ) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = "Normalized connectivity score in THP1 cells", y = NULL,
         title = "Top 30 reversal compounds \u2014 THP1",
         subtitle = paste0("Ranked by THP1-specific NCS | larger points = ",
                           "FDA-approved | approval status individually verified"),
         caption = paste0("Approval status verified compound-by-compound from ",
                          "regulatory sources; no aggregate MOA-class analysis ",
                          "performed (see Methods).")) +
    THEME_BASE +
    theme(legend.position = "top",
          legend.justification = "center",
          legend.box.just = "center",
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", size = P3B_LABEL_SIZE),
          axis.text.x = element_text(size = P3B_XTEXT_SIZE),
          axis.title.x = element_text(size = P3B_XTITLE_SIZE),
          legend.text = element_text(size = P3B_LABEL_SIZE - 2),
          legend.title = element_text(size = P3B_LABEL_SIZE - 1, face = "bold"))
}
save_panel(p_thp1, "fig3_panelB_thp1_hits", 7, 6)

fig3 <- p_ncs_hist + p_thp1 +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1.25)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(fig3, "fig3_combined", 14, 8)

# -----------------------------------------------------------------------------
# FIG 4 - Candidate robustness and matched-null comparisons (merged 05F)
# -----------------------------------------------------------------------------
build_candidate_panels <- function() {
  check <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
  # EDIT THESE NUMBERS to tune the exports for PowerPoint.
  BASE_SIZE <- 11
  AXIS_TEXT_SIZE <- 9.5
  AXIS_TITLE_SIZE <- 10
  LEGEND_TEXT_SIZE <- 9.5
  TITLE_SIZE <- 12
  SUBTITLE_SIZE <- 9
  AXIS_LINE_WIDTH <- 0.8
  TICK_LINE_WIDTH <- 0.65
  DATA_LABEL_SIZE <- 2.7           # geom_text size in mm; theme sizes are points
  ROOT <- here::here()
  SOURCE_DIR <- dirname(average_tie_path)
  MARKER_PATH <- file.path(ROOT, "data", "processed", "revision_marker_ifn_states.rds")
  inputs <- c(file.path(SOURCE_DIR, c("plot_priority_core.csv", "plot_priority_loo.csv",
                                    "plot_priority_null.csv")), MARKER_PATH)
  check(all(file.exists(inputs)), paste("Missing saved input(s):",
    paste(inputs[!file.exists(inputs)], collapse = "\n")))
  before <- unname(tools::md5sum(inputs))
  core <- utils::read.csv(inputs[1], check.names = FALSE, stringsAsFactors = FALSE)
  loo <- utils::read.csv(inputs[2], check.names = FALSE, stringsAsFactors = FALSE)
  null <- utils::read.csv(inputs[3], check.names = FALSE, stringsAsFactors = FALSE)
  donor <- as.data.frame(marker_state$donor_summary)
  require_columns <- function(x, cols, name) {
    check(all(cols %in% names(x)), paste(name, "is missing:",
      paste(setdiff(cols, names(x)), collapse = ", ")))
  }
  require_columns(core, c("pert", "compound", "query_id", "rank_all", "NCS"), "Core table")
  require_columns(loo, c("pert", "compound", "excluded_patient", "rank_all", "NCS"), "LOO table")
  require_columns(null, c("pert", "compound", "query_id", "null_NCS", "primary_NCS"), "Null table")
  require_columns(donor, c("cohort", "patient", "total_cells", "marker_ifn_high_n",
                           "marker_ifn_high_pct"), "Marker donor table")
  donor <- donor[donor$cohort == "canonical_nine_patient", , drop = FALSE]
  donor$patient <- as.integer(as.character(donor$patient))
  donor <- donor[order(donor$patient), , drop = FALSE]
  DONORS <- as.integer(c(1, 2, 3, 4, 5, 7, 8, 9, 10))
  check(nrow(donor) == 9L && identical(donor$patient, DONORS) &&
          sum(donor$total_cells) == 22197 && sum(donor$marker_ifn_high_n) == 3195 &&
          all(is.finite(donor$marker_ifn_high_pct)),
        "Figure 1D does not match the locked nine-donor 22,197/3,195 population.")
  check(all(donor$total_cells >= 30),
        "A donor has fewer than 30 candidates; review the denominator before hiding the yield legend.")
  check(all(abs(donor$marker_ifn_high_pct -
                  100 * donor$marker_ifn_high_n / donor$total_cells) < 0.011),
        "Saved marker percentages disagree with their numerators/denominators.")
  donor$patient_label <- factor(paste0("P", donor$patient), levels = paste0("P", DONORS))

  QUERY_LABELS <- c(
    PRIMARY_150 = "Primary\n150/150", PRIMARY_50 = "Primary\n50/50",
    PRIMARY_100 = "Primary\n100/100", PRIMARY_200 = "Primary\n200/200",
    C5_ALLREST_MAX150 = "C5 vs\nall rest", C5_MARKERLOW_MAX150 = "C5 vs\nmarker-low",
    MARKER_MAX150 = "Marker-high\nvs low", RES03_ALLREST_MAX150 = "Resolution\n0.3",
    RES07_ALLREST_MAX150 = "Resolution\n0.7", RES10_ALLREST_MAX150 = "Resolution\n1.0",
    P6_INCLUDED_MARKER_MAX150 = "Patient 6\nincluded")
  COMPOUNDS <- c("rosuvastatin", "roflumilast")
  COLORS <- c(rosuvastatin = "#0072B2", roflumilast = "#D55E00")
  DISPLAY <- c(rosuvastatin = "Rosuvastatin", roflumilast = "Roflumilast")
  for (name in c("core", "loo", "null")) {
    d <- get(name)
    check(setequal(tolower(d$compound), COMPOUNDS), paste(name, "has unexpected compounds."))
    d$compound <- factor(tolower(d$compound), levels = COMPOUNDS)
    assign(name, d)
  }
  check(nrow(core) == 22L && setequal(core$query_id, names(QUERY_LABELS)) &&
          !anyDuplicated(core[c("pert", "query_id")]) &&
          all(is.finite(core$rank_all) & core$rank_all >= 1 & core$rank_all <= 1376) &&
          all(is.finite(core$NCS)), "Core table must contain 2 compounds x 11 complete queries.")
  if ("percentile_all" %in% names(core)) {
    check(all(abs(core$percentile_all - (core$rank_all - 1) / 1375) < 1e-7),
          "Saved core ranks and percentiles disagree.")
  }
  check(nrow(loo) == 18L && setequal(as.character(loo$excluded_patient), as.character(DONORS)) &&
          !anyDuplicated(loo[c("pert", "excluded_patient")]) &&
          all(is.finite(loo$rank_all) & loo$rank_all >= 1 & loo$rank_all <= 1376) &&
          all(is.finite(loo$NCS)), "LOO table must contain 2 compounds x 9 omitted donors.")
  check(nrow(null) == 200L && all(table(null$compound) == 100L) &&
          !anyDuplicated(null[c("pert", "query_id")]) &&
          all(is.finite(null$null_NCS)) && all(is.finite(null$primary_NCS)),
        "Null table must contain 100 saved results for each compound.")
  core$query_label <- factor(QUERY_LABELS[core$query_id], levels = unname(QUERY_LABELS))
  loo$excluded_patient <- factor(as.character(loo$excluded_patient), levels = as.character(DONORS))

  # Explicit x/y settings prevent inherited child theme elements hiding the axes.
  publication_theme <- ggplot2::theme_minimal(base_size = BASE_SIZE, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text.x = ggplot2::element_text(colour = "black", face = "bold", size = AXIS_TEXT_SIZE),
      axis.text.y = ggplot2::element_text(colour = "black", face = "bold", size = AXIS_TEXT_SIZE),
      axis.title.x = ggplot2::element_text(colour = "black", face = "bold", size = AXIS_TITLE_SIZE),
      axis.title.y = ggplot2::element_text(colour = "black", face = "bold", size = AXIS_TITLE_SIZE),
      axis.line.x.bottom = ggplot2::element_line(colour = "black", linewidth = AXIS_LINE_WIDTH),
      axis.line.y.left = ggplot2::element_line(colour = "black", linewidth = AXIS_LINE_WIDTH),
      axis.ticks.x.bottom = ggplot2::element_line(colour = "black", linewidth = TICK_LINE_WIDTH),
      axis.ticks.y.left = ggplot2::element_line(colour = "black", linewidth = TICK_LINE_WIDTH),
      axis.ticks.length = grid::unit(2, "mm"),
      panel.grid.major = ggplot2::element_line(colour = "#E1E8F0", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = "black", face = "bold", size = TITLE_SIZE),
      plot.subtitle = ggplot2::element_text(colour = "black", size = SUBTITLE_SIZE),
      plot.caption = ggplot2::element_text(colour = "black", size = SUBTITLE_SIZE),
      plot.tag = ggplot2::element_text(colour = "black", face = "plain", size = 14),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(colour = "black", size = LEGEND_TEXT_SIZE),
      strip.text = ggplot2::element_text(colour = "black", face = "bold", size = AXIS_TITLE_SIZE),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 12, 8, 8))
  color_scale <- function() ggplot2::scale_colour_manual(
    values = COLORS, breaks = COMPOUNDS, labels = DISPLAY, name = NULL)

  p_core <- ggplot2::ggplot(core, ggplot2::aes(query_label, rank_all,
                                            colour = compound, group = compound)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = NCS == 0), size = 2.6, stroke = 0.8) +
    color_scale() +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
      labels = c(`FALSE` = "Nonzero NCS", `TRUE` = "NCS = 0"), name = NULL) +
    ggplot2::scale_y_reverse(breaks = c(1, 100, 300, 500, 700),
                            expand = ggplot2::expansion(mult = c(0.04, 0.06))) +
    ggplot2::labs(x = NULL, y = "THP1 average-tie rank (1,376 compounds)",
      title = "Candidate ranks across core analyses",
      subtitle = "Lower ranks indicate stronger negative connectivity; zeros retain their average tied rank") +
    publication_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

  p_loo_ncs <- ggplot2::ggplot(loo, ggplot2::aes(excluded_patient, NCS,
                                               colour = compound, group = compound)) +
    ggplot2::geom_hline(yintercept = -0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.6) + color_scale() +
    ggplot2::labs(x = "Excluded donor", y = "THP1 NCS",
      title = "Leave-one-donor-out connectivity",
      subtitle = "Fixed state labels; dashed line indicates the operational NCS threshold of -0.5") +
    publication_theme
  rank_breaks <- unique(c(1, pretty(range(loo$rank_all), n = 4)))
  rank_breaks <- rank_breaks[rank_breaks >= 1]
  p_loo_rank <- ggplot2::ggplot(loo, ggplot2::aes(excluded_patient, rank_all,
                                                colour = compound, group = compound)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.6) + color_scale() +
    ggplot2::scale_y_reverse(breaks = rank_breaks) +
    ggplot2::labs(x = "Excluded donor", y = "THP1 average-tie rank (1,376 compounds)",
      title = "Leave-one-donor-out ranks",
      subtitle = "One saved query per excluded donor; original state labels retained") + publication_theme

  # Identical bins and axes for the two null panels; zeros remain in the histogram.
  bin_edges <- seq(min(null$null_NCS), max(null$null_NCS), length.out = 21)
  check(diff(range(bin_edges)) > 0, "Null score range is zero; review the saved table.")
  hist_max <- max(vapply(COMPOUNDS, function(cmp) max(graphics::hist(
    null$null_NCS[null$compound == cmp], breaks = bin_edges,
    right = TRUE, include.lowest = TRUE, plot = FALSE)$counts), numeric(1)))
  xrange <- range(c(null$null_NCS, null$primary_NCS))
  make_null <- function(cmp) {
    d <- null[null$compound == cmp, , drop = FALSE]
    observed <- unique(d$primary_NCS)
    check(length(observed) == 1L, paste("Conflicting primary scores for", cmp))
    d$compound_label <- DISPLAY[[cmp]]
    ggplot2::ggplot(d, ggplot2::aes(null_NCS)) +
      ggplot2::geom_histogram(breaks = bin_edges, closed = "right", fill = "#98A8BE",
                             colour = "white", linewidth = 0.25) +
      ggplot2::geom_vline(xintercept = observed, colour = COLORS[[cmp]], linewidth = 0.85) +
      ggplot2::facet_wrap(~compound_label) +
      ggplot2::scale_x_continuous(limits = xrange + c(-0.04, 0.04) * diff(xrange)) +
      ggplot2::scale_y_continuous(limits = c(0, hist_max * 1.12),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::labs(x = "THP1 NCS", y = "Matched random signatures",
        title = "Expression-matched null comparison",
        subtitle = "Vertical line: observed primary NCS; 100 null queries") + publication_theme
  }
  p_null_rosu <- make_null("rosuvastatin")
  p_null_roflu <- make_null("roflumilast")

  check(identical(before, unname(tools::md5sum(inputs))), "A Figure 4 source changed during plotting.")
  list(core = p_core, loo_ncs = p_loo_ncs, loo_rank = p_loo_rank,
       null_roflumilast = p_null_roflu, null_rosuvastatin = p_null_rosu,
       theme = publication_theme, query_labels = QUERY_LABELS,
       source_dir = SOURCE_DIR, inputs = inputs, input_md5 = before,
       colors = COLORS, display = DISPLAY)
}
robust_panels <- build_candidate_panels()
save_panel(robust_panels$core, "fig4_panelA_core_rank_trajectories", 10, 4.2)
save_panel(robust_panels$loo_ncs, "fig4_panelB_loo_ncs_stability", 10, 3.5)
save_panel(robust_panels$null_roflumilast, "fig4_panelC_roflumilast_matched_null", 5, 3.8)
save_panel(robust_panels$null_rosuvastatin, "fig4_panelD_rosuvastatin_matched_null", 5, 3.8)
fig4 <- patchwork::wrap_plots(robust_panels$core, robust_panels$loo_ncs,
  robust_panels$null_roflumilast, robust_panels$null_rosuvastatin,
  design = "AA\nBB\nCD", heights = c(1.1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme(
    plot.tag = element_text(colour = "black", face = "bold", size = 14)))
save_combined(fig4, "fig4_combined", 10, 11)

# -----------------------------------------------------------------------------
# FIG 5A / 5B — Lead compound cell-line profiles
# -----------------------------------------------------------------------------
p_rosuva <- cellline_panel(BRD_ROSUVASTATIN, "Rosuvastatin",
                           "Gold highlight = THP1 (primary monocytic context)")
p_roflu  <- cellline_panel(BRD_ROFLUMILAST,  "Roflumilast",
                           "Gold highlight = THP1 (primary monocytic context)")

save_panel(p_rosuva, "fig5_panelA_rosuvastatin_celllines", W_CELLLINE, H_CELLLINE)
save_panel(p_roflu,  "fig5_panelB_roflumilast_celllines",  W_CELLLINE, H_CELLLINE)

fig5 <- (p_rosuva + p_roflu +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")) &
  theme(legend.position = "top", legend.justification = "center")
save_combined(fig5, "fig5_combined", 12, 6)

# =============================================================================
# SUPPLEMENTARY FIGURES: S1 through S8 (pathway last)
# =============================================================================

# -----------------------------------------------------------------------------
# S1 — PBMC annotation (A UMAP, B lineage dotplot, C subcluster composition)
# -----------------------------------------------------------------------------
# Section 2 already loaded, validated, annotated, and retained this
# object in memory. Do not read the 122,724-cell cache a second time here.
stopifnot(
  isTRUE(S1_CACHE_READY),
  inherits(pbmc_sce, "SingleCellExperiment"),
  sum(pbmc_sce$cell_type == "monocyte") == 22197L,
  sum(pbmc_sce$cell_type == "other") == 100527L
)

cat(S1_CACHE_NOTE, "\n")
{
  S1_CELLTYPE_COL <- "cell_type"
  S1_CELL_COLORS  <- c("monocyte" = PAL$disease, "other" = PAL$neutral)
  S1_MARKERS <- c("CD14","LYZ","S100A8","S100A9","VCAN","FCN1",
                  "CD3D","CD3E","CD19","MS4A1",
                  "GNLY","NKG7","FCER1A","CST3","FCGR3A","MS4A7")
  S1_PT_SIZE <- 0.25
  W_S1 <- 13; H_S1 <- 9

  if (S1_CACHE_READY)
    {
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
      theme(
        panel.grid = element_blank(),

        axis.text = element_text(
          color = "black",
          face = "bold",
          size = BASE_SIZE - 1.5
        ),

        axis.title = element_text(
          color = "black",
          face = "bold",
          size = BASE_SIZE - 1
        ),

        axis.ticks = element_line(
          color = "black",
          linewidth = 0.65
        ),

        axis.line = element_line(
          color = "black",
          linewidth = 0.8
        ),

        axis.ticks.length = grid::unit(2, "mm")
      )

    # ── B: lineage marker dotplot ──
    markers_present <- intersect(S1_MARKERS, rownames(pbmc_sce))
    sce_subset <- pbmc_sce[markers_present, ]
    seu_all <- Seurat::CreateSeuratObject(
      counts    = SingleCellExperiment::counts(sce_subset),
      meta.data = as.data.frame(SummarizedExperiment::colData(sce_subset))
    ) |> Seurat::NormalizeData(verbose = FALSE)
    rm(sce_subset); gc()

    pS1b <- Seurat::DotPlot(
      seu_all,
      features = markers_present,
      group.by = S1_CELLTYPE_COL,
      cols = c(PAL$neutral, PAL$disease),
      scale = FALSE
    ) +
      labs(
        title = "Canonical lineage marker expression",
        x = NULL,
        y = NULL
      ) +
      THEME_BASE +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          color = "black",
          face = "bold",
          size = BASE_SIZE - 2
        ),

        axis.text.y = element_text(
          color = "black",
          face = "bold",
          size = BASE_SIZE - 1.5
        ),

        axis.ticks = element_line(
          color = "black",
          linewidth = 0.65
        ),

        axis.line = element_line(
          color = "black",
          linewidth = 0.8
        ),

        axis.ticks.length = grid::unit(2, "mm")
      )
    rm(seu_all); gc()

    save_panel(pS1a, "suppfig_S1_panelA_pbmc_umap", 5.5, 5)
    save_panel(pS1b, "suppfig_S1_panelB_marker_dotplot", 7.5, 5)
  } else {
    stop(
      paste0(
        "S1A and S1B were not regenerated. ",
        S1_CACHE_NOTE,
        "\nThe script would otherwise copy the unchanged historical figures."
      ),
      call. = FALSE
    )
  }

  # ── C: subcluster composition per patient ──
  # MoCA labels are intentionally not carried into the figure: the locked
  # evidence set does not support a cognition association.
  patient_levels <- paste0("P", sort(unique(cell_meta$patient)))
  comp_df <- cell_meta |>
    dplyr::group_by(patient, subcluster) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      patient_label = factor(
        paste0("P", patient),
        levels = patient_levels),
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

  save_panel(pS1c, "suppfig_S1_panelC_subcluster_comp", 9,   4.5)

  if (S1_CACHE_READY) {
    figS1 <- (pS1a + pS1b + patchwork::plot_layout(widths = c(1, 1.35))) / pS1c +
      patchwork::plot_layout(heights = c(1, 0.9)) +
      patchwork::plot_annotation(tag_levels = "A")
    save_combined(figS1, "suppfig_S1_combined", W_S1, H_S1)
  } else {
    message("S1C regenerated")
  }
  if (length(S1_REUSE_RECORDS)) readr::write_csv(
    dplyr::bind_rows(S1_REUSE_RECORDS),
    file.path(FIGURES_DIR, "S1_reuse_status.csv"))
}

# -----------------------------------------------------------------------------
# S2 — State-definition, resolution, and Patient 6 sensitivity
# -----------------------------------------------------------------------------
{
  resolution_path <- file.path(
    PROCESSED_DIR,
    "revision_state_resolution_sensitivity.rds"
  )
  if (!file.exists(resolution_path)) {
    resolution_path <- find_revision_file(
      "revision_state_resolution_sensitivity.rds",
      roots = c(PROCESSED_DIR, REVISION_RESULTS_DIR)
    )
  }
  resolution_result <- read_figure_rds(resolution_path)
  resolution_summary <- as.data.frame(resolution_result$resolution_summary)
  if (length(setdiff(c("resolution", "selected_n_cells", "jaccard_with_primary_c5"),
                     colnames(resolution_summary))) > 0L) {
    stop("State-resolution sensitivity object is missing required fields.")
  }

  pS2a <- ggplot(
    resolution_summary,
    aes(x = factor(resolution), y = selected_n_cells)
  ) +
    geom_col(fill = PAL$drug, width = 0.68) +
    geom_hline(yintercept = 940, linetype = "dashed", color = PAL$disease) +
    geom_text(aes(label = selected_n_cells), vjust = -0.35, size = 3) +
    labs(
      x = "Clustering resolution",
      y = "Selected IFN-response-high cells",
      title = "Cells defining each resolution sensitivity contrast",
      subtitle = "Dashed line marks the submitted 940-cell cluster-5 population"
    ) +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank())

  pS2b <- ggplot(
    resolution_summary,
    aes(x = resolution, y = jaccard_with_primary_c5)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = PAL$textMute) +
    geom_line(color = PAL$drug, linewidth = 0.75) +
    geom_point(color = PAL$drug, size = 2.8) +
    scale_y_continuous(
      limits = c(0, 1.05),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Clustering resolution",
      y = "Jaccard overlap with submitted cluster 5",
      title = "Concordance across clustering resolutions",
      subtitle = "The marker-defined state remains the reportable definition"
    ) +
    THEME_BASE

  canonical_row <- marker_state$population_summary |>
    dplyr::filter(cohort == "canonical_nine_patient")
  uniform_row <- marker_state$population_summary |>
    dplyr::filter(cohort == "uniform_ten_patient")
  if (nrow(canonical_row) != 1L || nrow(uniform_row) != 1L) {
    stop("Marker-state population summary is missing canonical or uniform rows.")
  }
  canonical_pop <- data.frame(
    cohort = rep("Canonical\n9 donors", 2),
    state = c("IFN-response-low", "IFN-response-high"),
    n = c(canonical_row$marker_ifn_low_n, canonical_row$marker_ifn_high_n),
    stringsAsFactors = FALSE
  )
  uniform_pop <- data.frame(
    cohort = rep("Uniform reconstruction\n10 donors", 2),
    state = c("IFN-response-low", "IFN-response-high"),
    n = c(uniform_row$marker_ifn_low_n, uniform_row$marker_ifn_high_n),
    stringsAsFactors = FALSE
  )
  p6_assign <- as.data.frame(marker_state$uniform_ten_patient_assignments) |>
    dplyr::filter(as.character(patient) == "6")
  p6_pop <- data.frame(
    cohort = "Patient 6\nonly",
    state = c("IFN-response-low", "IFN-response-high"),
    n = c(
      sum(p6_assign$marker_ifn_state == "marker_ifn_low"),
      sum(p6_assign$marker_ifn_state == "marker_ifn_high")
    ),
    stringsAsFactors = FALSE
  )
  population_plot_df <- dplyr::bind_rows(canonical_pop, uniform_pop, p6_pop) |>
    dplyr::mutate(
      state = factor(state, levels = c("IFN-response-low", "IFN-response-high")),
      cohort = factor(cohort, levels = c(
        "Canonical\n9 donors", "Uniform reconstruction\n10 donors", "Patient 6\nonly"
      ))
    )

  population_totals <- population_plot_df |>
    dplyr::group_by(cohort) |>
    dplyr::summarise(total_cells = sum(n), .groups = "drop")
  pS2c <- ggplot(population_plot_df, aes(x = cohort, y = n, fill = state)) +
    geom_col(width = 0.68) +
    geom_text(data = population_totals,
              aes(x = cohort, y = total_cells,
                  label = format(total_cells, big.mark = ",", trim = TRUE)),
              inherit.aes = FALSE, vjust = -0.45, size = 3, color = "black") +
    scale_fill_manual(
      values = c("IFN-response-low" = PAL$resting, "IFN-response-high" = PAL$ifn),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = NULL,
      y = "Candidate monocytes",
      title = "Patient 6 inclusion and population reconstruction",
      subtitle = "Patient 6 is a sensitivity population, not the canonical analysis"
    ) +
    THEME_BASE +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "top"
    )

  save_panel(pS2a, "suppfig_S2_panelA_resolution_cells", 8, 3.5)
  save_panel(pS2b, "suppfig_S2_panelB_resolution_concordance", 3.5, 4.2)
  save_panel(pS2c, "suppfig_S2_panelC_patient6_reconstruction", 4.5, 4.5)
  figS2 <- (pS2a + pS2b) / pS2c +
    patchwork::plot_layout(heights = c(1, 1.1)) +
    patchwork::plot_annotation(tag_levels = "A")
  save_combined(figS2, "suppfig_S2_state_definition_sensitivity", 11, 9)
  message("Supplementary Figure S2 saved.")
}

# -----------------------------------------------------------------------------
# S3 — Pseudobulk heatmap
# -----------------------------------------------------------------------------
{
  S3_N_UP <- 30; S3_N_DOWN <- 30
  W_S3 <- 11; H_S3 <- 10

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

  sig_up      <- head(inph_signature$up,   S3_N_UP)
  sig_down    <- head(inph_signature$down, S3_N_DOWN)
  sig_present <- intersect(c(sig_up, sig_down), rownames(pb_logcpm))

  heat_scaled <- t(scale(t(as.matrix(pb_logcpm[sig_present, ]))))

  ann_col <- data.frame(State = pb_meta$state, MoCA_group = pb_meta$moca_group,
                        row.names = col_order)
  ann_row <- data.frame(Direction = ifelse(sig_present %in% sig_up, "UP", "DOWN"),
                        row.names = sig_present)
  stopifnot(all(unique(ann_col$MoCA_group) %in% names(ANN_COLORS$MoCA_group)))

  gap_col <- sum(pb_meta$state == "classical")

  heat_dirs <- figure_dirs("suppfig_S3_heatmap")
  heat_variants <- if (SAVE_TITLED) names(heat_dirs) else "clean"
  heat_wh <- SIZE_OVERRIDES[["suppfig_S3_heatmap"]]
  if (!is.null(heat_wh)) {
    stopifnot(is.numeric(heat_wh), length(heat_wh) == 2L, all(is.finite(heat_wh)), all(heat_wh > 0))
    W_S3 <- heat_wh[1]; H_S3 <- heat_wh[2]
  }
  for (variant in heat_variants) {
    dest <- heat_dirs[[variant]]
    for (ext in c("pdf", "png")) {
      out_path <- file.path(dest, paste0("suppfig_S3_heatmap.", ext))
      if (ext == "pdf") {
        open_pdf_device(
          filename = out_path,
          width = W_S3,
          height = H_S3
        )
      } else {
        grDevices::png(
          filename = out_path,
          width = W_S3 * COMBINED_DPI,
          height = H_S3 * COMBINED_DPI,
          res = COMBINED_DPI
        )
      }
      pheatmap::pheatmap(
        heat_scaled, color = HEAT_COLORS,
        annotation_col = ann_col, annotation_row = ann_row,
        annotation_colors = ANN_COLORS,
        cluster_rows = TRUE, cluster_cols = TRUE,
        cutree_cols = 2, cutree_rows = 2,
        show_rownames = TRUE, show_colnames = FALSE,
        fontsize_row = 6.5, fontsize = 9, border_color = NA, gaps_col = gap_col,
        main = if (variant == "with_titles") paste0(
          "Top ", S3_N_UP, " UP + ", S3_N_DOWN,
          " DOWN signature genes across 18 pseudobulk samples\n",
          "Row-scaled log-CPM | 9 iNPH patients \u00d7 2 cell states") else "",
        legend_breaks = c(-2, -1, 0, 1, 2),
        legend_labels = c("-2", "-1", "0", "1", "2\n(scaled\nexpression)"))
      dev.off()
    }
  }
  s3_paths <- unlist(lapply(heat_variants, function(v) {
    file.path(heat_dirs[[v]], paste0("suppfig_S3_heatmap.", c("pdf", "png")))
  }), use.names = FALSE)
  register_exports("suppfig_S3_heatmap", s3_paths)

  s3_info <- file.info(s3_paths)

  if (
    any(
      !file.exists(s3_paths) |
      is.na(s3_info$size) |
      s3_info$size <= 0
    )
  ) {
    stop("Supplementary Figure S3 export verification failed.")
  }

  message("Supplementary Figure S3 saved and verified (titled + clean).")
}

# -----------------------------------------------------------------------------
# S4 — Bulk vs subcluster scoring
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
    data.frame(score = score[seu$sub_res05 == "5"], group = "IFN-response-high cluster 5\n(940 cells)"))
  s3_df$group <- factor(s3_df$group,
                        levels = c("All iNPH monocytes\n(22,197 cells)", "IFN-response-high cluster 5\n(940 cells)"))

  bulk_mean <- round(mean(score), 3)
  ifn_mean  <- round(mean(score[seu$sub_res05 == "5"]), 3)
  cat("S4: bulk mean =", bulk_mean, "| IFN-activated mean =", ifn_mean, "\n")

  pS4 <- ggplot(s3_df, aes(x = group, y = score, fill = group)) +
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

  save_panel(pS4, "suppfig_S4_bulk_vs_subcluster", 5.5, 4.5)
  message("Supplementary Figure S4 saved.")
}

# -----------------------------------------------------------------------------
# S5A — Leave-one-out signature stability
# -----------------------------------------------------------------------------
{
  P4A_THRESHOLD <- 0.70
  P4A_SHAPES <- c("UP genes" = 16, "DOWN genes" = 17)
  P4A_YLIM   <- c(0.58, 0.86)

  patient_order <- loo_df |>
    dplyr::arrange(dplyr::desc(jaccard_mean)) |>
    dplyr::mutate(patient_label = paste0("P", dropped)) |>
    dplyr::pull(patient_label)

  p4a_df <- loo_df |>
    dplyr::mutate(
      patient_label = factor(paste0("P", dropped), levels = patient_order)) |>
    tidyr::pivot_longer(c(jaccard_up, jaccard_dn),
                        names_to = "set", values_to = "jaccard") |>
    dplyr::mutate(set = dplyr::recode(set, jaccard_up = "UP genes",
                                      jaccard_dn = "DOWN genes"))

  p4a_text <- loo_df |>
    dplyr::mutate(patient_label = factor(paste0("P", dropped), levels = patient_order),
                  label = formatC(jaccard_mean, digits = 3, format = "f"))

  p_loo <- ggplot(p4a_df, aes(x = patient_label, y = jaccard, shape = set)) +
    geom_hline(yintercept = P4A_THRESHOLD, linetype = "dashed",
               color = PAL$disease, linewidth = 0.75) +
    annotate("text", x = Inf, y = P4A_THRESHOLD + 0.013,
             label = paste0("Reference = ", P4A_THRESHOLD),
             hjust = 1.05, size = 2.8, color = PAL$disease, fontface = "italic") +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text(data = p4a_text,
              aes(x = patient_label, y = P4A_YLIM[1] + 0.004, label = label),
              color = PAL$textMute, size = 2.6, inherit.aes = FALSE) +
    scale_shape_manual(values = P4A_SHAPES, name = "Gene set") +
    scale_y_continuous(limits = P4A_YLIM,
                       breaks = seq(P4A_YLIM[1], P4A_YLIM[2], 0.05)) +
    labs(x = "Patient removed from analysis",
         y = "Jaccard similarity to full 9-patient signature",
         title = "Leave-one-out signature stability",
         subtitle = paste0("Overall mean: UP = 0.714 | DOWN = 0.733 | ",
                           "Combined = 0.724 | Reference = 0.70")) +
    THEME_BASE +
    theme(legend.position = "right",
          axis.text.x = element_text(size = BASE_SIZE - 1))
}
save_panel(p_loo, "suppfig_S5_panelA_loo_stability", W_LOO, H_LOO)

# -----------------------------------------------------------------------------
# S5B — Per-gene persistence vs DE significance
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
    labs(
      x = expression("Shrunken log"[2]~"fold change (apeglm)"),
      y = expression(-log[10]~"(FDR)"),
      title = "Differential expression: submitted cluster 5 vs classical resting",
      subtitle = paste0(
        "DESeq2 | ~ patient + state | cluster 5 vs cluster 0 | ",
        "1,231 genes at FDR < 0.05 | 9 patients"
      )
    ) +
    THEME_BASE +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(color = "black", face = "bold"))
}
save_panel(p_persist, "suppfig_S5_panelB_gene_persistence", W_PERSIST, H_PERSIST)

figS5 <- p_loo / p_persist + patchwork::plot_layout(heights = c(1, 1.4)) +
  patchwork::plot_annotation(tag_levels = "A")
save_combined(figS5, "suppfig_S5_conditional_signature_stability", 11, 10)

# -----------------------------------------------------------------------------
# S6 — Baricitinib cell-line profile
# -----------------------------------------------------------------------------
{
  bari_rank <- compound_named |>
    dplyr::filter(!is.na(thp1_NCS)) |>
    dplyr::arrange(thp1_NCS) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::filter(pert == BRD_BARICITINIB) |>
    dplyr::pull(rank)

  pS6 <- cellline_panel(
    BRD_BARICITINIB, "Baricitinib",
    subtitle_extra = paste0(
      "Post hoc mechanistically aligned secondary candidate; ranked ",
      bari_rank, " by THP1 NCS and identified outside the prespecified top-30 review."),
    highlight = PAL$purple)

  save_panel(pS6, "suppfig_S6_baricitinib_celllines", W_CELLLINE, H_CELLLINE)
  cat("S6: baricitinib rank =", bari_rank, "\n")
  message("Supplementary Figure S6 saved.")
}

# -----------------------------------------------------------------------------
# S7 — Regenerate individual panels from saved 05B/05C tables
# -----------------------------------------------------------------------------
{
  check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  csf_path <- file.path(PROCESSED_DIR, "revision_csf_ifn_contextualization.rds")
  audit_parent <- file.path(REVISION_RESULTS_DIR, "csf_cell_count_heterogeneity_audit")
  audit_paths <- list.files(audit_parent,
                            pattern = "^revision_csf_cell_count_heterogeneity_audit[.]rds$",
                            recursive = TRUE, full.names = TRUE)
  audit_paths <- audit_paths[file.exists(file.path(dirname(audit_paths), "_SUCCESS.txt"))]
  check(length(audit_paths) > 0L, "No completed 05C audit found for S7.")
  audit_path <- audit_paths[which.max(file.info(audit_paths)$mtime)]
  input_paths_s7 <- c(csf_path, audit_path)
  input_md5_s7 <- unname(tools::md5sum(input_paths_s7))
  check(!anyNA(input_md5_s7), "S7 input checksum unavailable.")
  csf <- read_figure_rds(csf_path)
  audit <- read_figure_rds(audit_path)
  check(identical(csf$analysis_id, "revision_05B_csf_monocyte_ifn_contextualization") &&
          identical(audit$analysis_id, "revision_05C_csf_cell_count_heterogeneity_audit"),
        "Unexpected S7 input objects.")
  check(identical(audit$source_05B_md5, input_md5_s7[[1]]),
        "05C audit does not match the current 05B result.")
  paired <- as.data.frame(csf$paired_donor_summary)
  donor <- as.data.frame(audit$donor_audit)
  sensitivity <- as.data.frame(audit$minimum_cell_sensitivity)

  required_paired <- c(
    "patient", "n_cells_CSF", "median_transfer_score_PBMC",
    "median_transfer_score_CSF"
  )
  required_donor <- c(
    "patient", "n_cells", "marker_ifn_high_pct",
    "marker_ifn_high_wilson_lower_pct",
    "marker_ifn_high_wilson_upper_pct"
  )
  check(
    all(required_paired %in% names(paired)) &&
      all(required_donor %in% names(donor)),
    "Unexpected 05B/05C donor schema."
  )

  paired$patient <- as.character(paired$patient)
  donor$patient <- as.character(donor$patient)
  paired <- paired[order(as.numeric(paired$patient)), , drop = FALSE]
  donor <- donor[order(as.numeric(donor$patient)), , drop = FALSE]
  paired$yield_group <- ifelse(
    paired$n_cells_CSF >= 30L,
    ">=30 CSF candidates",
    "<30 CSF candidates"
  )
  donor$yield_group <- ifelse(
    donor$n_cells >= 30L,
    ">=30 CSF candidates",
    "<30 CSF candidates"
  )

  paired_long <- rbind(
    data.frame(
      patient = paired$patient,
      compartment = "PBMC",
      compartment_order = 1,
      median_transfer_score = paired$median_transfer_score_PBMC,
      n_CSF_candidates = paired$n_cells_CSF,
      yield_group = paired$yield_group,
      stringsAsFactors = FALSE
    ),
    data.frame(
      patient = paired$patient,
      compartment = "Ventricular CSF",
      compartment_order = 2,
      median_transfer_score = paired$median_transfer_score_CSF,
      n_CSF_candidates = paired$n_cells_CSF,
      yield_group = paired$yield_group,
      stringsAsFactors = FALSE
    )
  )
  paired_long$compartment <- factor(
    paired_long$compartment,
    levels = c("PBMC", "Ventricular CSF")
  )
  paired_long$yield_group <- factor(
    paired_long$yield_group,
    levels = c(">=30 CSF candidates", "<30 CSF candidates")
  )
  donor$yield_group <- factor(
    donor$yield_group,
    levels = c(">=30 CSF candidates", "<30 CSF candidates")
  )


  check(nrow(paired) == 9L && nrow(donor) == 10L &&
          !anyDuplicated(paired$patient) && !anyDuplicated(donor$patient),
        "S7 expected nine unique paired donors and ten CSF donors.")
  palette <- c(">=30 CSF candidates" = "#0072B2", "<30 CSF candidates" = "#D55E00")
  base_theme <- THEME_BASE + ggplot2::theme(
    legend.position = "bottom", legend.title = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank())
  panel_a <- ggplot2::ggplot(
    paired_long,
    ggplot2::aes(
      x = compartment,
      y = median_transfer_score,
      group = patient,
      colour = yield_group
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.35) +
    ggplot2::geom_line(alpha = 0.75, linewidth = 0.65) +
    ggplot2::geom_point(size = 2.4) +
    ggrepel::geom_text_repel(
      data = paired_long[paired_long$compartment == "Ventricular CSF", ],
      ggplot2::aes(label = paste0("P", patient)),
      direction = "y",
      hjust = 0,
      nudge_x = 0.04,
      box.padding = 0.15,
      point.padding = 0.1,
      min.segment.length = 0,
      segment.alpha = 0.45,
      seed = 42,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.12, 0.24))) +
    ggplot2::labs(
      title = "Paired transfer of the blood-derived program",
      x = NULL,
      y = "Median transfer score"
    ) +
    base_theme +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 24, 5.5, 5.5)) +
    ggplot2::coord_cartesian(clip = "off")

  panel_b <- ggplot2::ggplot(
    donor,
    ggplot2::aes(
      x = factor(patient, levels = patient),
      y = n_cells,
      fill = yield_group
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 30,
      linetype = "dashed",
      colour = "#555555",
      linewidth = 0.5
    ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_point(ggplot2::aes(colour = yield_group), size = 1.5, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_cells),
      vjust = -0.35,
      size = 3
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_log10(
      breaks = c(1, 10, 30, 100, 1000, 3000),
      labels = c("1", "10", "30", "100", "1,000", "3,000"),
      expand = ggplot2::expansion(mult = c(0.02, 0.14))
    ) +
    ggplot2::labs(
      title = "Unequal CSF candidate-cell recovery",
      x = "Patient",
      y = "CSF candidates"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")

  panel_c <- ggplot2::ggplot(
    donor,
    ggplot2::aes(
      x = factor(patient, levels = patient),
      y = marker_ifn_high_pct,
      colour = yield_group
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = marker_ifn_high_wilson_lower_pct,
        ymax = marker_ifn_high_wilson_upper_pct
      ),
      width = 0.15,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::labs(
      title = "Descriptive IFN-high fractions and uncertainty",
      x = "Patient",
      y = "CSF IFN-high cells"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")


  save_panel(panel_a, "suppfig_S7_panelA_paired_transfer", 8, 2.8)
  save_panel(panel_b, "suppfig_S7_panelB_cell_recovery", 5.5, 4.5)
  save_panel(panel_c, "suppfig_S7_panelC_IFN_fraction", 5.5, 4.5)
  figS7 <- panel_a / (panel_b + panel_c) +
    patchwork::plot_annotation(tag_levels = "A")
  save_combined(figS7, "suppfig_S7_csf_donor_transparency", 12, 9)
  check(identical(input_md5_s7, unname(tools::md5sum(input_paths_s7))),
        "An S7 input changed during figure export.")
  message("S7 individual and combined PNG/PDF exports complete; saved data unchanged.")
}

# -----------------------------------------------------------------------------
# S8 — Proposed mechanism schematic
# -----------------------------------------------------------------------------
{
  FS_DRUG  <- 3.8;  FS_DRUG_S <- 2.7
  FS_TGT   <- 3.2;  FS_INTER  <- 3.0
  FS_CONV  <- 3.8;  FS_CONV_S <- 2.5
  FS_OUT   <- 4.0;  FS_OUT_S  <- 2.6
  FS_NOTE  <- 3
  W_MECH   <- 7;    H_MECH    <- 9

  AR <- grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")

  pS8 <- ggplot() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(4, 4, 4, 4),
          aspect.ratio    = 11.2 / 10) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 11.2), expand = FALSE) +

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
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 6.6, ymax = 7.35,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 2.2, y = 6.975, label = "Isoprenoid synthesis \u2193",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 2.2, y = 6.6, yend = 6.12,
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +
    annotate("rect", xmin = 0.15, xmax = 4.25, ymin = 5.2, ymax = 6.05,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 2.2, y = 5.77,
             label = "RAC1/RHOA membrane anchoring \u2193",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("text", x = 2.2, y = 5.43, label = "JAK activation \u2193",
             size = FS_INTER - 0.2, color = PAL$textMute, hjust = 0.5) +
    annotate("segment", x = 2.2, xend = 4.05, y = 5.2, yend = 4.55,
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +

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
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 6.6, ymax = 7.35,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 7.8, y = 6.975, label = "Intracellular cAMP \u2191",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 7.8, y = 6.6, yend = 6.12,
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +
    annotate("rect", xmin = 5.75, xmax = 9.85, ymin = 5.2, ymax = 6.05,
             fill = PAL$intFill, color = PAL$intBorder, linewidth = 0.45) +
    annotate("text", x = 7.8, y = 5.77, label = "PKA activation \u2191",
             size = FS_INTER, color = PAL$text, hjust = 0.5) +
    annotate("text", x = 7.8, y = 5.43,
             label = "NF-\u03baB \u2193  \u00b7  STAT1 phosphorylation \u2193",
             size = FS_INTER - 0.2, color = PAL$textMute, hjust = 0.5) +
    annotate("segment", x = 7.8, xend = 5.95, y = 5.2, yend = 4.55,
             color = PAL$textSoft, linewidth = 0.7, linetype = "dashed", arrow = AR) +

    # ── CONVERGENCE ──
    annotate("rect", xmin = 3.5, xmax = 6.5, ymin = 3.1, ymax = 4.5,
             fill = PAL$intFill, color = PAL$intBorder) +
    annotate("text", x = 5.0, y = 4.1, label = "JAK / STAT signalling",
             size = FS_CONV, fontface = "bold", color = "black", hjust = 0.5) +
    annotate("text", x = 5.0, y = 3.6,
             label = "STAT1 \u2193  \u00b7  NF-\u03baB \u2193  \u00b7  IRF \u2193",
             size = FS_CONV_S, color = PAL$textMute, hjust = 0.5) +
    annotate("text", x = 5.0, y = 3.22, label = "(predicted shared endpoint)",
             size = FS_CONV_S - 0.3, color = PAL$textMute,
             fontface = "italic", hjust = 0.5) +
    annotate("segment", x = 5.0, xend = 5.0, y = 3.1, yend = 2.62,
             color = PAL$textSoft, linewidth = 1.0, linetype = "dashed", arrow = AR) +

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
                            "Arrows (\u2192) = hypothesized downstream activation  \u00b7  ",
                            "\u2191/\u2193 = relative change  \u00b7  untested"),
             size = FS_NOTE, color = PAL$textMute,
             fontface = "italic", hjust = 0.5)

  save_panel(pS8, "suppfig_S8_mechanism", W_MECH, H_MECH)
  message("Supplementary Figure S8 saved.")
}

# =============================================================================
# SECTION FINAL · FIGURE ORDER MANIFESTS
# =============================================================================
main_order <- c(
  "fig1_combined", "fig1_panelA_umap", "fig1_panelB_isg_violin",
  "fig1_panelC_moca_box", "fig1_panelD_marker_frequency", "fig1_panelE_volcano",
  "fig2_combined", "fig2_panelA_cross_disease_population",
  "fig2_panelB_cross_disease_burden", "fig2_panelC_cross_disease_concordance",
  "fig3_combined", "fig3_panelA_thp1_ncs_zero_aware", "fig3_panelB_thp1_hits",
  "fig4_combined", "fig4_panelA_core_rank_trajectories", "fig4_panelB_loo_ncs_stability",
  "fig4_panelC_roflumilast_matched_null", "fig4_panelD_rosuvastatin_matched_null",
  "fig5_combined", "fig5_panelA_rosuvastatin_celllines", "fig5_panelB_roflumilast_celllines"
)
supp_order <- c(
  "suppfig_S1_combined", "suppfig_S1_panelA_pbmc_umap", "suppfig_S1_panelB_marker_dotplot",
  "suppfig_S1_panelC_subcluster_comp", "suppfig_S2_state_definition_sensitivity",
  "suppfig_S2_panelA_resolution_cells", "suppfig_S2_panelB_resolution_concordance",
  "suppfig_S2_panelC_patient6_reconstruction", "suppfig_S3_heatmap",
  "suppfig_S4_bulk_vs_subcluster", "suppfig_S5_conditional_signature_stability",
  "suppfig_S5_panelA_loo_stability", "suppfig_S5_panelB_gene_persistence",
  "suppfig_S6_baricitinib_celllines", "suppfig_S7_csf_donor_transparency",
  "suppfig_S7_panelA_paired_transfer", "suppfig_S7_panelB_cell_recovery",
  "suppfig_S7_panelC_IFN_fraction", "suppfig_S8_mechanism"
)
manifest <- data.frame(
  order = integer(), group = character(), figure = character(), stem = character(),
  role = character(), stringsAsFactors = FALSE
)
add_order <- function(group, stems, figure_labels, role = "figure") {
  stopifnot(length(stems) == length(figure_labels))
  data.frame(order = seq_along(stems), group = group,
             figure = figure_labels, stem = stems, role = role,
             stringsAsFactors = FALSE)
}
manifest <- dplyr::bind_rows(
  add_order("main", main_order,
            c("Figure 1", paste0("Figure 1", LETTERS[1:5]),
              "Figure 2", paste0("Figure 2", LETTERS[1:3]),
              "Figure 3", paste0("Figure 3", LETTERS[1:2]),
              "Figure 4", paste0("Figure 4", LETTERS[1:4]),
              "Figure 5", paste0("Figure 5", LETTERS[1:2]))),
  add_order("supplementary", supp_order,
            c("Supplementary Figure S1", paste0("Supplementary Figure S1", LETTERS[1:3]),
              "Supplementary Figure S2", paste0("Supplementary Figure S2", LETTERS[1:3]),
              "Supplementary Figure S3", "Supplementary Figure S4",
              "Supplementary Figure S5", paste0("Supplementary Figure S5", LETTERS[1:2]),
              "Supplementary Figure S6", "Supplementary Figure S7",
              paste0("Supplementary Figure S7", LETTERS[1:3]), "Supplementary Figure S8"))
)
manifest$root <- unname(OUTPUT_ROOTS[manifest$group])
manifest$clean_dir <- file.path(manifest$root, "clean")
manifest$titled_dir <- file.path(manifest$root, "with_titles")
readr::write_csv(manifest, file.path(FIGURES_DIR, "figure_order_manifest.csv"))
for (grp in names(OUTPUT_ROOTS)) {
  readr::write_csv(manifest[manifest$group == grp, , drop = FALSE],
                   file.path(OUTPUT_ROOTS[[grp]], "figure_order_manifest.csv"))
}

exports <- if (length(EXPORT_RECORDS)) dplyr::bind_rows(EXPORT_RECORDS) else data.frame()
if (nrow(exports)) {
  readr::write_csv(exports, file.path(FIGURES_DIR, "export_inventory.csv"))
  if (any(!file.exists(exports$path) | file.info(exports$path)$size <= 0)) {
    stop("At least one recorded figure export is missing or empty.")
  }
}
message("FINAL FIGURE EXPORT COMPLETE")
