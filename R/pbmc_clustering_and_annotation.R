# =============================================================================
# PBMC QC, clustering, and monocyte-annotation utilities
# Script: R/pbmc_clustering_and_annotation.R
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
# Provide reusable functions for gene-symbol harmonization, historical PBMC
# quality control, Harmony integration, Louvain clustering, UMAP generation,
# and optional marker-based monocyte-cluster diagnostics.
#
# Scope
# -----
# The PBMC reconstruction is used for Supplementary Figure S1. The canonical
# 22,197-cell classical-monocyte population used for downstream analyses is
# defined by the released canonical monocyte object and is not redefined by
# the marker-based diagnostic helper below.
# =============================================================================


#' Harmonize gene identifiers to gene symbols
#'
#' Replaces row names with the `Symbol` field in `rowData` when available.
#' Missing or empty symbols retain the original identifier, and duplicated
#' symbols are made unique.
#'
#' @param sce SingleCellExperiment object.
#' @return SingleCellExperiment with harmonized row names.

fix_gene_symbols <- function(sce) {
  symbols <- SummarizedExperiment::rowData(sce)$Symbol
  if (is.null(symbols)) {
    logger::log_warn("No 'Symbol' column in rowData; rownames unchanged.")
    return(sce)
  }
  new_names <- ifelse(is.na(symbols) | symbols == "", rownames(sce), symbols)
  new_names <- make.unique(new_names)
  rownames(sce) <- new_names
  logger::log_info("Converted rownames to gene symbols ({sum(symbols == new_names, na.rm = TRUE)} direct matches).")
  sce
}

#' QC filter, normalize, integrate by patient, and cluster
#'
#' Reconstructs the historical PBMC preprocessing workflow used for
#' Supplementary Figure S1. Quality control excludes cells with fewer than
#' `min_genes_per_cell` detected genes or mitochondrial content above
#' `max_pct_mt`. No total-count MAD filter is applied.
#'
#' After quality control, the function performs Seurat normalization,
#' variable-feature selection, scaling, principal-component analysis,
#' Harmony integration by patient, shared-nearest-neighbor graph
#' construction, Louvain clustering, and UMAP embedding.
#'
#' @param sce SingleCellExperiment with `patient` in colData.
#' @param min_genes_per_cell minimum detected genes per cell.
#' @param max_pct_mt maximum mitochondrial percentage per cell.
#' @param resolution Seurat clustering resolution.
#' @param n_pcs number of principal components used for Harmony,
#'   neighborhood construction, and UMAP.
#' @param clustering_algorithm algorithm passed to
#'   Seurat::FindClusters(); 1 corresponds to Louvain.
#' @param clustering_seed random seed passed to
#'   Seurat::FindClusters().
#' @param umap_seed random seed passed to Seurat::RunUMAP().
#' @return SingleCellExperiment containing the reconstructed PBMC
#'   clustering and dimensional-reduction results.
run_qc_and_clustering <- function(
    sce,
    min_genes_per_cell = 200,
    max_pct_mt = 20,
    resolution = 0.5,
    n_pcs = 30,
    clustering_algorithm = 1,
    clustering_seed = 0,
    umap_seed = 42
) {
  logger::log_info("Starting QC + clustering on {ncol(sce)} cells.")

  seu <- Seurat::CreateSeuratObject(
    counts    = SingleCellExperiment::counts(sce),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))
  )

  seu[["percent_mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

  before <- ncol(seu)
  seu <- subset(seu,
                subset = nFeature_RNA >= min_genes_per_cell &
                  percent_mt <= max_pct_mt)
  logger::log_info("QC filter: {before} \u2192 {ncol(seu)} cells ",
                   "(min_genes={min_genes_per_cell}, max_pct_mt={max_pct_mt}).")

  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)

  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop(
      "The harmony package is required for the reported analysis. ",
      "Run renv::restore() before executing this function."
    )
  }

  logger::log_info("Running Harmony integration by patient.")
  seu <- harmony::RunHarmony(seu, group.by.vars = "patient", verbose = FALSE)
  reduction_for_clustering <- "harmony"

  seu <- Seurat::FindNeighbors(
    seu, reduction = reduction_for_clustering,
    dims = seq_len(n_pcs), verbose = FALSE
  )

  seu <- Seurat::FindClusters(
    seu, resolution = resolution,
    algorithm = clustering_algorithm,
    random.seed = clustering_seed,
    verbose = FALSE
  )

  seu <- Seurat::RunUMAP(
    seu,
    reduction = reduction_for_clustering,
    dims = seq_len(n_pcs),
    umap.method = "uwot",
    metric = "cosine",
    seed.use = umap_seed,
    verbose = FALSE
  )

  logger::log_info("Clustering done: {length(unique(seu$seurat_clusters))} clusters.")

  sce_out <- Seurat::as.SingleCellExperiment(seu)
  sce_out
}


#' Identify monocyte-enriched clusters by marker expression
#'
#' Calculates a per-cell score from canonical monocyte markers and summarizes
#' the score by cluster. By default, clusters whose mean marker score is at
#' least `relative_threshold` times the highest cluster score are labeled as
#' monocyte-enriched.
#'
#' This helper is intended as a PBMC-level diagnostic and annotation utility.
#' It does not redefine the canonical 22,197-cell classical-monocyte population
#' used for the manuscript's downstream analyses; that population is defined
#' by the released canonical monocyte object.
#'
#' @param sce clustered SingleCellExperiment returned by
#'   run_qc_and_clustering().
#' @param markers character vector of canonical monocyte marker genes.
#' @param relative_threshold minimum fraction of the highest cluster-level
#'   marker score required for automatic monocyte-cluster assignment.
#' @param clusters_override optional explicit cluster identifier(s) to label
#'   as monocyte-enriched instead of using the relative-threshold rule.
#' @return SingleCellExperiment with `marker_score` and `cell_type` added to
#'   colData and monocyte-cluster diagnostics stored in metadata.
identify_monocytes <- function(sce,
                               markers = c("CD14", "LYZ", "S100A8", "S100A9",
                                           "FCN1", "VCAN", "CD68", "CST3"),
                               relative_threshold = 0.75,
                               clusters_override = NULL) {
  logger::log_info("Identifying monocyte clusters via marker score.")

  available_markers <- intersect(markers, rownames(sce))
  if (length(available_markers) < 3) {
    stop("Fewer than 3 monocyte markers found in dataset. Check gene symbols.")
  }

  # Per-cell normalized expression of marker panel
  if (
    !"logcounts" %in%
      SummarizedExperiment::assayNames(sce)
  ) {
    sce <- scuttle::logNormCounts(sce)
  }

  norm_counts <- SingleCellExperiment::logcounts(sce)
  marker_score <- Matrix::colMeans(norm_counts[available_markers, ])
  sce$marker_score <- marker_score

  # Mean marker score per cluster
  cluster_col <- if ("seurat_clusters" %in% colnames(SummarizedExperiment::colData(sce))) {
    "seurat_clusters"
  } else if ("ident" %in% colnames(SummarizedExperiment::colData(sce))) {
    "ident"
  } else {
    stop("No cluster column found in colData.")
  }

  cluster_scores <- tapply(marker_score,
                           SummarizedExperiment::colData(sce)[[cluster_col]],
                           mean)

  if (!is.null(clusters_override)) {
    unknown_override <- setdiff(
      as.character(clusters_override),
      names(cluster_scores)
    )

    if (length(unknown_override) > 0L) {
      stop(
        "Unknown override cluster(s): ",
        paste(unknown_override, collapse = ", ")
      )
    }

    monocyte_clusters <- as.character(clusters_override)
  } else {
    top_score <- max(cluster_scores)
    score_cutoff <- relative_threshold * top_score
    monocyte_clusters <- names(cluster_scores)[cluster_scores >= score_cutoff]
    if (length(monocyte_clusters) == 0) {
      stop("No clusters passed relative_threshold = ", relative_threshold,
           ". Highest cluster score: ", round(top_score, 2))
    }
  }

  logger::log_info(
    "Cluster marker scores: {paste(names(cluster_scores), round(cluster_scores, 2), sep='=', collapse=', ')}"
  )

  if (is.null(clusters_override)) {
    logger::log_info(
      "Monocyte clusters selected at relative threshold {relative_threshold}: {paste(monocyte_clusters, collapse=', ')}"
    )
  } else {
    logger::log_info(
      "Monocyte clusters selected by explicit override: {paste(monocyte_clusters, collapse=', ')}"
    )
  }

  # Label cells
  sce$cell_type <- ifelse(
    as.character(SummarizedExperiment::colData(sce)[[cluster_col]]) %in% monocyte_clusters,
    "monocyte", "other"
  )
  logger::log_info("Total monocytes identified: {sum(sce$cell_type == 'monocyte')} / {ncol(sce)} cells.")

  # Preserve the evidence used to assign the monocyte cluster
  S4Vectors::metadata(sce)$monocyte_cluster_scores <- cluster_scores
  S4Vectors::metadata(sce)$monocyte_clusters <- monocyte_clusters
  S4Vectors::metadata(sce)$monocyte_relative_threshold <- relative_threshold

  sce
}
