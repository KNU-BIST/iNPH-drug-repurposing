# =============================================================================
# R/clustering.R — Phase 1.5: QC, normalization, clustering, monocyte ID
# =============================================================================

#' Convert Ensembl IDs in rownames to HGNC symbols
#'
#' Cell Ranger output uses Ensembl IDs as rownames with symbols stored in
#' rowData()$Symbol. Downstream tools (Seurat, AUCell, signature search, GSEA)
#' all expect symbols. This function does the conversion safely — handling
#' duplicates and empty/NA symbols by falling back to the Ensembl ID.
#'
#' @param sce a SingleCellExperiment.
#' @return SingleCellExperiment with symbols as rownames.
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

#' QC filter, normalize, integrate by patient, cluster
#'
#' Standard Seurat workflow with patient-level integration via Harmony.
#' Returns a SingleCellExperiment with cluster assignments in colData.
#'
#' @param sce SingleCellExperiment with `patient` column in colData.
#' @param min_genes_per_cell minimum detected genes per cell.
#' @param max_pct_mt maximum percent mitochondrial per cell.
#' @param resolution Seurat clustering resolution (0.5 → ~10-15 clusters).
#' @param n_pcs principal components for clustering.
run_qc_and_clustering <- function(sce,
                                  min_genes_per_cell = 200,
                                  max_pct_mt         = 20,
                                  resolution         = 0.5,
                                  n_pcs              = 30) {
  logger::log_info("Starting QC + clustering on {ncol(sce)} cells.")

  # Convert to Seurat
  seu <- Seurat::CreateSeuratObject(
    counts    = SingleCellExperiment::counts(sce),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))
  )

  # Mitochondrial percentage
  seu[["percent_mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

  # QC filter
  before <- ncol(seu)
  seu <- subset(seu,
                subset = nFeature_RNA >= min_genes_per_cell &
                  percent_mt   <= max_pct_mt)
  logger::log_info("QC filter: {before} → {ncol(seu)} cells ",
                   "(min_genes={min_genes_per_cell}, max_pct_mt={max_pct_mt}).")

  # Standard normalization
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = n_pcs, verbose = FALSE)

  # Integrate across patients (batch correction)
  if (requireNamespace("harmony", quietly = TRUE)) {
    logger::log_info("Running Harmony integration by patient.")
    seu <- harmony::RunHarmony(seu, group.by.vars = "patient", verbose = FALSE)
    reduction_for_clustering <- "harmony"
  } else {
    logger::log_warn("harmony not installed — clustering on raw PCA. ",
                     "Install with: install.packages('harmony')")
    reduction_for_clustering <- "pca"
  }

  # Clustering
  seu <- Seurat::FindNeighbors(seu, reduction = reduction_for_clustering,
                               dims = seq_len(n_pcs), verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = FALSE)
  seu <- Seurat::RunUMAP(seu, reduction = reduction_for_clustering,
                         dims = seq_len(n_pcs), verbose = FALSE)

  logger::log_info("Clustering done: {length(unique(seu$seurat_clusters))} clusters.")

  # Convert back to SCE so downstream pipeline keeps working
  sce_out <- Seurat::as.SingleCellExperiment(seu)
  sce_out
}

#' Identify monocyte-positive clusters by marker gene expression
#'
#' Scores each cluster on a panel of canonical monocyte markers and labels
#' as monocytes any cluster whose score is within `relative_threshold` of
#' the single highest-scoring cluster. This is adaptive to whatever real
#' separation exists in the data, rather than a fixed absolute cutoff -
#' and, since clustering resolution/numbering is not seeded elsewhere in
#' this pipeline, it identifies the monocyte-enriched cluster(s) by score
#' pattern rather than by a specific cluster ID, so it survives a future
#' rebuild landing on a different cluster number for the same biology.
#'
#' @param sce clustered SCE from run_qc_and_clustering().
#' @param markers character vector of monocyte marker genes.
#' @param relative_threshold clusters scoring at least this fraction of the
#'   top cluster's score are labeled monocytes (default 0.75 - i.e. within
#'   25% of the top score joins the same tier; a real gap below this is
#'   treated as a different population).
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
  norm_counts <- SingleCellExperiment::logcounts(sce)
  if (is.null(norm_counts)) {
    # Fall back to raw counts log-normalized inline
    norm_counts <- log1p(SingleCellExperiment::counts(sce))
  }
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
    monocyte_clusters <- as.character(clusters_override)
    logger::log_info("Using override cluster list: {paste(monocyte_clusters, collapse=',')}")
  } else {
    top_score <- max(cluster_scores)
    score_cutoff <- relative_threshold * top_score
    monocyte_clusters <- names(cluster_scores)[cluster_scores >= score_cutoff]
    if (length(monocyte_clusters) == 0) {
      stop("No clusters passed relative_threshold = ", relative_threshold,
           ". Highest cluster score: ", round(top_score, 2))
    }
  }

  logger::log_info("Cluster marker scores: {paste(names(cluster_scores), round(cluster_scores, 2), sep='=', collapse=', ')}")
  logger::log_info("Monocyte clusters (>= {relative_threshold} of top score {round(max(cluster_scores), 2)}): {paste(monocyte_clusters, collapse=', ')}")
  # Label cells
  sce$cell_type <- ifelse(
    as.character(SummarizedExperiment::colData(sce)[[cluster_col]]) %in% monocyte_clusters,
    "monocyte", "other"
  )
  logger::log_info("Total monocytes identified: {sum(sce$cell_type == 'monocyte')} / {ncol(sce)} cells.")

  sce
}
