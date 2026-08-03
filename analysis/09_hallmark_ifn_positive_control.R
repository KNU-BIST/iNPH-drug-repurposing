# =============================================================================
# HALLMARK IFN-ALPHA SENSITIVITY ANALYSIS - replaces the invalid RA query
# =============================================================================
#
# STATUS: per Dr. Nila 2026-07-30/31 - KEEP as supplementary material, framed
# as an "orthogonal interferon-pathway sensitivity analysis," NOT as a
# positive control or independent validation. Report all tested compounds
# transparently. Baricitinib's non-significant result here does NOT
# invalidate its significant recovery in the primary iNPH query (different
# 55-gene up-only query; LINCS connectivity is compound/dose/context-
# dependent) - but the inconsistency must be stated plainly. NEVER write that
# baricitinib was "consistently recovered across independent queries."
#
# WHY THIS REPLACES 09_ra_positive_control.R
# --------------------------------------------
# The archived phase6_positive_control_ra.rds query (19 up / 14 down genes)
# was found circular in a 2026-06-23 session - all 19 UP genes overlapped the
# iNPH signature itself - and did not correspond to any real published RA
# paper despite being captioned "published RA monocyte transcriptomics."
# Retired completely; do not resurrect.
#
# SOURCE (verified by direct fetch, not reconstructed from memory)
# -------------------------------------------------------------------------
# HALLMARK_INTERFERON_ALPHA_RESPONSE, MSigDB systematic name M5911.
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html
# Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P.
# PMID 26771021. License: CC-BY-4.0.
# 97 genes with author-provided Entrez IDs used directly - no separate
# symbol->Entrez conversion step, avoiding the coverage-loss pattern seen
# converting the iNPH signature in script 03.
#
# DESIGN: UP-ONLY QUERY. A Hallmark response signature is one-directional;
# there is no natural "down" companion the way the iNPH signature pairs
# activated-vs-resting.
#
# TWO BUGS FIXED THIS SESSION (both confirmed against live output, not assumed)
# --------------------------------------------------------------------------
# 1. result(hallmark_result) fails ("no applicable method") when gess_lincs()
#    is run with an empty downset - the object is confirmed S4/gessResult
#    with a fully populated @result slot (136,460 rows verified), but the
#    result() GENERIC cannot dispatch on it. Fixed: direct slot access,
#    hallmark_result@result, throughout.
# 2. The original compound-naming join used ONLY the narrow GSE70138 file
#    (data/raw/GSE70138_Broad_LINCS_pert_info.txt.gz), which does not contain
#    methotrexate or hydroxychloroquine by name at all. This silently
#    mislabeled hydroxychloroquine as "not profiled in THP1" when it actually
#    has a real, SIGNIFICANT result (NCS +0.898, FDR 0.00327) - confirmed by
#    checking the broader compoundinfo_beta.txt file and looking up its BRD
#    ID directly against hallmark_table. Fixed: all compound identification
#    below uses compoundinfo_beta.txt (broader coverage), with per-BRD-ID
#    direct lookup against the primary result table rather than a fragile
#    name-based join - so a naming gap can never again silently hide a real
#    result the way it did here.
#
# =============================================================================

library(here); library(dplyr); library(tibble); library(ggplot2)
library(signatureSearch); library(ExperimentHub); library(readr); library(stringr)
library(purrr)

# -----------------------------------------------------------------------------
# Section 0 . Environment (same path-safety pattern as 03)
# -----------------------------------------------------------------------------
ext_root <- Sys.getenv("INPH_EXTERNAL_DATA")
eh_cache <- if (nzchar(ext_root)) file.path(ext_root, "experiment_hub_cache") else
  here::here("data", "raw", "experiment_hub_cache")
dir.create(eh_cache, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(eh_cache)) {
  eh_cache <- here::here("data", "raw", "experiment_hub_cache")
  dir.create(eh_cache, showWarnings = FALSE, recursive = TRUE)
  stopifnot(dir.exists(eh_cache))
}
setExperimentHubOption("CACHE", eh_cache)
PROCESSED_DIR <- here::here("data", "processed")

# -----------------------------------------------------------------------------
# Section 1 . The authoritative 97-gene Hallmark IFN-alpha set
# -----------------------------------------------------------------------------
hallmark_ifn_alpha <- tibble::tribble(
  ~symbol,     ~entrez,
  "ADAR",      "103",     "B2M",       "567",     "BATF2",     "116071",
  "BST2",      "684",     "C1S",       "716",     "CASP1",     "834",
  "CASP8",     "841",     "CCRL2",     "9034",    "CD47",      "961",
  "CD74",      "972",     "CMPK2",     "129607",  "CNP",       "1267",
  "CSF1",      "1435",    "CXCL10",    "3627",    "CXCL11",    "6373",
  "DDX60",     "55601",   "DHX58",     "79132",   "EIF2AK2",   "5610",
  "ELF1",      "1997",    "EPSTI1",    "94240",   "MVB12A",    "93343",
  "TENT5A",    "55603",   "CMTR1",     "23070",   "GBP2",      "2634",
  "GBP4",      "115361",  "GMPR",      "2766",    "HERC6",     "55008",
  "HLA-C",     "3107",    "IFI27",     "3429",    "IFI30",     "10437",
  "IFI35",     "3430",    "IFI44",     "10561",   "IFI44L",    "10964",
  "IFIH1",     "64135",   "IFIT2",     "3433",    "IFIT3",     "3437",
  "IFITM1",    "8519",    "IFITM2",    "10581",   "IFITM3",    "10410",
  "IL15",      "3600",    "IL4R",      "3566",    "IL7",       "3574",
  "IRF1",      "3659",    "IRF2",      "3660",    "IRF7",      "3665",
  "IRF9",      "10379",   "ISG15",     "9636",    "ISG20",     "3669",
  "LAMP3",     "27074",   "LAP3",      "51056",   "LGALS3BP",  "3959",
  "LPAR6",     "10161",   "LY6E",      "4061",    "MOV10",     "4343",
  "MX1",       "4599",    "NCOA7",     "135112",  "NMI",       "9111",
  "NUB1",      "51667",   "OAS1",      "4938",    "OASL",      "8638",
  "OGFR",      "11054",   "PARP12",    "64761",   "PARP14",    "54625",
  "PARP9",     "83666",   "PLSCR1",    "5359",    "PNPT1",     "87178",
  "HELZ2",     "85441",   "PROCR",     "10544",   "PSMA3",     "5684",
  "PSMB8",     "5696",    "PSMB9",     "5698",    "PSME1",     "5720",
  "PSME2",     "5721",    "RIPK2",     "8767",    "RNF31",     "55072",
  "RSAD2",     "91543",   "RTP4",      "64108",   "SAMD9",     "54809",
  "SAMD9L",    "219285",  "SELL",      "6402",    "SLC25A28",  "81894",
  "SP110",     "3431",    "STAT2",     "6773",    "TAP1",      "6890",
  "TDRD7",     "23424",   "TMEM140",   "55281",   "TRAFD1",    "10906",
  "TRIM14",    "9830",    "TRIM21",    "6737",    "TRIM25",    "7706",
  "TRIM26",    "7726",    "TRIM5",     "85363",   "TXNIP",     "10628",
  "UBA7",      "7318",    "UBE2L6",    "9246",    "USP18",     "11274",
  "WARS1",     "7453"
)
stopifnot(nrow(hallmark_ifn_alpha) == 97)
cat("Hallmark IFN-alpha set loaded: 97 genes (Liberzon et al. 2015, M5911)\n\n")

# -----------------------------------------------------------------------------
# Section 2 . Live overlap check against the ACTUAL iNPH signature
# -----------------------------------------------------------------------------
inph_signature <- readRDS(file.path(PROCESSED_DIR, "inph_ifn_signature.rds"))
overlap_genes    <- intersect(hallmark_ifn_alpha$symbol, inph_signature$up)
nonoverlap_genes <- setdiff(hallmark_ifn_alpha$symbol, inph_signature$up)

cat("Overlap with iNPH signature UP set (computed live, not assumed):\n")
cat("  Overlapping (excluded from query):", length(overlap_genes), "of 97 (expect 42)\n")
cat("  Non-overlapping (used in query):  ", length(nonoverlap_genes), "of 97 (expect 55)\n\n")

query_genes <- hallmark_ifn_alpha |> dplyr::filter(symbol %in% nonoverlap_genes)

# -----------------------------------------------------------------------------
# Section 3 . Access LINCS2 + run the query ONCE
# -----------------------------------------------------------------------------
eh          <- ExperimentHub(ask = FALSE)
lincs2_path <- eh[["EH7297"]]

qsig_hallmark <- qSig(
  query = list(upset = query_genes$entrez, downset = character(0)),
  gess_method = "LINCS", refdb = lincs2_path
)

cat("Running gess_lincs() - UP-only query, silent, 5-30 min...\n")
hallmark_result <- gess_lincs(qsig_hallmark, sortby = "NCS", tau = TRUE, workers = 1)
cat("Scoring complete.\n")

# BUG FIX 1: direct slot access, NOT result(). See header note.
hallmark_table <- hallmark_result@result
cat("Total signatures scored:", nrow(hallmark_table), "(expect 136,460)\n\n")

# -----------------------------------------------------------------------------
# Section 4 . Compound ranking - BRD-only, THP1-primary (diagnostic; same
# rationale as the main analysis: THP1 is the only monocytic LINCS line)
# -----------------------------------------------------------------------------
hallmark_compounds <- hallmark_table |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert)) |>
  dplyr::group_by(pert) |>
  dplyr::summarise(
    n_cell_lines = dplyr::n_distinct(cell),
    thp1_NCS     = { v <- NCS[toupper(cell) == "THP1"];      if (length(v)) round(v[1], 4) else NA_real_ },
    thp1_FDR     = { v <- WTCS_FDR[toupper(cell) == "THP1"]; if (length(v)) v[1]           else NA_real_ },
    .groups      = "drop"
  ) |>
  dplyr::arrange(thp1_NCS)

cat("Top 15 Hallmark-query reversal compounds by THP1 NCS (diagnostic):\n")
print(hallmark_compounds |> head(15), n = 15)
cat("\n")

# -----------------------------------------------------------------------------
# Section 5 . Compound identification - BROADER file, per-BRD-ID direct
# lookup against hallmark_table. BUG FIX 2: see header note - this is what
# prevents a naming-join gap from silently hiding a real result again.
# -----------------------------------------------------------------------------
compoundinfo_beta <- readr::read_tsv(
  here::here("data", "raw", "compoundinfo_beta.txt"), show_col_types = FALSE)

lookup_brd_ids <- function(drug_name) {
  compoundinfo_beta |>
    dplyr::filter(tolower(cmap_name) == tolower(drug_name)) |>
    dplyr::pull(pert_id) |>
    unique()
}

thp1_lookup <- function(brd_id) {
  hallmark_table |> dplyr::filter(pert == brd_id, toupper(cell) == "THP1")
}

# -----------------------------------------------------------------------------
# Section 6 . Final verdict table - curated comparator drugs, BRD IDs
# resolved programmatically every run (never hardcoded from a single
# debugging session), one row per (drug, BRD ID) pair actually found.
# -----------------------------------------------------------------------------
comparator_drugs <- tibble::tribble(
  ~Compound,             ~Role,
  "rosuvastatin",        "Primary lead (iNPH query)",
  "roflumilast",         "Primary lead (iNPH query)",
  "baricitinib",         "Secondary candidate (iNPH query, post hoc)",
  "tofacitinib",         "JAK inhibitor (approved)",
  "ruxolitinib",         "JAK inhibitor (approved)",
  "upadacitinib",        "JAK inhibitor (approved)",
  "filgotinib",          "JAK inhibitor (approved)",
  "methotrexate",        "Conventional DMARD",
  "sulfasalazine",       "Conventional DMARD",
  "leflunomide",         "Conventional DMARD",
  "hydroxychloroquine",  "Conventional DMARD"
)

# Column structure matches Dr. Nila's exact spec (2026-07/08-31): compound
# identity, THP1 profiling status, NCS, FDR, direction of connectivity, and
# interpretation as FIVE SEPARATE columns - not a single combined verdict
# string. Absent compounds get ONE neutral status phrase, no reasoning
# attached (her instruction: "state only that they were not identified in
# the available THP1 reference data, without assigning an unverified
# reason" - this deliberately excludes any speculation, e.g. about why
# upadacitinib might be absent from the screening library; that reasoning
# was only ever raised in chat discussion, never in this table, and stays
# out of it here).
NOT_IDENTIFIED <- "Not identified in the available THP1 reference data"

build_verdict_row <- function(Compound, Role) {
  ids <- lookup_brd_ids(Compound)
  if (length(ids) == 0)
    return(tibble::tibble(Compound = Compound, Role = Role, BRD_ID = NA_character_,
                          THP1_Profiling_Status = NOT_IDENTIFIED,
                          NCS = NA_real_, FDR = NA_real_,
                          Direction = NA_character_, Interpretation = NA_character_))

  purrr::map_dfr(ids, function(id) {
    hit <- thp1_lookup(id)
    if (nrow(hit) == 0) {
      tibble::tibble(Compound = Compound, Role = Role, BRD_ID = id,
                     THP1_Profiling_Status = NOT_IDENTIFIED,
                     NCS = NA_real_, FDR = NA_real_,
                     Direction = NA_character_, Interpretation = NA_character_)
    } else {
      ncs <- round(hit$NCS[1], 4); fdr <- hit$WTCS_FDR[1]
      direction <- if (ncs < 0) "Reversal (negative)" else "Activating (positive)"
      interpretation <- dplyr::case_when(
        ncs < -0.5 & fdr < 0.05 ~ "Significant reversal",
        ncs < 0                 ~ "Negative, not significant",
        ncs >= 0 & fdr < 0.05   ~ "Significant positive connectivity",
        TRUE                    ~ "Positive connectivity, not significant"
      )
      tibble::tibble(Compound = Compound, Role = Role, BRD_ID = id,
                     THP1_Profiling_Status = "Profiled",
                     NCS = ncs, FDR = fdr,
                     Direction = direction, Interpretation = interpretation)
    }
  })
}

final_table <- purrr::pmap_dfr(comparator_drugs, build_verdict_row) |>
  dplyr::mutate(Compound = stringr::str_to_sentence(Compound)) |>
  dplyr::distinct()

cat("=== Supplementary Table: Hallmark IFN-alpha sensitivity analysis ===\n")
cat("(orthogonal interferon-pathway check, NOT positive control / validation)\n\n")
print(final_table, n = Inf)
cat("\n")

# -----------------------------------------------------------------------------
# Section 7 . Save
# -----------------------------------------------------------------------------
saveRDS(hallmark_result,    file.path(PROCESSED_DIR, "hallmark_ifn_positive_control_result.rds"))
saveRDS(hallmark_compounds, file.path(PROCESSED_DIR, "hallmark_ifn_positive_control_compounds.rds"))
saveRDS(list(overlap = overlap_genes, nonoverlap = nonoverlap_genes, query_entrez = query_genes),
        file.path(PROCESSED_DIR, "hallmark_overlap_check.rds"))

readr::write_csv(final_table, file.path(PROCESSED_DIR, "hallmark_ifn_alpha_pathway_sensitivity_analysis.csv"))

cat("=== Hallmark IFN-alpha sensitivity analysis complete ===\n")
cat("Citation for Methods/Supplementary Table legend: Liberzon et al. 2015\n")
cat("(PMID 26771021), MSigDB HALLMARK_INTERFERON_ALPHA_RESPONSE (M5911).\n")
cat("55 of 97 genes used (42 overlapping with the iNPH signature excluded);\n")
cat("up-only query; THP1-primary ranking; same WTCS/tau=TRUE methodology as\n")
cat("the primary LINCS2 screen (03_lincs_connectivity.R).\n")
cat("\nFraming per Dr. Nila: orthogonal interferon-pathway sensitivity\n")
cat("analysis, NOT a positive control or validation. Report all compounds\n")
cat("transparently. Do not describe baricitinib as consistently recovered\n")
cat("across independent queries.\n")
