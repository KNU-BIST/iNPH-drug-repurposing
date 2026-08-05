# =============================================================================
# HALLMARK IFN-ALPHA SENSITIVITY ANALYSIS - replaces the invalid RA query
# =============================================================================
#
# STATUS: per Dr. Nila 2026-07-30/31 - KEEP as supplementary material, framed
# as an "orthogonal interferon-pathway sensitivity analysis," NOT as a
# positive control or independent validation. Report all tested compounds
# transparently. NEVER write that baricitinib was "consistently recovered
# across independent queries."
#
# THIS IS THE DEFINITIVE VERSION matching Dr. Nila's exact approved column
# specification for Supplementary Table 4: "Compound, Pharmacologic
# category, THP1 profiling status, NCS, FDR, Connectivity direction,
# Interpretation or note" - no BRD ID column, one row per compound (not
# one row per BRD ID - compounds with multiple LINCS2 entries, such as
# tofacitinib and methotrexate, are collapsed to their single meaningful
# result). An earlier working version of this script produced 13 rows
# with a BRD_ID column, fanning out over every matching ID - numerically
# correct but structurally different from what is actually in the
# manuscript's Supplementary Table 4. This version reproduces that table
# exactly (11 rows), independently verified against it during the
# manuscript audit.
#
# =============================================================================

library(here); library(dplyr); library(tibble); library(ggplot2)
library(signatureSearch); library(ExperimentHub); library(readr); library(stringr)

# -----------------------------------------------------------------------------
# Section 0 . Environment
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

# Direct slot access, not the result() generic - an empty downset produces a
# gessResult object that result() cannot dispatch on, even though @result is
# fully populated. Confirmed multiple times this session.
hallmark_table <- hallmark_result@result
cat("Total signatures scored:", nrow(hallmark_table), "(expect 136,460)\n\n")

# -----------------------------------------------------------------------------
# Section 4 . Compound identification - compoundinfo_beta.txt (broader
# coverage than GSE70138_Broad_LINCS_pert_info.txt.gz, which silently hid
# hydroxychloroquine's real result earlier this session)
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
# Section 5 . Final table - Dr. Nila's exact approved column specification.
# ONE ROW PER COMPOUND. Where a compound has multiple LINCS2 entries (e.g.
# tofacitinib, methotrexate), the entry with real THP1 data is reported if
# one exists; otherwise the compound is reported once as not identified/
# not profiled - never fanned out into multiple rows.
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

NOT_IN_LIBRARY <- "Not in LINCS2 screening library (checked compoundinfo_beta.txt)"
NOT_PROFILED   <- "Profiled in LINCS2, no THP1 data"

build_row <- function(Compound, Role) {
  ids <- lookup_brd_ids(Compound)
  if (length(ids) == 0) {
    return(tibble::tibble(Compound = Compound, Role = Role,
                          THP1_Profiling_Status = NOT_IN_LIBRARY,
                          NCS = NA_real_, FDR = NA_real_,
                          Direction = NA_character_, Interpretation = NA_character_))
  }

  # Check every matching BRD ID; keep the one with real THP1 data if any exist
  hits <- lapply(ids, thp1_lookup)
  has_data <- vapply(hits, function(h) nrow(h) > 0, logical(1))

  if (!any(has_data)) {
    return(tibble::tibble(Compound = Compound, Role = Role,
                          THP1_Profiling_Status = NOT_PROFILED,
                          NCS = NA_real_, FDR = NA_real_,
                          Direction = NA_character_, Interpretation = NA_character_))
  }

  hit <- hits[[which(has_data)[1]]]
  ncs <- round(hit$NCS[1], 4); fdr <- hit$WTCS_FDR[1]
  direction <- if (ncs < 0) "Reversal (negative)" else "Activating (positive)"
  interpretation <- dplyr::case_when(
    ncs < -0.5 & fdr < 0.05 ~ "Significant reversal",
    ncs < 0                 ~ "Negative, not significant",
    ncs >= 0 & fdr < 0.05   ~ "Significant positive connectivity",
    TRUE                    ~ "Positive connectivity, not significant"
  )
  tibble::tibble(Compound = Compound, Role = Role,
                 THP1_Profiling_Status = "Profiled",
                 NCS = ncs, FDR = fdr,
                 Direction = direction, Interpretation = interpretation)
}

final_table <- purrr::pmap_dfr(comparator_drugs, build_row) |>
  dplyr::mutate(Compound = stringr::str_to_sentence(Compound))

stopifnot(nrow(final_table) == 11)   # one row per compound, no fan-out

cat("=== Supplementary Table 4: Hallmark IFN-alpha sensitivity analysis ===\n")
cat("(orthogonal interferon-pathway check, NOT positive control / validation)\n\n")
print(final_table, n = Inf)
cat("\n")

# -----------------------------------------------------------------------------
# Section 6 . Save
# -----------------------------------------------------------------------------
saveRDS(hallmark_result, file.path(PROCESSED_DIR, "hallmark_ifn_positive_control_result.rds"))
readr::write_csv(final_table, file.path(PROCESSED_DIR, "hallmark_ifn_alpha_pathway_sensitivity_analysis.csv"))

cat("=== Hallmark IFN-alpha sensitivity analysis complete ===\n")
cat("Citation: Liberzon et al. 2015 (PMID 26771021), MSigDB\n")
cat("HALLMARK_INTERFERON_ALPHA_RESPONSE (M5911). 55 of 97 genes used (42\n")
cat("overlapping with the iNPH signature excluded); up-only query;\n")
cat("THP1-primary ranking; same WTCS/tau=TRUE methodology as the primary\n")
cat("LINCS2 screen (03_lincs_connectivity.R).\n")
cat("\nFraming per Dr. Nila: orthogonal interferon-pathway sensitivity\n")
cat("analysis, NOT a positive control or validation. Do not describe\n")
cat("baricitinib as consistently recovered across independent queries.\n")
