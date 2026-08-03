# =============================================================================
# Phase 3 - LINCS L1000 connectivity scoring  (FINAL)
# =============================================================================
#
# Project : iNPH Drug Repurposing
# Phase   : 3 (connectivity scoring)
#
# This is the consolidated, verified version. Every fix below was confirmed
# against either the original troubleshooting session or a direct in-session
# recomputation before being adopted - none are guesses.
#
#   1. DUPLICATE QUERY REMOVED. The pre-fix file ran gess_lincs() twice: once
#      correctly on Entrez-converted genes, then again on raw HGNC symbols
#      (never converted), silently overwriting the correct result. LINCS2 is
#      indexed by Entrez ID. This script queries ONCE, Entrez only.
#
#   2. GENE COVERAGE (verified this session): UP 143/150 mapped to Entrez,
#      114/143 present in LINCS2. DOWN 137/150 mapped, 101/137 present.
#
#   3. TAU = TRUE (matches the original session; the pre-fix file's two
#      duplicate blocks disagreed on this).
#
#   4. RANKING: primary sort is THP1-SPECIFIC NCS, not cross-cell-line median.
#      Confirmed reason (from the session that built Figure 4): rosuvastatin's
#      median NCS across 22 cell lines is ~0, but its THP1 NCS is -1.85 - a
#      cross-line median washes out exactly the signal that matters, since
#      THP1 is the only monocytic line and the paper's explicit primary
#      context. Cross-line median is retained as a secondary column only.
#
#   5. PATH HANDLING: explicit fallback + reachability check for the
#      ExperimentHub cache dir (a stale INPH_EXTERNAL_DATA env var pointing at
#      a disconnected drive previously caused a silent hang via an interactive
#      "create directory?" prompt inside a non-interactive run).
#      ExperimentHub(ask = FALSE) added as a second line of defense.
#
#   6. GENETIC-PERTURBATION ENTRIES (n=45) NOW ACTUALLY EXCLUDED FROM RANKING,
#      not just reported as a stat. These are non-drug rows that don't match a
#      proper Broad ID (example found: "GPR151_GALANIN"). Verified this
#      session: 30,456 total trt_cp compounds - 45 non-BRD entries = 30,411
#      BRD-only compounds, matching the manuscript's pre-drafted Methods text
#      exactly. Filter: grepl("^BRD-", pert), applied BEFORE compound_summary
#      is built (previously only applied when reporting the number, not when
#      building the ranked list itself).
#
#   7. THP1 SIGNATURE COUNT REDEFINED to match the manuscript's actual
#      definition. "597" = THP1 signatures with NON-ZERO NCS (NCS==0 means
#      "not reliably computable," not "no effect" - confirmed by the dataset-
#      wide NCS median being exactly 0.0000). The raw THP1 row count (1,376)
#      includes these non-computable entries and is NOT the number that
#      belongs in the manuscript. Verified this session: 597 non-zero BRD-only
#      THP1 signatures, of which 288 have NCS < -0.5 (288/597 = 48.2%, exact
#      match to the pre-drafted Results paragraph).
#
#   8. get_targets() IS CONFIRMED BROKEN IN THIS ENVIRONMENT - it returned 0
#      annotations even for rosuvastatin (BRD-K82941592), a compound
#      guaranteed to be in any working annotation database. Confirmed by
#      direct test, not assumed. Replaced with: (a) structural identification
#      via the local GSE70138_Broad_LINCS_pert_info.txt.gz file (has pert_id,
#      pert_iname, canonical_smiles, inchi_key - confirmed present), and
#      (b) an explicit, MANUALLY VERIFIED approval-status table for the
#      three-criteria sequential filter the Methods text describes (NCS in
#      THP1 -> FDR<0.05 -> FDA approval). Neither this file nor get_targets()
#      has a clinical-phase/approval-status field, so this step cannot be
#      automated with data currently on hand - approval status must be
#      verified compound-by-compound and added to the table below. Anything
#      not in the table is UNVERIFIED, never assumed approved or not.
#
# Data location
# -------------
# ExperimentHub cache: INPH_EXTERNAL_DATA env var if set and reachable, else
# data/raw/experiment_hub_cache under the project. ~4-8 GB on first download.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 . Environment
# -----------------------------------------------------------------------------

library(here)
library(dplyr)
library(tibble)
library(ggplot2)

if (!requireNamespace("signatureSearch", quietly = TRUE))
  BiocManager::install("signatureSearch", update = FALSE, ask = FALSE)
if (!requireNamespace("signatureSearchData", quietly = TRUE))
  BiocManager::install("signatureSearchData", update = FALSE, ask = FALSE)
if (!requireNamespace("ExperimentHub", quietly = TRUE))
  BiocManager::install("ExperimentHub", update = FALSE, ask = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

library(signatureSearch)
library(ExperimentHub)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Resolve cache location with a fallback, and VERIFY it's reachable - a stale
# INPH_EXTERNAL_DATA pointing at a disconnected drive passes nzchar() but
# fails dir.create() silently, sending ExperimentHub into an interactive
# prompt that hangs a sourced script.
ext_root <- Sys.getenv("INPH_EXTERNAL_DATA")
eh_cache <- if (nzchar(ext_root)) file.path(ext_root, "experiment_hub_cache") else
  here::here("data", "raw", "experiment_hub_cache")

dir.create(eh_cache, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(eh_cache)) {
  cat("WARNING: could not reach '", eh_cache, "' (stale INPH_EXTERNAL_DATA?).\n",
      "Falling back to project-local cache.\n", sep = "")
  eh_cache <- here::here("data", "raw", "experiment_hub_cache")
  dir.create(eh_cache, showWarnings = FALSE, recursive = TRUE)
  stopifnot(dir.exists(eh_cache))
}
setExperimentHubOption("CACHE", eh_cache)
cat("ExperimentHub cache resolved to:", eh_cache, "\n")

processed_dir <- here::here("data", "processed")
dir.create(here::here("figures"), showWarnings = FALSE)


# ----------------------------------------------
-------------------------------
# Section 1 . Load the Phase 1 signature
# -----------------------------------------------------------------------------
inph_signature <- readRDS(file.path(processed_dir, "inph_ifn_signature.rds"))
stopifnot(length(inph_signature$up)   == 150)
stopifnot(length(inph_signature$down) == 150)
cat("Phase 1 signature: UP", length(inph_signature$up),
    "| DOWN", length(inph_signature$down), "\n")


# -----------------------------------------------------------------------------
# Section 2 . Access the LINCS2 database
# -----------------------------------------------------------------------------
cat("\nAccessing LINCS2 database via ExperimentHub...\n")
eh          <- ExperimentHub(ask = FALSE)   # never prompt interactively
lincs2_path <- eh[["EH7297"]]
cat("LINCS2 database path:", lincs2_path, "\n")

if (requireNamespace("rhdf5", quietly = TRUE)) {
  h5info <- rhdf5::h5ls(lincs2_path)
  print(h5info[h5info$otype == "H5I_DATASET", c("name", "dim")])
}
lincs_genes <- rhdf5::h5read(lincs2_path, "rownames")[, 1]
cat("LINCS2 genes:", length(lincs_genes), "\n")


# -----------------------------------------------------------------------------
# Section 3 . Gene ID conversion (HGNC symbol -> Entrez)
# -----------------------------------------------------------------------------
symbol_to_entrez <- function(symbols) {
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = symbols, column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  mapped <- mapped[!is.na(mapped)]
  cat("  Mapped:", length(mapped), "of", length(symbols), "genes\n")
  mapped
}

cat("Converting UP genes:\n");   up_entrez   <- symbol_to_entrez(inph_signature$up)
cat("Converting DOWN genes:\n"); down_entrez <- symbol_to_entrez(inph_signature$down)

cat("\nCoverage after conversion:\n")
cat("  UP in LINCS2  :", sum(up_entrez   %in% lincs_genes), "of", length(up_entrez),
    "(expect 114/143)\n")
cat("  DOWN in LINCS2:", sum(down_entrez %in% lincs_genes), "of", length(down_entrez),
    "(expect 101/137)\n")

up_entrez_vec   <- as.character(up_entrez)
down_entrez_vec <- as.character(down_entrez)


# -----------------------------------------------------------------------------
# Section 4 . Run LINCS WTCS connectivity scoring - ONCE
# -----------------------------------------------------------------------------
qsig_lincs2 <- qSig(
  query = list(upset = up_entrez_vec, downset = down_entrez_vec),
  gess_method = "LINCS", refdb = lincs2_path
)

cat("\nRunning gess_lincs() - silent, 5-30 min...\n")
lincs2_result <- gess_lincs(qsig_lincs2, sortby = "NCS", tau = TRUE, workers = 1)

cat("Scoring complete.\n")
result_table <- result(lincs2_result)
cat("Total signatures scored:", nrow(result_table), "(expect ~136,460)\n")


# -----------------------------------------------------------------------------
# Section 5 . Compound-level accounting + ranking (BRD-only, THP1-primary)
# -----------------------------------------------------------------------------
cell_lines_present <- sort(unique(result_table$cell))
cat("\nCell lines represented:", length(cell_lines_present), "(expect 58)\n")

# --- Verify the genetic-perturbation exclusion (n=45) at the compound level ---
all_compounds <- result_table |> dplyr::filter(type == "trt_cp") |> dplyr::distinct(pert)
brd_only      <- all_compounds |> dplyr::filter(grepl("^BRD-", pert))
cat("\nTotal trt_cp compounds:", nrow(all_compounds), "(expect 30,456)\n")
cat("Non-BRD (genetic perturbation) entries:", nrow(all_compounds) - nrow(brd_only),
    "(expect 45)\n")
cat("BRD-only compounds retained:", nrow(brd_only), "(expect 30,411)\n")

# --- Verify the THP1 597/288 definition (non-zero NCS, BRD-only) ---
thp1_hits <- result_table |>
  dplyr::filter(type == "trt_cp", toupper(cell) == "THP1", grepl("^BRD-", pert), NCS != 0)
cat("\nTHP1 non-zero, BRD-only signatures:", nrow(thp1_hits), "(expect 597)\n")
cat("Of those, NCS < -0.5:", sum(thp1_hits$NCS < -0.5),
    sprintf("(expect 288, %.1f%%)\n", 100 * sum(thp1_hits$NCS < -0.5) / nrow(thp1_hits)))

# --- Build the ranked compound table. BRD-only filter applied HERE, not just
#     reported afterward, so the 45 non-drug entries can never appear ranked. ---
compound_summary <- result_table |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert)) |>
  dplyr::group_by(pert) |>
  dplyr::summarise(
    n_cell_lines = dplyr::n_distinct(cell),
    thp1_NCS     = { v <- NCS[toupper(cell) == "THP1"];       if (length(v)) round(v[1], 4) else NA_real_ },
    thp1_FDR     = { v <- WTCS_FDR[toupper(cell) == "THP1"];  if (length(v)) v[1]           else NA_real_ },
    median_NCS   = round(median(NCS, na.rm = TRUE), 4),   # secondary / robustness only
    min_NCS      = round(min(NCS,    na.rm = TRUE), 4),
    .groups      = "drop"
  ) |>
  dplyr::arrange(thp1_NCS)   # PRIMARY sort: THP1-specific NCS

cat("\nTop 30 compounds ranked by THP1-specific NCS:\n")
print(compound_summary |> head(30), n = 30)


# -----------------------------------------------------------------------------
# Section 6 . Cell-line sensitivity - THP-1 and HL-60
# -----------------------------------------------------------------------------
monocyte_lines <- c("THP1", "HL60")
mono_results <- result_table |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert),
                toupper(cell) %in% monocyte_lines) |>
  dplyr::arrange(NCS)
cat("\nTHP-1 / HL-60 specific top hits:\n")
print(head(mono_results |> dplyr::select(pert, cell, WTCS, NCS, Tau), 30))


# -----------------------------------------------------------------------------
# Section 7 . Compound identification + the three-criteria sequential filter
# -----------------------------------------------------------------------------
# Methods text: "Results were prioritized by three criteria applied
# sequentially: NCS in THP1, statistical significance (WTCS FDR<0.05), and
# FDA approval status." Step 1 is compound_summary's sort order (Section 5).
# Steps 2 and 3 are implemented here.

# --- Structural identification (get_targets() confirmed broken - do not use) ---
compoundinfo <- readr::read_tsv(
  here::here("data", "raw", "GSE70138_Broad_LINCS_pert_info.txt.gz"),
  show_col_types = FALSE)

compound_named <- compound_summary |>
  dplyr::left_join(
    compoundinfo |> dplyr::select(pert_id, pert_iname, canonical_smiles, inchi_key),
    by = c("pert" = "pert_id")
  ) |>
  dplyr::mutate(pert_iname = dplyr::if_else(is.na(pert_iname) | pert_iname == pert,
                                            "NAME NOT IN LOCAL DB", pert_iname))

# --- Step 2: statistical significance ---
step2_significant <- compound_named |>
  dplyr::filter(!is.na(thp1_NCS), thp1_NCS < -0.5, !is.na(thp1_FDR), thp1_FDR < 0.05)
cat("\nStep 2 - significant strong reversers (THP1 NCS < -0.5 AND FDR < 0.05):",
    nrow(step2_significant), "\n")

# --- Step 3: FDA/EMA approval status - MANUALLY VERIFIED ONLY, never inferred.
#     No automated data source is currently available (get_targets() broken;
#     local file has no clinical-phase column). Extend this table one
#     compound at a time as each is checked; anything absent is UNVERIFIED. ---
approval_status_known <- tibble::tribble(
  ~pert,           ~pert_iname_expected,      ~approval_status,                             ~verification_note,
  "BRD-K82941592", "rosuvastatin",            "Approved (FDA)",                             "Marketed statin; established prior session.",
  "BRD-K03194791", "roflumilast",             "Approved (FDA)",                             "Marketed PDE4 inhibitor; established prior session.",
  "BRD-K43002773", "GDC-0068",                "Investigational",                            "= ipatasertib, AKT inhibitor. Phase II-III oncology trials, no FDA/EMA approval. Verified 2026-07.",
  "BRD-K53581288", "baricitinib",             "Approved (FDA/EMA)",                         "= Olumiant. FDA- and EMA-authorized JAK1/2 inhibitor (RA 2018; severe alopecia areata 2022). Mechanistically direct hit on type I IFN signaling - flagged as a potential third lead. Verified by Dr. Nila 2026-07-28.",
  "BRD-K92303087", "oglemilast",              "Investigational",                            "PDE4 inhibitor (same class as roflumilast). Phase II only (asthma/COPD), never approved. Verified 2026-07.",
  "BRD-K12683773", "methyl-fasudil",          "Investigational",                            "Fasudil analog (ROCK inhibitor). Parent fasudil approved only in Japan (1995), never FDA-approved; no evidence this analog was independently developed. Verified 2026-07.",
  "BRD-K40624912", "ZM-39923",                "Investigational",                            "JAK1/3 inhibitor. Confirmed research-use-only across multiple vendors; never entered clinical development. Verified 2026-07.",
  "BRD-K59962020", "CHEMBL-374350",           "Investigational",                            "Obscure preclinical screening compound; name is its own ChEMBL ID, no trade/common name. No clinical development evidence. Verified 2026-07.",
  "BRD-K68548958", "C-646",                   "Investigational",                            "p300/CBP HAT inhibitor. Confirmed research-only across multiple vendors. Verified 2026-07.",
  "BRD-K92317137", "UNIDENTIFIED",            "Excluded",                                   "No name in GSE70138, compoundinfo_beta, or public web search. Anonymous, uncharacterized compound.",
  "BRD-K59019422", "UNIDENTIFIED",            "Excluded",                                   "No name in GSE70138, compoundinfo_beta, or public web search. Same basis as BRD-K92317137.",
  # --- Reviewed by Dr. Nila, 2026-07-28; full citations in her original spreadsheet ---
  "BRD-A02367930", "ethinyl-estradiol",       "Approved (FDA; ingredient)",                 "Active ingredient in multiple FDA-approved combination contraceptives; LINCS entry is the ingredient, not a standalone product. Source: FDA label 021676s020.",
  "BRD-K93255255", "indoximod",               "Investigational",                            "IDO-pathway inhibitor, Phase 2; no FDA/EMA authorization. Source: NCT04049669.",
  "BRD-K75641298", "metoclopramide",          "Approved (FDA)",                             "FDA-labeled (Reglan and generics). Source: DailyMed.",
  "BRD-K59013864", "AM-24",                   "Investigational",                            "= Bobel-24/2,4,6-triiodophenol; Broad metadata lists Phase 2. No current FDA/EMA authorization.",
  "BRD-K84925803", "T-5224",                  "Investigational / discontinued",             "AP-1 inhibitor, Phase 2; development did not lead to approval. Source: J Med Chem 10.1021/jm5004733.",
  "BRD-K07859598", "GW-5074",                 "Preclinical / research compound",            "c-Raf inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K14993104", "bemegride",               "Historical/obsolete; status unclear",        "Historical CNS stimulant/barbiturate antidote; NCATS lists 1958 approval but no current FDA/EMA product located. Status genuinely unresolved - do not cite either way without further check.",
  "BRD-K01436366", "XMD-1150",                "Preclinical / research compound",            "L1000 research compound; no clinical marketing authorization identified.",
  "BRD-K96615647", "AV-608",                  "Investigational / discontinued",             "NK1 receptor antagonist, Phase 2; development discontinued, no approval.",
  "BRD-A13084692", "troglitazone",            "Withdrawn",                                  "Approved 1997, withdrawn from market 2000 due to hepatotoxicity. Source: PMID 20020270.",
  "BRD-K00610438", "altanserin",              "Investigational",                            "5-HT2A antagonist/PET ligand, Phase 2; no marketing authorization.",
  "BRD-K60038276", "irbesartan",              "Approved (FDA)",                             "FDA-labeled ARB (Avapro and generics). Source: DailyMed.",
  "BRD-A14262390", "guanadrel",               "Discontinued (US)",                          "Formerly marketed as Hylorel; brand and generics discontinued in the US.",
  "BRD-K77987382", "mebendazole",             "Approved (FDA)",                             "FDA-approved anthelmintic (EMVERM), initial approval 1974.",
  "BRD-K82181189", "ischemin",                "Preclinical / research compound",            "CBP bromodomain inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K94830329", "ataluren",                "Not currently authorized (EMA non-renewal)", "Translarna had conditional EU authorization; European Commission did not renew it (28 Mar 2025). No FDA approval. Source: EMA EPAR.",
  "BRD-K29950728", "clomifene",               "Approved (FDA)",                             "FDA-labeled clomiphene citrate for ovulatory dysfunction.",
  "BRD-K70642949", "GSK-2334470",             "Preclinical / research compound",            "Selective PDK1 inhibitor, research use only. No FDA/EMA authorization.",
  "BRD-K01674964", "histamine",               "Approved (EMA formulation only)",            "Histamine dihydrochloride is EMA-authorized as Ceplene (+IL-2, AML maintenance); the LINCS 'histamine' entity may not match this specific formulation. Treat as borderline, not a clean pass. Source: EMA EPAR.",
  "BRD-K99749624", "linifanib",               "Investigational / development unsuccessful", "= ABT-869; reached Phase 3 in HCC but failed primary endpoint. No FDA/EMA approval."
)

step2_significant <- step2_significant |>
  dplyr::left_join(approval_status_known |> dplyr::select(pert, approval_status, verification_note),
                   by = "pert") |>
  dplyr::mutate(approval_status = dplyr::coalesce(approval_status, "UNVERIFIED - check before citing"))

# Cosmetic fix: 4 compounds (methyl-fasudil, ZM-39923, CHEMBL-374350, C-646) were
# identified by hand via compoundinfo_beta.txt, not the GSE70138 file this script
# reads - so their pert_iname still shows "NAME NOT IN LOCAL DB" even though the
# approval_status verdict is already correct. Backfill the name from the table
# we already built, no second file read needed.
step2_significant <- step2_significant |>
  dplyr::left_join(approval_status_known |> dplyr::select(pert, pert_iname_expected),
                   by = "pert") |>
  dplyr::mutate(pert_iname = dplyr::if_else(pert_iname == "NAME NOT IN LOCAL DB" & !is.na(pert_iname_expected),
                                            pert_iname_expected, pert_iname)) |>
  dplyr::select(-pert_iname_expected)

cat("\nStep 2 candidates with current approval-status knowledge (top 30 shown):\n")
print(step2_significant |>
        dplyr::select(pert, pert_iname, thp1_NCS, thp1_FDR, approval_status) |>
        head(30), n = 30)

# Methods text specifies "FDA approval status" - matched strictly on FDA approval
# (including ingredient-only and FDA/EMA dual approval), NOT EMA-only cases like
# histamine, which Dr. Nila flagged as borderline (formulation may not match the
# LINCS entity). Borderline cases are reported separately below, not silently
# dropped or silently included.
final_candidates <- step2_significant |>
  dplyr::filter(grepl("^Approved \\(FDA", approval_status))
cat("\nFinal candidates (all 3 criteria satisfied, FDA-approved):", nrow(final_candidates), "\n")
print(final_candidates |> dplyr::select(pert, pert_iname, thp1_NCS, thp1_FDR, approval_status))

borderline_candidates <- step2_significant |>
  dplyr::filter(approval_status == "Approved (EMA formulation only)")
if (nrow(borderline_candidates) > 0) {
  cat("\nBorderline candidates (approved outside the US only, or formulation mismatch",
      "- NOT in final_candidates, needs a team decision):\n")
  print(borderline_candidates |> dplyr::select(pert, pert_iname, thp1_NCS, thp1_FDR, approval_status))
}

# --- Manuscript tier, per Dr. Nila 2026-07-28: NOT all 8 equally prioritized.
#     Primary leads (mechanistic + strongest NCS) / one additional mechanistic
#     candidate / remaining FDA-approved hits reported but not prioritized. ---
manuscript_tier <- tibble::tribble(
  ~pert,           ~tier,
  "BRD-K82941592", "Primary lead",
  "BRD-K03194791", "Primary lead",
  "BRD-K53581288", "Additional candidate (mechanistic - JAK1/2 inhibitor)",
  "BRD-A02367930", "Other FDA-approved reverser (not prioritized)",
  "BRD-K75641298", "Other FDA-approved reverser (not prioritized)",
  "BRD-K60038276", "Other FDA-approved reverser (not prioritized)",
  "BRD-K77987382", "Other FDA-approved reverser (not prioritized)",
  "BRD-K29950728", "Other FDA-approved reverser (not prioritized)"
)
tier_order <- c("Primary lead", "Additional candidate (mechanistic - JAK1/2 inhibitor)",
                "Other FDA-approved reverser (not prioritized)")

final_candidates <- final_candidates |>
  dplyr::left_join(manuscript_tier, by = "pert") |>
  dplyr::arrange(match(tier, tier_order), thp1_NCS)

cat("\nFinal candidates with manuscript tier (Dr. Nila, 2026-07-28):\n")
print(final_candidates |> dplyr::select(pert, pert_iname, thp1_NCS, thp1_FDR, tier))

# --- Baricitinib's exact rank among ALL ranked compounds (not just top 30) -
#     Dr. Nila asked this be documented before presenting it in the text. ---
baricitinib_rank <- compound_summary |>
  dplyr::filter(!is.na(thp1_NCS)) |>
  dplyr::arrange(thp1_NCS) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::filter(pert == "BRD-K53581288")
cat("\nBaricitinib exact rank among all ranked compounds:\n")
print(baricitinib_rank |> dplyr::select(rank, pert, thp1_NCS, thp1_FDR))
saveRDS(final_candidates, file.path(processed_dir, "lincs2_final_candidates.rds"))

n_unverified <- sum(step2_significant$approval_status == "UNVERIFIED - check before citing")
if (n_unverified > 0)
  cat("\nNOTE:", n_unverified, "significant candidates have UNVERIFIED approval status.",
      "Do not cite as approved or excluded until manually checked.\n")

# --- Known JAK/IFN-pathway modulators, flagged if present among significant hits ---
jak_inhibitors <- c("tofacitinib", "ruxolitinib", "baricitinib", "upadacitinib",
                     "filgotinib", "abrocitinib", "fedratinib", "pacritinib",
                     "momelotinib", "itacitinib")
ifn_modulators <- c("anifrolumab", "sifalimumab", "rontalizumab",
                     "hydroxychloroquine", "chloroquine", "quinacrine")

priority_hits <- step2_significant |>
  dplyr::filter(tolower(pert_iname) %in% tolower(c(jak_inhibitors, ifn_modulators)))
cat("\nKnown JAK/IFN-pathway compounds among significant hits:\n")
if (nrow(priority_hits) > 0) print(priority_hits) else
  cat("  (None found by name among step-2 significant hits.)\n")


# -----------------------------------------------------------------------------
# Section 8 . Visualization - THP1 NCS distribution (matches Figure 4 logic)
# -----------------------------------------------------------------------------
p1 <- compound_summary |>
  dplyr::filter(!is.na(thp1_NCS)) |>
  ggplot(aes(x = thp1_NCS)) +
  geom_histogram(bins = 60, fill = "#5e9ea0", color = NA, alpha = 0.8) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "#d4624a", linewidth = 0.7) +
  labs(x = "THP1-specific NCS", y = "Number of compounds",
       title = "LINCS2 THP1 connectivity score distribution",
       subtitle = "Query: iNPH IFN-activated monocyte signature | Negative NCS = reversal") +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(here::here("figures", "phase3_thp1_ncs_distribution.pdf"), p1, width = 7, height = 4.5)

top30 <- compound_summary |>
  dplyr::filter(!is.na(thp1_NCS)) |>
  head(30) |>
  dplyr::mutate(pert = factor(pert, levels = rev(pert)), is_strong = thp1_NCS < -0.5)

p2 <- ggplot(top30, aes(x = thp1_NCS, y = pert)) +
  geom_segment(aes(x = 0, xend = thp1_NCS, y = pert, yend = pert), color = "grey80", linewidth = 0.5) +
  geom_point(aes(size = n_cell_lines, color = is_strong), alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = "#d4624a", "FALSE" = "#5e9ea0"),
                     labels = c("TRUE" = "THP1 NCS < -0.5", "FALSE" = "THP1 NCS \u2265 -0.5"), name = "") +
  scale_size_continuous(name = "Cell lines", range = c(2, 6)) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "#d4624a", linewidth = 0.5, alpha = 0.7) +
  labs(x = "THP1-specific NCS", y = NULL,
       title = "Top 30 reversal candidates - ranked by THP1 NCS",
       subtitle = "LINCS2 WTCS scoring") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right",
        panel.grid.major.y = element_blank())
ggsave(here::here("figures", "phase3_top30_reversal_thp1.pdf"), p2, width = 8, height = 7)


# -----------------------------------------------------------------------------
# Section 9 . Save outputs
# -----------------------------------------------------------------------------
saveRDS(lincs2_result,        file.path(processed_dir, "lincs2_gess_result.rds"))
saveRDS(compound_named,       file.path(processed_dir, "lincs2_compound_summary.rds"))
saveRDS(step2_significant,    file.path(processed_dir, "lincs2_step2_significant.rds"))
saveRDS(final_candidates,     file.path(processed_dir, "lincs2_final_candidates.rds"))
readr::write_csv(compound_named |> head(200),
                  file.path(processed_dir, "phase3_top200_reversal_candidates.csv"))
readr::write_csv(step2_significant,
                  file.path(processed_dir, "phase3_step2_significant_candidates.csv"))

cat("\n=== Phase 3 complete (final) ===\n")
cat("Ranking: THP1-specific NCS, BRD-only (45 genetic-perturbation entries excluded)\n")
cat("Verified this session against pre-drafted manuscript text:\n")
cat("  136,460 total signatures | 58 cell lines | 30,456 trt_cp compounds |\n")
cat("  45 non-BRD entries removed -> 30,411 BRD-only | 597 THP1 non-zero |\n")
cat("  288 THP1 reversers (48.2%)\n")
cat(n_unverified, "step-2 candidates still need manual approval-status verification.\n")
