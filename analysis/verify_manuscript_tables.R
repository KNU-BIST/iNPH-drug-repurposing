# =============================================================================
# MANUSCRIPT TABLE VERIFICATION - recompute every number fresh, compare to
# what is currently sitting in the rebuilt Word tables.
# =============================================================================
#
# PURPOSE: Claude drafted Tables 1-4 across two Word documents this session.
# Table 2 and Table 4's numbers came from your own pasted console output -
# high confidence. Table 1 and Table 3's NCS/WTCS/FDR came from the OLD
# Table 3 already in your uploaded document, cross-checked against partial
# numbers from earlier in the project, but NOT re-derived fresh from a
# script run this session. This script closes that gap: it recomputes every
# checkable number directly from your local processed/ objects, so you can
# compare against the actual docx tables yourself rather than trust a
# second-hand transcription.
#
# WHAT THIS CANNOT CHECK: regulatory status and pharmacologic class are
# external facts, not derived from your data - no script can verify them.
# Check those against primary sources (FDA labels, DailyMed, EMA EPARs)
# directly. See the printed reminder at the end for exactly which cells in
# Table 3 carry the least independent verification.
#
# HOW TO USE: run this, then open Manuscript_v8_draft.docx and
# Supplementary_Tables_v5_draft.docx side by side and check every printed
# "computed" value against the corresponding table cell. Any mismatch means
# the docx needs correcting, not the other way around - this script reads
# directly from your saved objects, which are the ground truth.
# =============================================================================

suppressMessages({library(dplyr); library(readr)})
PROCESSED_DIR <- here::here("data", "processed")
cat("=================== TABLE VERIFICATION ===================\n\n")

# -----------------------------------------------------------------------------
# SECTION A - Table 2 (Supplementary): cross-disease v4 benchmark
# -----------------------------------------------------------------------------
cat("---- TABLE 2: cross-disease v4 benchmark ----\n\n")
bench <- readRDS(file.path(PROCESSED_DIR, "phase2_benchmark_v4.rds"))
A <- bench$engine_A_primary; B <- bench$engine_B_sensitivity

cat("AUCell population-level (docx rows 1-3):\n")
cat(sprintf("  iNPH n=%d median=%.4f | SLE-hi n=%d median=%.4f p=%.3g | HC-hi n=%d median=%.4f p=%.3g\n",
            A$n_pop["inph_c5"], A$med_pop["inph_c5"], A$n_pop["sle_hi"], A$med_pop["sle_hi"], A$p_c5_sle,
            A$n_pop["hc_hi"], A$med_pop["hc_hi"], A$p_c5_hc))
cat("AUCell burden (docx rows 4-6):\n")
cat(sprintf("  iNPH n=%d %.1f%% | SLE n=%d %.1f%% p=%.3g | HC n=%d %.1f%% p=%.3g\n",
            A$burden_n["inph"], A$burden_med["inph"], A$burden_n["sle"], A$burden_med["sle"], A$p_burden_sle,
            A$burden_n["hc"], A$burden_med["hc"], A$p_burden_hc))
cat("AUCell overlap (docx row 7):\n")
cat(sprintf("  %.1f%% of cluster 5 is IFN-high | shared threshold = %.4f\n\n",
            A$overlap["pct_c5_is_ifnhi"], A$thr_shared))

cat("Mean-difference population-level (docx rows 8-10):\n")
cat(sprintf("  iNPH n=%d median=%.4f | SLE-hi n=%d median=%.4f p=%.3g | HC-hi n=%d median=%.4f p=%.3g\n",
            B$n_pop["inph_c5"], B$med_pop["inph_c5"], B$n_pop["sle_hi"], B$med_pop["sle_hi"], B$p_c5_sle,
            B$n_pop["hc_hi"], B$med_pop["hc_hi"], B$p_c5_hc))
cat("Mean-difference burden (docx rows 11-13):\n")
cat(sprintf("  iNPH n=%d %.1f%% | SLE n=%d %.1f%% p=%.3g | HC n=%d %.1f%% p=%.3g\n",
            B$burden_n["inph"], B$burden_med["inph"], B$burden_n["sle"], B$burden_med["sle"], B$p_burden_sle,
            B$burden_n["hc"], B$burden_med["hc"], B$p_burden_hc))
cat("Mean-difference overlap (docx row 14):\n")
cat(sprintf("  %.1f%% of cluster 5 is IFN-high | shared threshold = %.4f\n\n",
            B$overlap["pct_c5_is_ifnhi"], B$thr_shared))

# -----------------------------------------------------------------------------
# SECTION B - Table 4 (Supplementary): Hallmark sensitivity analysis
# -----------------------------------------------------------------------------
cat("---- TABLE 4: Hallmark IFN-alpha sensitivity analysis ----\n\n")
hallmark_result <- readRDS(file.path(PROCESSED_DIR, "hallmark_ifn_positive_control_result.rds"))
hallmark_table <- hallmark_result@result
hallmark_compounds <- readRDS(file.path(PROCESSED_DIR, "hallmark_ifn_positive_control_compounds.rds"))

compoundinfo_beta <- readr::read_tsv(here::here("data", "raw", "compoundinfo_beta.txt"), show_col_types = FALSE)

check_drug <- function(name) {
  ids <- compoundinfo_beta |> dplyr::filter(tolower(cmap_name) == tolower(name)) |>
    dplyr::pull(pert_id) |> unique()
  if (length(ids) == 0) { cat(sprintf("  %-20s NOT FOUND in compoundinfo_beta\n", name)); return(invisible()) }
  for (id in ids) {
    hit <- hallmark_table |> dplyr::filter(pert == id, toupper(cell) == "THP1")
    if (nrow(hit) == 0) cat(sprintf("  %-20s %-16s no THP1 row\n", name, id))
    else cat(sprintf("  %-20s %-16s NCS=%.4f FDR=%.4g\n", name, id, hit$NCS[1], hit$WTCS_FDR[1]))
  }
}

for (d in c("rosuvastatin","roflumilast","baricitinib","tofacitinib","ruxolitinib",
            "upadacitinib","filgotinib","methotrexate","sulfasalazine",
            "leflunomide","hydroxychloroquine")) check_drug(d)
cat("\n")

# -----------------------------------------------------------------------------
# SECTION C - Table 1 + Table 3: LINCS2 compound rankings
# THIS IS THE GAP-CLOSING SECTION - recomputes NCS/WTCS/FDR/rank fresh from
# the raw LINCS2 result table, not from the old docx Table 3.
# -----------------------------------------------------------------------------
cat("---- TABLE 1 / TABLE 3: LINCS2 rankings (fresh from lincs2_gess_result.rds) ----\n\n")

res <- readRDS(file.path(PROCESSED_DIR, "lincs2_gess_result.rds")) |> signatureSearch::result()

# Full ranking, exactly matching 03_lincs_connectivity.R's method
ranked <- res |>
  dplyr::filter(!is.na(NCS), type == "trt_cp", grepl("^BRD-", pert)) |>
  dplyr::group_by(pert) |>
  dplyr::summarise(
    thp1_NCS  = { v <- NCS[toupper(cell) == "THP1"];       if (length(v)) v[1] else NA_real_ },
    thp1_WTCS = { v <- WTCS[toupper(cell) == "THP1"];      if (length(v)) v[1] else NA_real_ },
    thp1_FDR  = { v <- WTCS_FDR[toupper(cell) == "THP1"];  if (length(v)) v[1] else NA_real_ },
    .groups = "drop"
  ) |>
  dplyr::filter(!is.na(thp1_NCS)) |>
  dplyr::arrange(thp1_NCS) |>
  dplyr::mutate(rank = dplyr::row_number())

# Every BRD ID that appears in Table 1 or Table 3, with the docx's claimed
# rank/NCS/WTCS/FDR alongside - PASTE THESE NEXT TO THE ACTUAL TABLE CELLS.
check_ids <- tibble::tribble(
  ~docx_rank, ~pert,           ~docx_compound,
  1,  "BRD-K82941592", "Rosuvastatin",
  2,  "BRD-K43002773", "GDC-0068",
  3,  "BRD-K12683773", "Methyl-fasudil",
  4,  "BRD-K92317137", "BRD-K92317137",
  5,  "BRD-K92303087", "Oglemilast",
  6,  "BRD-A02367930", "Ethinyl estradiol",
  7,  "BRD-K93255255", "Indoximod",
  8,  "BRD-K75641298", "Metoclopramide",
  9,  "BRD-K59013864", "AM-24",
  10, "BRD-K40624912", "ZM-39923",
  11, "BRD-K84925803", "T-5224",
  12, "BRD-K07859598", "GW-5074",
  13, "BRD-K14993104", "Bemegride",
  14, "BRD-K03194791", "Roflumilast",
  15, "BRD-K01436366", "XMD-1150",
  16, "BRD-K96615647", "AV-608",
  17, "BRD-A13084692", "Troglitazone",
  18, "BRD-K00610438", "Altanserin",
  19, "BRD-K60038276", "Irbesartan",
  20, "BRD-A14262390", "Guanadrel",
  21, "BRD-K77987382", "Mebendazole",
  22, "BRD-K82181189", "Ischemin",
  23, "BRD-K59962020", "CHEMBL-374350",
  24, "BRD-K94830329", "Ataluren",
  25, "BRD-K59019422", "BRD-K59019422",
  26, "BRD-K29950728", "Clomifene",
  27, "BRD-K70642949", "GSK-2334470",
  28, "BRD-K01674964", "Histamine",
  29, "BRD-K99749624", "Linifanib",
  30, "BRD-K68548958", "C-646",
  39, "BRD-K53581288", "Baricitinib"
)

verify_df <- check_ids |>
  dplyr::left_join(ranked, by = "pert")

cat(sprintf("%-12s %-20s %-16s %6s %8s %8s %10s %10s\n",
            "docx_rank", "docx_compound", "BRD_ID", "rank", "NCS", "WTCS", "FDR", "MATCH?"))
for (i in seq_len(nrow(verify_df))) {
  r <- verify_df[i, ]
  rank_match <- !is.na(r$rank) && r$rank == r$docx_rank
  cat(sprintf("%-12d %-20s %-16s %6s %8.4f %8.4f %10.3g %10s\n",
              r$docx_rank, r$docx_compound, r$pert,
              ifelse(is.na(r$rank), "NA", r$rank), r$thp1_NCS, r$thp1_WTCS, r$thp1_FDR,
              ifelse(rank_match, "OK", "CHECK")))
}

n_mismatch <- sum(verify_df$rank != verify_df$docx_rank, na.rm = TRUE)
cat(sprintf("\n%d of %d rank mismatches. Any row marked CHECK needs the docx corrected to match this output.\n\n",
            n_mismatch, nrow(verify_df)))

# -----------------------------------------------------------------------------
# SECTION D - What this script CANNOT verify
# -----------------------------------------------------------------------------
cat("---- WHAT WAS NOT INDEPENDENTLY VERIFIED THIS SESSION ----\n\n")
cat("Regulatory status and pharmacologic class are external facts, not in\n")
cat("your data - no script can check them. Confidence levels as drafted:\n\n")
cat("HIGH (sourced with citations, per your and Dr. Nila's review):\n")
cat("  All 8 Table 1 compounds' regulatory status.\n")
cat("  Table 3 ranks 1, 6, 8, 14, 17, 19, 21, 24, 26, 28, 39 (had citations).\n\n")
cat("LOWER (pharmacologic class from general knowledge / old CLUE column,\n")
cat("not independently re-sourced this session) - check these against\n")
cat("primary references before treating as final:\n")
cat("  Ranks 9 (AM-24), 13 (bemegride), 15 (XMD-1150) - class assignments\n")
cat("  are the least certain in the table.\n")
cat("  Ranks 23 (CHEMBL-374350) - marked 'not verified' already, correctly.\n\n")

cat("=================== VERIFICATION COMPLETE ===================\n")
