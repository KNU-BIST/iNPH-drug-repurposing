# =============================================================================
# Manuscript-table verification report
# Script: 09_verify_manuscript_tables.R
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
# Verify the computationally derived values used in the manuscript and
# supplementary tables against the finalized analysis outputs.
#
# Scope
# -----
# - cross-disease contextualization values used in Supplementary Table 2;
# - Hallmark IFN-alpha sensitivity values used in Supplementary Table 4;
# - LINCS2 ranks, NCS, WTCS, and FDR values used in Tables 1 and 3;
# - consistency of the full 597-compound export with the underlying LINCS2
#   result object;
# - membership and metrics of the eight reviewed FDA-approved reversers;
# - regulatory-review boundary for the 138 significant strong reversers.
#
# Important
# ---------
# Regulatory status and pharmacologic class are external annotations. This
# script verifies that the saved computational tables contain the finalized
# annotations, but it does not independently establish regulatory truth.
# Those annotations must remain traceable to the cited external sources.
#
# This script is intentionally the final numerical verification step. It does
# not regenerate any upstream analysis.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 . Environment, helpers, and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "readr",
  "tibble",
  "signatureSearch"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
  library(readr)
  library(tibble)
})

PROCESSED_DIR <- here::here(
  "data",
  "processed"
)

input_paths <- c(
  contextualization =
    file.path(PROCESSED_DIR, "phase2_contextualization.rds"),
  hallmark_result =
    file.path(PROCESSED_DIR, "hallmark_ifn_sensitivity_result.rds"),
  hallmark_table =
    file.path(
      PROCESSED_DIR,
      "hallmark_ifn_alpha_pathway_sensitivity_analysis.csv"
    ),
  lincs_result =
    file.path(PROCESSED_DIR, "lincs2_gess_result.rds"),
  compound_summary =
    file.path(PROCESSED_DIR, "lincs2_compound_summary.rds"),
  significant_reversers =
    file.path(PROCESSED_DIR, "lincs2_significant_thp1_reversers.rds"),
  final_candidates =
    file.path(PROCESSED_DIR, "lincs2_final_candidates.rds"),
  full_thp1_export =
    file.path(PROCESSED_DIR, "lincs2_thp1_nonzero_597_compounds.csv")
)

input_info <- file.info(
  unname(input_paths)
)

bad_inputs <- names(input_paths)[
  !file.exists(unname(input_paths)) |
    is.na(input_info$size) |
    input_info$size <= 0
]

if (length(bad_inputs) > 0L) {
  stop(
    "Required verification input(s) missing or empty: ",
    paste(bad_inputs, collapse = ", ")
  )
}

assert_equal <- function(actual, expected, label) {
  if (
    length(actual) != 1L ||
    is.na(actual) ||
    actual != expected
  ) {
    stop(
      label,
      " mismatch: expected ",
      expected,
      ", observed ",
      paste(actual, collapse = ", ")
    )
  }

  invisible(TRUE)
}

assert_close <- function(actual, expected, tolerance, label) {
  if (
    length(actual) != 1L ||
    is.na(actual) ||
    abs(actual - expected) > tolerance
  ) {
    stop(
      label,
      " mismatch: expected ",
      expected,
      " +/- ",
      tolerance,
      ", observed ",
      paste(actual, collapse = ", ")
    )
  }

  invisible(TRUE)
}

assert_true <- function(condition, label) {
  if (
    length(condition) != 1L ||
    is.na(condition) ||
    !isTRUE(condition)
  ) {
    stop(
      label,
      " failed."
    )
  }

  invisible(TRUE)
}

cat(
  "=================== MANUSCRIPT TABLE VERIFICATION ===================\n\n"
)


# -----------------------------------------------------------------------------
# Section 1 . Supplementary Table 2: cross-disease contextualization
# -----------------------------------------------------------------------------

cat(
  "---- SUPPLEMENTARY TABLE 2: CROSS-DISEASE CONTEXTUALIZATION ----\n\n"
)

contextualization <- readRDS(
  input_paths[["contextualization"]]
)

A <- contextualization$engine_A_primary
B <- contextualization$engine_B_sensitivity

if (
  is.null(A) ||
  is.null(B)
) {
  stop(
    "phase2_contextualization.rds is missing the primary or sensitivity engine."
  )
}

cat("AUCell primary engine - population-level activity:\n")
cat(
  sprintf(
    "  iNPH-c5 n=%d median=%.4f | SLE-high n=%d median=%.4f p=%.3g | HC-high n=%d median=%.4f p=%.3g\n",
    A$n_pop["inph_c5"],
    A$med_pop["inph_c5"],
    A$n_pop["sle_hi"],
    A$med_pop["sle_hi"],
    A$p_c5_sle,
    A$n_pop["hc_hi"],
    A$med_pop["hc_hi"],
    A$p_c5_hc
  )
)

cat("AUCell primary engine - donor-level IFN-high burden:\n")
cat(
  sprintf(
    "  iNPH n=%d %.1f%% | SLE n=%d %.1f%% p=%.3g | HC n=%d %.1f%% p=%.3g\n",
    A$burden_n["inph"],
    A$burden_med["inph"],
    A$burden_n["sle"],
    A$burden_med["sle"],
    A$p_burden_sle,
    A$burden_n["hc"],
    A$burden_med["hc"],
    A$p_burden_hc
  )
)

cat(
  sprintf(
    "AUCell cluster-5 concordance: %.1f%% | shared threshold %.4f\n\n",
    A$overlap["pct_c5_is_ifnhi"],
    A$thr_shared
  )
)

cat("Mean-difference sensitivity engine - population-level activity:\n")
cat(
  sprintf(
    "  iNPH-c5 n=%d median=%.4f | SLE-high n=%d median=%.4f p=%.3g | HC-high n=%d median=%.4f p=%.3g\n",
    B$n_pop["inph_c5"],
    B$med_pop["inph_c5"],
    B$n_pop["sle_hi"],
    B$med_pop["sle_hi"],
    B$p_c5_sle,
    B$n_pop["hc_hi"],
    B$med_pop["hc_hi"],
    B$p_c5_hc
  )
)

cat("Mean-difference sensitivity engine - donor-level IFN-high burden:\n")
cat(
  sprintf(
    "  iNPH n=%d %.1f%% | SLE n=%d %.1f%% p=%.3g | HC n=%d %.1f%% p=%.3g\n",
    B$burden_n["inph"],
    B$burden_med["inph"],
    B$burden_n["sle"],
    B$burden_med["sle"],
    B$p_burden_sle,
    B$burden_n["hc"],
    B$burden_med["hc"],
    B$p_burden_hc
  )
)

cat(
  sprintf(
    "Mean-difference cluster-5 concordance: %.1f%% | shared threshold %.4f\n\n",
    B$overlap["pct_c5_is_ifnhi"],
    B$thr_shared
  )
)

# Locked manuscript-level values.
assert_close(
  unname(A$med_pop["inph_c5"]),
  0.0307,
  0.0002,
  "Table 2 AUCell iNPH-c5 median"
)
assert_close(
  unname(A$med_pop["sle_hi"]),
  0.0702,
  0.0002,
  "Table 2 AUCell SLE-high median"
)
assert_close(
  unname(A$med_pop["hc_hi"]),
  0.0252,
  0.0002,
  "Table 2 AUCell HC-high median"
)
assert_close(
  A$p_c5_sle,
  6.42e-06,
  2e-07,
  "Table 2 AUCell iNPH-c5 vs SLE-high p"
)
assert_close(
  A$p_c5_hc,
  0.213,
  0.003,
  "Table 2 AUCell iNPH-c5 vs HC-high p"
)
assert_close(
  unname(A$burden_med["hc"]),
  25.0,
  0.1,
  "Table 2 AUCell HC burden"
)
assert_close(
  unname(A$burden_med["sle"]),
  91.5,
  0.1,
  "Table 2 AUCell SLE burden"
)
assert_close(
  unname(A$burden_med["inph"]),
  13.5,
  0.1,
  "Table 2 AUCell iNPH burden"
)
assert_close(
  A$p_burden_sle,
  1.29e-05,
  3e-07,
  "Table 2 AUCell iNPH vs SLE burden p"
)
assert_close(
  A$p_burden_hc,
  0.0128,
  0.0003,
  "Table 2 AUCell iNPH vs HC burden p"
)
assert_close(
  unname(A$overlap["pct_c5_is_ifnhi"]),
  95.5,
  0.1,
  "Table 2 AUCell cluster-5 concordance"
)
assert_close(
  A$thr_shared,
  0.0973,
  0.0002,
  "Table 2 AUCell shared threshold"
)

assert_close(
  unname(B$med_pop["inph_c5"]),
  0.0725,
  0.0003,
  "Table 2 mean-difference iNPH-c5 median"
)
assert_close(
  unname(B$med_pop["sle_hi"]),
  0.1886,
  0.0003,
  "Table 2 mean-difference SLE-high median"
)
assert_close(
  unname(B$med_pop["hc_hi"]),
  0.0599,
  0.0003,
  "Table 2 mean-difference HC-high median"
)
assert_close(
  B$p_c5_sle,
  1.48e-06,
  1e-07,
  "Table 2 mean-difference iNPH-c5 vs SLE-high p"
)
assert_close(
  B$p_c5_hc,
  0.356,
  0.004,
  "Table 2 mean-difference iNPH-c5 vs HC-high p"
)
assert_close(
  unname(B$burden_med["hc"]),
  25.0,
  0.1,
  "Table 2 mean-difference HC burden"
)
assert_close(
  unname(B$burden_med["sle"]),
  91.8,
  0.1,
  "Table 2 mean-difference SLE burden"
)
assert_close(
  unname(B$burden_med["inph"]),
  30.1,
  0.1,
  "Table 2 mean-difference iNPH burden"
)
assert_close(
  B$p_burden_sle,
  0.00098,
  0.00003,
  "Table 2 mean-difference iNPH vs SLE burden p"
)
assert_close(
  B$p_burden_hc,
  0.491,
  0.004,
  "Table 2 mean-difference iNPH vs HC burden p"
)
assert_close(
  unname(B$overlap["pct_c5_is_ifnhi"]),
  99.5,
  0.1,
  "Table 2 mean-difference cluster-5 concordance"
)
assert_close(
  B$thr_shared,
  0.2517,
  0.0002,
  "Table 2 mean-difference shared threshold"
)


# -----------------------------------------------------------------------------
# Section 2 . Supplementary Table 4: Hallmark IFN-alpha sensitivity analysis
# -----------------------------------------------------------------------------

cat(
  "---- SUPPLEMENTARY TABLE 4: HALLMARK IFN-ALPHA SENSITIVITY ----\n\n"
)

hallmark_result <- readRDS(
  input_paths[["hallmark_result"]]
)

if (!"result" %in% methods::slotNames(hallmark_result)) {
  stop(
    "hallmark_ifn_sensitivity_result.rds is missing the expected result slot."
  )
}

hallmark_result_table <- methods::slot(
  hallmark_result,
  "result"
)

assert_equal(
  nrow(hallmark_result_table),
  136460L,
  "Hallmark LINCS2 signature count"
)

hallmark_table <- readr::read_csv(
  input_paths[["hallmark_table"]],
  show_col_types = FALSE
)

required_hallmark_columns <- c(
  "Compound",
  "Role",
  "THP1_Profiling_Status",
  "NCS",
  "FDR",
  "Direction",
  "Interpretation"
)

missing_hallmark_columns <- setdiff(
  required_hallmark_columns,
  colnames(hallmark_table)
)

if (length(missing_hallmark_columns) > 0L) {
  stop(
    "Hallmark Supplementary Table 4 export is missing required column(s): ",
    paste(missing_hallmark_columns, collapse = ", ")
  )
}

if (
  nrow(hallmark_table) != 11L ||
  anyDuplicated(hallmark_table$Compound)
) {
  stop(
    "Hallmark Supplementary Table 4 must contain exactly 11 unique compounds."
  )
}

expected_hallmark_profiled <- tibble::tribble(
  ~Compound,             ~BRD_ID,          ~NCS,     ~FDR,       ~Interpretation,
  "Rosuvastatin",        "BRD-K82941592",  -2.2427,  0.0001043,  "Significant reversal",
  "Roflumilast",         "BRD-K03194791",  -1.2500,  0.0001043,  "Significant reversal",
  "Baricitinib",         "BRD-K53581288",  -0.7772,  0.08445,    "Negative, not significant",
  "Tofacitinib",         "BRD-K36280065",  -1.0656,  0.0002255,  "Significant reversal",
  "Ruxolitinib",         "BRD-K53972329",  -1.2454,  0.0001043,  "Significant reversal",
  "Leflunomide",         "BRD-K78692225",   1.2074,  0.0001043,  "Significant positive connectivity",
  "Hydroxychloroquine",  "BRD-A99117172",   0.9130,  0.002171,   "Significant positive connectivity"
)

for (i in seq_len(nrow(expected_hallmark_profiled))) {

  expected_row <- expected_hallmark_profiled[i, ]

  exported_row <- hallmark_table |>
    dplyr::filter(
      Compound == expected_row$Compound
    )

  if (nrow(exported_row) != 1L) {
    stop(
      "Expected exactly one Supplementary Table 4 row for ",
      expected_row$Compound,
      "."
    )
  }

  if (exported_row$THP1_Profiling_Status != "Profiled") {
    stop(
      expected_row$Compound,
      " should be marked as Profiled in Supplementary Table 4."
    )
  }

  assert_close(
    exported_row$NCS,
    expected_row$NCS,
    0.002,
    paste0(
      "Table 4 ",
      expected_row$Compound,
      " NCS"
    )
  )

  fdr_tolerance <- if (expected_row$FDR > 0.01) {
    0.001
  } else {
    0.00001
  }

  assert_close(
    exported_row$FDR,
    expected_row$FDR,
    fdr_tolerance,
    paste0(
      "Table 4 ",
      expected_row$Compound,
      " FDR"
    )
  )

  if (exported_row$Interpretation != expected_row$Interpretation) {
    stop(
      "Table 4 ",
      expected_row$Compound,
      " interpretation mismatch."
    )
  }

  raw_row <- hallmark_result_table |>
    dplyr::filter(
      pert == expected_row$BRD_ID,
      toupper(cell) == "THP1"
    )

  if (nrow(raw_row) != 1L) {
    stop(
      "Expected exactly one raw Hallmark THP1 row for ",
      expected_row$Compound,
      " (",
      expected_row$BRD_ID,
      ")."
    )
  }

  assert_close(
    raw_row$NCS,
    exported_row$NCS,
    0.00011,
    paste0(
      "Raw/exported Hallmark NCS for ",
      expected_row$Compound
    )
  )

  assert_close(
    raw_row$WTCS_FDR,
    exported_row$FDR,
    1e-12,
    paste0(
      "Raw/exported Hallmark FDR for ",
      expected_row$Compound
    )
  )
}

NOT_IN_LIBRARY <- "Not in LINCS2 screening library (checked compoundinfo_beta.txt)"
NOT_PROFILED <- "Profiled in LINCS2, no THP1 data"

expected_unprofiled <- c(
  "Upadacitinib" = NOT_IN_LIBRARY,
  "Filgotinib" = NOT_PROFILED,
  "Methotrexate" = NOT_PROFILED,
  "Sulfasalazine" = NOT_PROFILED
)

for (compound_name in names(expected_unprofiled)) {

  exported_row <- hallmark_table |>
    dplyr::filter(
      Compound == compound_name
    )

  if (nrow(exported_row) != 1L) {
    stop(
      "Expected exactly one Supplementary Table 4 row for ",
      compound_name,
      "."
    )
  }

  if (
    exported_row$THP1_Profiling_Status !=
      expected_unprofiled[[compound_name]]
  ) {
    stop(
      "Table 4 ",
      compound_name,
      " profiling-status mismatch."
    )
  }

  if (
    !is.na(exported_row$NCS) ||
    !is.na(exported_row$FDR)
  ) {
    stop(
      "Table 4 ",
      compound_name,
      " should not contain THP1 NCS/FDR values."
    )
  }
}

cat("Verified current Supplementary Table 4 values:\n")
print(
  hallmark_table,
  n = Inf
)
cat(
  "\nQuery definition verified upstream by Script 07: ",
  "97 Hallmark genes; 42 overlap with iNPH UP, 1 with iNPH DOWN, ",
  "43 unique overlaps excluded, 54 non-overlapping genes retained, ",
  "49 represented in LINCS2.\n\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 3 . Tables 1 and 3: independent LINCS2 THP1 ranking
# -----------------------------------------------------------------------------

cat(
  "---- TABLES 1 / 3: LINCS2 THP1 RANKING ----\n\n"
)

lincs_result <- readRDS(
  input_paths[["lincs_result"]]
) |>
  signatureSearch::result()

required_lincs_columns <- c(
  "pert",
  "type",
  "cell",
  "NCS",
  "WTCS",
  "WTCS_FDR"
)

missing_lincs_columns <- setdiff(
  required_lincs_columns,
  colnames(lincs_result)
)

if (length(missing_lincs_columns) > 0L) {
  stop(
    "LINCS2 result is missing required column(s): ",
    paste(missing_lincs_columns, collapse = ", ")
  )
}

assert_equal(
  nrow(lincs_result),
  136460L,
  "Primary LINCS2 signature count"
)

thp1_raw <- lincs_result |>
  dplyr::filter(
    type == "trt_cp",
    toupper(cell) == "THP1",
    grepl("^BRD-", pert),
    !is.na(NCS)
  )

thp1_counts <- thp1_raw |>
  dplyr::count(
    pert,
    name = "n_thp1_rows"
  )

duplicate_thp1 <- thp1_counts |>
  dplyr::filter(
    n_thp1_rows != 1L
  )

if (nrow(duplicate_thp1) > 0L) {
  stop(
    "Expected exactly one THP1 row per BRD compound with THP1 data; ",
    "duplicates found for: ",
    paste(
      duplicate_thp1$pert,
      collapse = ", "
    )
  )
}

ranked_raw <- thp1_raw |>
  dplyr::filter(
    NCS != 0
  ) |>
  dplyr::arrange(
    NCS
  ) |>
  dplyr::transmute(
    rank = dplyr::row_number(),
    pert,
    thp1_NCS = NCS,
    thp1_WTCS = WTCS,
    thp1_FDR = WTCS_FDR
  )

if (
  nrow(ranked_raw) != 597L ||
  anyDuplicated(ranked_raw$pert) ||
  !all(ranked_raw$rank == seq_len(597L))
) {
  stop(
    "Independent raw LINCS2 THP1 ranking failed the finalized 597-compound validation."
  )
}


# -----------------------------------------------------------------------------
# Section 4 . Verify Script 08 full-ranking export against raw LINCS2 results
# -----------------------------------------------------------------------------

cat(
  "---- FULL 597-COMPOUND EXPORT CONSISTENCY ----\n\n"
)

full_export <- readr::read_csv(
  input_paths[["full_thp1_export"]],
  show_col_types = FALSE
)

required_export_columns <- c(
  "rank",
  "pert",
  "pert_iname",
  "n_cell_lines",
  "thp1_NCS",
  "thp1_WTCS",
  "thp1_FDR"
)

missing_export_columns <- setdiff(
  required_export_columns,
  colnames(full_export)
)

if (length(missing_export_columns) > 0L) {
  stop(
    "Script 08 export is missing required column(s): ",
    paste(missing_export_columns, collapse = ", ")
  )
}

if (
  nrow(full_export) != 597L ||
  anyDuplicated(full_export$pert) ||
  anyNA(full_export$rank) ||
  !all(full_export$rank == seq_len(597L))
) {
  stop(
    "Script 08 full THP1 export failed its 597-row/rank validation."
  )
}

raw_export_check <- full_export |>
  dplyr::select(
    rank,
    pert,
    export_NCS = thp1_NCS,
    export_WTCS = thp1_WTCS,
    export_FDR = thp1_FDR
  ) |>
  dplyr::left_join(
    ranked_raw |>
      dplyr::rename(
        raw_rank = rank,
        raw_NCS = thp1_NCS,
        raw_WTCS = thp1_WTCS,
        raw_FDR = thp1_FDR
      ),
    by = "pert"
  )

if (anyNA(raw_export_check$raw_rank)) {
  stop(
    "One or more Script 08 compounds are absent from the independently rebuilt raw ranking."
  )
}

assert_true(
  all(
    raw_export_check$rank ==
      raw_export_check$raw_rank
  ),
  "Script 08 ranks match the independently rebuilt raw ranking"
)

assert_true(
  isTRUE(
    all.equal(
      raw_export_check$export_NCS,
      raw_export_check$raw_NCS,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  ),
  "Script 08 NCS values match raw LINCS2"
)

assert_true(
  isTRUE(
    all.equal(
      raw_export_check$export_WTCS,
      raw_export_check$raw_WTCS,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  ),
  "Script 08 WTCS values match raw LINCS2"
)

assert_true(
  isTRUE(
    all.equal(
      raw_export_check$export_FDR,
      raw_export_check$raw_FDR,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  ),
  "Script 08 FDR values match raw LINCS2"
)

cat(
  "All 597 Script 08 rows match the independently rebuilt raw THP1 ranking.\n\n"
)


# -----------------------------------------------------------------------------
# Section 5 . Top 30 plus baricitinib rank verification
# -----------------------------------------------------------------------------

check_ids <- tibble::tribble(
  ~expected_rank, ~pert,           ~display_name,
  1L,  "BRD-K82941592", "Rosuvastatin",
  2L,  "BRD-K43002773", "GDC-0068",
  3L,  "BRD-K12683773", "Methyl-fasudil",
  4L,  "BRD-K92317137", "BRD-K92317137",
  5L,  "BRD-K92303087", "Oglemilast",
  6L,  "BRD-A02367930", "Ethinyl estradiol",
  7L,  "BRD-K93255255", "Indoximod",
  8L,  "BRD-K75641298", "Metoclopramide",
  9L,  "BRD-K59013864", "AM-24",
  10L, "BRD-K40624912", "ZM-39923",
  11L, "BRD-K84925803", "T-5224",
  12L, "BRD-K07859598", "GW-5074",
  13L, "BRD-K14993104", "Bemegride",
  14L, "BRD-K03194791", "Roflumilast",
  15L, "BRD-K01436366", "XMD-1150",
  16L, "BRD-K96615647", "AV-608",
  17L, "BRD-A13084692", "Troglitazone",
  18L, "BRD-K00610438", "Altanserin",
  19L, "BRD-K60038276", "Irbesartan",
  20L, "BRD-A14262390", "Guanadrel",
  21L, "BRD-K77987382", "Mebendazole",
  22L, "BRD-K82181189", "Ischemin",
  23L, "BRD-K59962020", "CHEMBL-374350",
  24L, "BRD-K94830329", "Ataluren",
  25L, "BRD-K59019422", "BRD-K59019422",
  26L, "BRD-K29950728", "Clomifene",
  27L, "BRD-K70642949", "GSK-2334470",
  28L, "BRD-K01674964", "Histamine",
  29L, "BRD-K99749624", "Linifanib",
  30L, "BRD-K68548958", "C-646",
  39L, "BRD-K53581288", "Baricitinib"
)

verify_df <- check_ids |>
  dplyr::left_join(
    ranked_raw,
    by = "pert"
  )

if (
  nrow(verify_df) != 31L ||
  anyNA(verify_df$rank)
) {
  stop(
    "Top-30 plus baricitinib verification set is incomplete."
  )
}

verify_df <- verify_df |>
  dplyr::mutate(
    rank_match =
      rank == expected_rank
  )

n_rank_mismatch <- sum(
  !verify_df$rank_match
)

if (n_rank_mismatch != 0L) {
  stop(
    n_rank_mismatch,
    " rank mismatch(es) detected in the top-30 plus baricitinib verification set."
  )
}

cat(
  sprintf(
    "%-8s %-20s %-16s %6s %9s %9s %11s %8s\n",
    "Expected",
    "Compound",
    "BRD_ID",
    "Rank",
    "NCS",
    "WTCS",
    "FDR",
    "Match"
  )
)

for (i in seq_len(nrow(verify_df))) {

  r <- verify_df[i, ]

  cat(
    sprintf(
      "%-8d %-20s %-16s %6d %9.4f %9.4f %11.3g %8s\n",
      r$expected_rank,
      r$display_name,
      r$pert,
      r$rank,
      r$thp1_NCS,
      r$thp1_WTCS,
      r$thp1_FDR,
      ifelse(
        r$rank_match,
        "OK",
        "CHECK"
      )
    )
  )
}

cat(
  "\n0 of 31 rank mismatches.\n\n"
)


# -----------------------------------------------------------------------------
# Section 6 . Table 1: eight reviewed FDA-approved reversers
# -----------------------------------------------------------------------------

cat(
  "---- TABLE 1: EIGHT REVIEWED FDA-APPROVED REVERSERS ----\n\n"
)

final_candidates <- readRDS(
  input_paths[["final_candidates"]]
)

significant_reversers <- readRDS(
  input_paths[["significant_reversers"]]
)

expected_table1 <- tibble::tribble(
  ~rank, ~pert,           ~compound,             ~NCS,    ~WTCS,   ~FDR,
  1L,   "BRD-K82941592", "rosuvastatin",        -1.845,  -0.487,  0.000614,
  6L,   "BRD-A02367930", "ethinyl estradiol",   -1.477,  -0.390,  0.000614,
  8L,   "BRD-K75641298", "metoclopramide",      -1.442,  -0.380,  0.000614,
  14L,  "BRD-K03194791", "roflumilast",         -1.375,  -0.363,  0.000614,
  19L,  "BRD-K60038276", "irbesartan",          -1.312,  -0.346,  0.000614,
  21L,  "BRD-K77987382", "mebendazole",         -1.299,  -0.343,  0.000614,
  26L,  "BRD-K29950728", "clomifene",           -1.282,  -0.338,  0.000618,
  39L,  "BRD-K53581288", "baricitinib",         -1.231,  -0.325,  0.000817
)

if (
  nrow(final_candidates) != 8L ||
  !setequal(
    final_candidates$pert,
    expected_table1$pert
  )
) {
  stop(
    "Final candidate object does not contain the finalized eight-compound Table 1 set."
  )
}

table1_raw <- expected_table1 |>
  dplyr::left_join(
    ranked_raw,
    by = "pert",
    suffix = c("_expected", "_raw")
  )

if (
  nrow(table1_raw) != 8L ||
  anyNA(table1_raw$rank_raw)
) {
  stop(
    "One or more Table 1 compounds are missing from the raw THP1 ranking."
  )
}

for (i in seq_len(nrow(table1_raw))) {

  r <- table1_raw[i, ]

  assert_equal(
    r$rank_raw,
    r$rank_expected,
    paste0(
      "Table 1 ",
      r$compound,
      " rank"
    )
  )

  assert_close(
    r$thp1_NCS,
    r$NCS,
    0.002,
    paste0(
      "Table 1 ",
      r$compound,
      " NCS"
    )
  )

  assert_close(
    r$thp1_WTCS,
    r$WTCS,
    0.002,
    paste0(
      "Table 1 ",
      r$compound,
      " WTCS"
    )
  )

  assert_close(
    r$thp1_FDR,
    r$FDR,
    0.00001,
    paste0(
      "Table 1 ",
      r$compound,
      " FDR"
    )
  )
}

assert_equal(
  nrow(significant_reversers),
  138L,
  "Significant strong THP1 reverser count"
)

if (!"approval_status" %in% colnames(significant_reversers)) {
  stop(
    "Significant-reverser object is missing approval_status."
  )
}

n_unreviewed <- sum(
  significant_reversers$approval_status ==
    "UNVERIFIED - outside reviewed set"
)

assert_equal(
  n_unreviewed,
  107L,
  "Significant reversers outside the reviewed regulatory-status set"
)

cat("Verified Table 1 computational values:\n")
print(
  table1_raw |>
    dplyr::transmute(
      rank = rank_raw,
      pert,
      compound,
      NCS = thp1_NCS,
      WTCS = thp1_WTCS,
      FDR = thp1_FDR
    )
)

cat(
  "\nSignificant strong THP1 reversers: ",
  nrow(significant_reversers),
  " | outside reviewed regulatory-status set: ",
  n_unreviewed,
  "\n\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 7 . External-annotation reminder
# -----------------------------------------------------------------------------

cat(
  "---- EXTERNAL ANNOTATION REVIEW ----\n\n"
)

cat(
  "Regulatory status and pharmacologic class are not established by the ",
  "LINCS2 connectivity objects. The saved candidate annotations remain ",
  "restricted to the prespecified top-30 review plus post hoc baricitinib ",
  "and should remain traceable to the cited regulatory or clinical sources.\n\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 8 . Final verification status
# -----------------------------------------------------------------------------

cat(
  "=====================================================================\n"
)
cat(
  "  MANUSCRIPT TABLE VERIFICATION PASSED\n"
)
cat(
  "=====================================================================\n"
)
cat(
  "  Supplementary Table 2 contextualization values : PASS\n"
)
cat(
  "  Supplementary Table 4 Hallmark values           : PASS\n"
)
cat(
  "  Raw LINCS2 non-zero THP1 ranking                : 597 compounds\n"
)
cat(
  "  Script 08 full-ranking export                   : 597/597 consistent\n"
)
cat(
  "  Top-30 plus baricitinib rank mismatches         : 0/31\n"
)
cat(
  "  Reviewed FDA-approved reversers                 : 8\n"
)
cat(
  "  Significant strong reversers                    : 138\n"
)
cat(
  "  Significant reversers outside reviewed set      : 107\n"
)
cat(
  "  Priority ranks                                  : 1 / 14 / 39\n"
)
cat(
  "=====================================================================\n"
)
