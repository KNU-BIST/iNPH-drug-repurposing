# =============================================================================
# Export the complete non-zero THP1 compound ranking
# Script: 08_export_full_thp1_compound_ranking.R
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
# Export the complete 597-compound THP1 ranking used by the finalized LINCS2
# analysis.
#
# Scope
# -----
# Connectivity metrics and compound identifiers only. Regulatory status and
# pharmacologic class were reviewed only for the prespecified top 30 compounds
# plus post hoc baricitinib and are intentionally not extrapolated to all 597
# compounds.
#
# Input
# -----
# data/processed/lincs2_compound_summary.rds
#
# Output
# ------
# data/processed/lincs2_thp1_nonzero_597_compounds.csv
#
# Reproducibility notes
# ---------------------
# - Compound-name resolution is performed upstream in
#   03_run_lincs2_connectivity.R using both LINCS metadata sources.
# - This script does not re-read or re-join raw compound metadata, avoiding a
#   second independent name-resolution path that could drift from Script 03.
# - Ranking is based on full-precision THP1 NCS in ascending order.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0 . Environment and input validation
# -----------------------------------------------------------------------------

required_packages <- c(
  "here",
  "dplyr",
  "readr"
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
})

PROCESSED_DIR <- here::here(
  "data",
  "processed"
)

COMPOUND_SUMMARY_PATH <- file.path(
  PROCESSED_DIR,
  "lincs2_compound_summary.rds"
)

if (
  !file.exists(COMPOUND_SUMMARY_PATH) ||
  is.na(file.info(COMPOUND_SUMMARY_PATH)$size) ||
  file.info(COMPOUND_SUMMARY_PATH)$size <= 0
) {
  stop(
    "Required Script 03 output is missing or empty: ",
    COMPOUND_SUMMARY_PATH
  )
}

cat(
  "Compound-summary MD5: ",
  unname(
    tools::md5sum(
      COMPOUND_SUMMARY_PATH
    )
  ),
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 1 . Load and validate the finalized Script 03 compound summary
# -----------------------------------------------------------------------------

compound_summary <- readRDS(
  COMPOUND_SUMMARY_PATH
)

required_columns <- c(
  "pert",
  "pert_iname",
  "n_cell_lines",
  "thp1_NCS",
  "thp1_WTCS",
  "thp1_FDR"
)

missing_columns <- setdiff(
  required_columns,
  colnames(compound_summary)
)

if (length(missing_columns) > 0L) {
  stop(
    "lincs2_compound_summary.rds is missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

if (
  nrow(compound_summary) != 30411L ||
  anyDuplicated(compound_summary$pert)
) {
  stop(
    "Script 03 compound summary does not match the finalized analysis: ",
    "expected 30,411 unique BRD compounds."
  )
}

if (
  any(
    is.na(
      compound_summary$n_cell_lines
    )
  ) ||
  any(
    compound_summary$n_cell_lines < 1L
  )
) {
  stop(
    "Invalid cell-line counts detected in the Script 03 compound summary."
  )
}

cat(
  "Validated Script 03 compound summary: ",
  nrow(compound_summary),
  " unique BRD compounds.\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# Section 2 . Build the complete non-zero THP1 ranking
# -----------------------------------------------------------------------------

ranked_thp1 <- compound_summary |>
  dplyr::filter(
    !is.na(thp1_NCS),
    thp1_NCS != 0
  ) |>
  dplyr::arrange(
    thp1_NCS
  ) |>
  dplyr::mutate(
    rank = dplyr::row_number()
  ) |>
  dplyr::select(
    rank,
    pert,
    pert_iname,
    n_cell_lines,
    thp1_NCS,
    thp1_WTCS,
    thp1_FDR
  )

if (
  nrow(ranked_thp1) != 597L ||
  anyDuplicated(ranked_thp1$pert) ||
  !identical(
    ranked_thp1$rank,
    seq_len(597L)
  ) ||
  anyNA(ranked_thp1$thp1_NCS) ||
  any(ranked_thp1$thp1_NCS == 0)
) {
  stop(
    "Complete THP1 ranking failed validation: expected exactly 597 unique ",
    "compounds with sequential ranks and non-zero, non-missing THP1 NCS."
  )
}

if (
  anyNA(ranked_thp1$thp1_WTCS) ||
  anyNA(ranked_thp1$thp1_FDR)
) {
  stop(
    "Missing THP1 WTCS or WTCS FDR values detected among the 597 ranked compounds."
  )
}


# -----------------------------------------------------------------------------
# Section 3 . Critical rank checks
# -----------------------------------------------------------------------------

critical_expected <- c(
  "BRD-K82941592" = 1L,
  "BRD-K03194791" = 14L,
  "BRD-K53581288" = 39L
)

critical_rows <- ranked_thp1 |>
  dplyr::filter(
    pert %in%
      names(
        critical_expected
      )
  )

if (nrow(critical_rows) != 3L) {
  stop(
    "One or more manuscript-priority compounds are missing from the full THP1 ranking."
  )
}

critical_observed <- setNames(
  critical_rows$rank,
  critical_rows$pert
)

if (
  !all(
    critical_observed[
      names(critical_expected)
    ] ==
      critical_expected
  )
) {
  stop(
    "Priority-compound ranks differ from the finalized analysis."
  )
}

cat(
  "\nCritical THP1 rank check:\n"
)

print(
  critical_rows |>
    dplyr::arrange(
      rank
    )
)


# -----------------------------------------------------------------------------
# Section 4 . Descriptive name-coverage check
# -----------------------------------------------------------------------------

n_identified <- sum(
  !is.na(ranked_thp1$pert_iname) &
    ranked_thp1$pert_iname != "" &
    ranked_thp1$pert_iname != "UNIDENTIFIED"
)

n_unidentified <- nrow(ranked_thp1) -
  n_identified

cat(
  "\nFull THP1 non-zero-connectivity ranking: ",
  nrow(ranked_thp1),
  " compounds.\n",
  sep = ""
)

cat(
  "Named compounds: ",
  n_identified,
  " | Unidentified: ",
  n_unidentified,
  "\n",
  sep = ""
)

cat(
  "\nTop five compounds by full-precision THP1 NCS:\n"
)

print(
  ranked_thp1 |>
    dplyr::slice_head(
      n = 5L
    )
)


# -----------------------------------------------------------------------------
# Section 5 . Export and verify
# -----------------------------------------------------------------------------

OUTPUT_PATH <- file.path(
  PROCESSED_DIR,
  "lincs2_thp1_nonzero_597_compounds.csv"
)

readr::write_csv(
  ranked_thp1,
  OUTPUT_PATH
)

if (
  !file.exists(OUTPUT_PATH) ||
  is.na(file.info(OUTPUT_PATH)$size) ||
  file.info(OUTPUT_PATH)$size <= 0
) {
  stop(
    "Failed to create the complete THP1 ranking CSV: ",
    OUTPUT_PATH
  )
}

export_check <- readr::read_csv(
  OUTPUT_PATH,
  show_col_types = FALSE
)

# readr may infer an integer-valued CSV column such as rank as numeric/double.
# Therefore validate rank values rather than requiring identical R storage mode.
rank_values_ok <- (
  length(export_check$rank) == 597L &&
    !anyNA(export_check$rank) &&
    all(export_check$rank == seq_len(597L)) &&
    all(export_check$rank == floor(export_check$rank))
)

compound_order_ok <- identical(
  as.character(export_check$pert),
  as.character(ranked_thp1$pert)
)

ncs_values_ok <- isTRUE(
  all.equal(
    export_check$thp1_NCS,
    ranked_thp1$thp1_NCS,
    tolerance = 1e-12,
    check.attributes = FALSE
  )
)

wtcs_values_ok <- isTRUE(
  all.equal(
    export_check$thp1_WTCS,
    ranked_thp1$thp1_WTCS,
    tolerance = 1e-12,
    check.attributes = FALSE
  )
)

fdr_values_ok <- isTRUE(
  all.equal(
    export_check$thp1_FDR,
    ranked_thp1$thp1_FDR,
    tolerance = 1e-12,
    check.attributes = FALSE
  )
)

if (
  nrow(export_check) != 597L ||
  anyDuplicated(export_check$pert) ||
  !rank_values_ok ||
  !compound_order_ok ||
  !ncs_values_ok ||
  !wtcs_values_ok ||
  !fdr_values_ok
) {
  stop(
    "Post-write verification of the complete THP1 ranking CSV failed."
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "  COMPLETE THP1 COMPOUND RANKING EXPORT FINISHED\n"
)

cat(
  "============================================================\n"
)

cat(
  "  Ranked non-zero THP1 compounds : ",
  nrow(ranked_thp1),
  "\n",
  sep = ""
)

cat(
  "  Named compounds                : ",
  n_identified,
  "\n",
  sep = ""
)

cat(
  "  Unidentified compounds         : ",
  n_unidentified,
  "\n",
  sep = ""
)

cat(
  "  Priority ranks                 : rosuvastatin 1 | roflumilast 14 | ",
  "baricitinib 39\n",
  sep = ""
)

cat(
  "============================================================\n"
)

cat(
  "\nSaved and verified: ",
  OUTPUT_PATH,
  "\n",
  sep = ""
)

cat(
  "NOTE: this 597-row export contains connectivity metrics and compound ",
  "identifiers only. Regulatory-status and pharmacologic-class annotations ",
  "are restricted to the top-30 review plus post hoc baricitinib in the ",
  "finalized Script 03 outputs.\n",
  sep = ""
)
