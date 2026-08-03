# =============================================================================
# Generate the full 597-compound repository CSV
# =============================================================================
# The Table 1 footnote and Supplementary Table 3 caption both now state the
# full list of 597 compounds with non-zero THP1 connectivity is "available
# as a repository CSV" (per Dr. Nila's correction, 2026-08-XX). That file did
# not previously exist - phase3_top200_reversal_candidates.csv only has 200
# rows, and phase3_step2_significant_candidates.csv only has the 138 that
# passed both significance criteria. This closes that gap.
#
# SCOPE: connectivity metrics only (NCS, WTCS, FDR, rank, cell-line count,
# compound name where identifiable). Regulatory status and pharmacologic
# class are NOT included here - those were only manually verified for the
# top 30 + baricitinib (Table 3), and extending that column to all 597 would
# imply a review that was never done. This CSV makes no claim beyond what
# 03_lincs_connectivity.R actually computed.
#
# Reads the already-saved LINCS2 result object - does NOT re-run gess_lincs().
# =============================================================================

library(here); library(dplyr); library(readr)
PROCESSED_DIR <- here::here("data", "processed")

res <- readRDS(file.path(PROCESSED_DIR, "lincs2_gess_result.rds")) |> signatureSearch::result()

compoundinfo_beta <- readr::read_tsv(
  here::here("data", "raw", "compoundinfo_beta.txt"), show_col_types = FALSE)
compoundinfo_gse <- readr::read_tsv(
  here::here("data", "raw", "GSE70138_Broad_LINCS_pert_info.txt.gz"), show_col_types = FALSE)

thp1_full <- res |>
  dplyr::filter(type == "trt_cp", toupper(cell) == "THP1", grepl("^BRD-", pert), !is.na(NCS), NCS != 0) |>
  dplyr::select(pert, thp1_NCS = NCS, thp1_WTCS = WTCS, thp1_FDR = WTCS_FDR)

stopifnot(nrow(thp1_full) == 597)  # matches every prior verified run - stop if this ever changes

n_cell_lines <- res |>
  dplyr::filter(type == "trt_cp", grepl("^BRD-", pert), !is.na(NCS)) |>
  dplyr::group_by(pert) |>
  dplyr::summarise(n_cell_lines = dplyr::n_distinct(cell), .groups = "drop")

# Both reference files carry duplicate rows for at least some pert_ids
# (confirmed earlier this project - e.g. hydroxychloroquine had two identical
# BRD-A99117172 rows in compoundinfo_beta). An un-deduplicated left_join fans
# out the left side on any such duplicate, which is exactly what happened
# here (597 -> 1,321 rows, GDC-0068 tripled). Deduplicate on pert_id first.
compoundinfo_gse_dedup <- compoundinfo_gse |>
  dplyr::select(pert_id, pert_iname) |>
  dplyr::distinct(pert_id, .keep_all = TRUE)
compoundinfo_beta_dedup <- compoundinfo_beta |>
  dplyr::select(pert_id, cmap_name) |>
  dplyr::distinct(pert_id, .keep_all = TRUE)

named <- thp1_full |>
  dplyr::left_join(n_cell_lines, by = "pert") |>
  dplyr::left_join(compoundinfo_gse_dedup, by = c("pert" = "pert_id")) |>
  dplyr::left_join(compoundinfo_beta_dedup, by = c("pert" = "pert_id")) |>
  dplyr::mutate(
    pert_iname = dplyr::case_when(
      !is.na(pert_iname) & pert_iname != pert ~ pert_iname,
      !is.na(cmap_name) & cmap_name != pert    ~ cmap_name,
      TRUE                                      ~ "UNIDENTIFIED"
    )
  ) |>
  dplyr::select(-cmap_name) |>
  dplyr::arrange(thp1_NCS) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::select(rank, pert, pert_iname, n_cell_lines, thp1_NCS, thp1_WTCS, thp1_FDR)

# Hard guard on the FINAL output, not just the intermediate step - this is
# what should have caught the bug immediately instead of needing a manual check.
stopifnot(nrow(named) == 597)
stopifnot(!anyDuplicated(named$pert))

cat("Full THP1 non-zero-connectivity compound list: ", nrow(named), " rows (expect 597)\n", sep = "")
cat("Identified by name: ", sum(named$pert_iname != "UNIDENTIFIED"), " of ", nrow(named), "\n", sep = "")
cat("Rank 1-5 check:\n")
print(head(named, 5))

out_path <- file.path(PROCESSED_DIR, "supp_table3_full_597_compounds.csv")
readr::write_csv(named, out_path)
cat("\nSaved:", out_path, "\n")
cat("NOTE: connectivity metrics only - no regulatory status or pharmacologic\n")
cat("class columns. Those exist only for the top 30 + baricitinib (see the\n")
cat("main Supplementary Table 3), not for all 597 rows here.\n")
