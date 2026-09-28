# Verify integrated ranks and candidate donor-omission results
# from the archived Supplementary Data 2 workbook.
# Run from the repository root.

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the readxl package before running this check.")
}

workbook <- file.path(
  "reproducibility_inputs",
  "Supplementary_Data_2_Candidate_Robustness.xlsx"
)
stopifnot(file.exists(workbook))

all_results <- as.data.frame(
  readxl::read_excel(workbook, sheet = "All query results")
)
saved_summary <- as.data.frame(
  readxl::read_excel(workbook, sheet = "All-compound summary")
)
saved_loo <- as.data.frame(
  readxl::read_excel(workbook, sheet = "Candidate donor omissions")
)

required_results <- c(
  "query_id", "query_role", "pert", "pert_iname",
  "NCS", "WTCS_FDR", "rank_all"
)
required_summary <- c(
  "pert", "pert_iname", "objective_robustness_rank",
  "core_query_n", "operational_hit_n",
  "loo_query_n", "loo_operational_hit_n"
)
required_loo <- c(
  "pert", "pert_iname", "excluded_patient", "query_id",
  "NCS", "WTCS_FDR", "rank_all"
)

stopifnot(
  all(required_results %in% names(all_results)),
  all(required_summary %in% names(saved_summary)),
  all(required_loo %in% names(saved_loo))
)

# The primary query plus the 10 biological sensitivity queries.
core <- all_results[
  all_results$query_role %in% c("primary", "biological_sensitivity"),
]
core$query_id <- as.character(core$query_id)
core$pert <- as.character(core$pert)

core_query_ids <- unique(core$query_id)
if (length(core_query_ids) != 11L) {
  stop("Expected 11 core queries; found ", length(core_query_ids), ".")
}
if (!all(table(core$query_id) == 1376L)) {
  stop("Each core query must contain exactly 1,376 compounds.")
}
if (anyDuplicated(paste(core$query_id, core$pert, sep = "::"))) {
  stop("Duplicate compound/query rows found in the core results.")
}
if (anyNA(core$NCS) || anyNA(core$WTCS_FDR)) {
  stop("Missing NCS or FDR values found in the core results.")
}
if (!"PRIMARY_150" %in% core$query_id) {
  stop("PRIMARY_150 is missing from the core queries.")
}

# Recalculate average-tie ranks and percentiles within each query.
core$percentile_recalculated <- ave(
  core$NCS,
  core$query_id,
  FUN = function(x) {
    (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  }
)

# Rebuild the integrated summary from the 11 core queries.
median_percentile <- aggregate(
  percentile_recalculated ~ pert, data = core, FUN = median
)
names(median_percentile)[2] <- "median_core_percentile"

negative_fraction <- aggregate(
  I(NCS < 0) ~ pert, data = core, FUN = mean
)
names(negative_fraction)[2] <- "inverse_proportion"

primary <- core[
  core$query_id == "PRIMARY_150",
  c("pert", "percentile_recalculated")
]
names(primary)[2] <- "primary_percentile"

core$operational_hit <- as.integer(
  core$NCS < -0.5 & core$WTCS_FDR < 0.05
)
hit_count <- aggregate(operational_hit ~ pert, data = core, FUN = sum)

query_count <- aggregate(
  query_id ~ pert, data = core,
  FUN = function(x) length(unique(x))
)
names(query_count)[2] <- "core_query_n"

recalculated <- Reduce(
  function(x, y) merge(x, y, by = "pert", all = TRUE),
  list(median_percentile, negative_fraction, primary, hit_count, query_count)
)

# Match the finalizer's deterministic ordering.
ordering <- order(
  recalculated$median_core_percentile,
  -recalculated$inverse_proportion,
  recalculated$primary_percentile,
  recalculated$pert,
  method = "radix"
)
recalculated$rank_recalculated <- integer(nrow(recalculated))
recalculated$rank_recalculated[ordering] <- seq_len(nrow(recalculated))

comparison <- merge(
  recalculated,
  saved_summary[
    , c("pert", "objective_robustness_rank",
        "core_query_n", "operational_hit_n")
  ],
  by = "pert",
  all = TRUE
)

if (nrow(comparison) != 1376L ||
    anyNA(comparison$objective_robustness_rank) ||
    any(comparison$rank_recalculated !=
        comparison$objective_robustness_rank) ||
    any(comparison$core_query_n.x != comparison$core_query_n.y) ||
    any(comparison$operational_hit != comparison$operational_hit_n)) {
  stop("Integrated ranking or core-query counts do not reproduce.")
}

# Check the highlighted compounds and their nine donor omissions.
targets <- c("rosuvastatin", "roflumilast")
target_summary <- saved_summary[
  tolower(saved_summary$pert_iname) %in% targets,
]
target_summary <- target_summary[
  match(targets, tolower(target_summary$pert_iname)),
]

if (nrow(target_summary) != 2L ||
    any(target_summary$objective_robustness_rank != c(1, 15)) ||
    any(target_summary$operational_hit_n != c(11, 10)) ||
    any(target_summary$loo_query_n != 9) ||
    any(target_summary$loo_operational_hit_n != 9)) {
  stop("Highlighted-compound summary does not match the reported values.")
}

candidate_loo <- saved_loo[
  tolower(saved_loo$pert_iname) %in% targets,
]
if (nrow(candidate_loo) != 18L ||
    anyDuplicated(paste(candidate_loo$pert, candidate_loo$query_id, sep = "::"))) {
  stop("Expected 18 unique candidate donor-omission rows.")
}

loo_counts <- table(tolower(candidate_loo$pert_iname))
if (!all(loo_counts[targets] == 9L)) {
  stop("Each highlighted compound must have nine donor-omission rows.")
}

candidate_loo$operational_hit <- as.integer(
  candidate_loo$NCS < -0.5 & candidate_loo$WTCS_FDR < 0.05
)
loo_hits <- tapply(
  candidate_loo$operational_hit,
  tolower(candidate_loo$pert_iname),
  sum
)
if (!all(loo_hits[targets] == 9L)) {
  stop("Donor-omission operational-hit counts do not reproduce.")
}

# Confirm the donor-omission sheet agrees with its rows in All query results.
all_loo <- all_results[
  grepl("^LOO_DROP_", all_results$query_id) &
    tolower(all_results$pert_iname) %in% targets,
]
key_saved <- paste(candidate_loo$pert, candidate_loo$query_id, sep = "::")
key_all <- paste(all_loo$pert, all_loo$query_id, sep = "::")
match_rows <- match(key_saved, key_all)

if (nrow(all_loo) != 18L || anyNA(match_rows)) {
  stop("Donor-omission rows do not match All query results.")
}
for (column in c("NCS", "WTCS_FDR", "rank_all")) {
  if (!isTRUE(all.equal(
    as.numeric(candidate_loo[[column]]),
    as.numeric(all_loo[[column]][match_rows]),
    tolerance = 1e-12
  ))) {
    stop("Donor-omission mismatch in column: ", column)
  }
}

cat("PASS: all 1,376 integrated ranks reproduced from saved query tables.\n")
cat("PASS: rosuvastatin rank 1 and 11/11 core operational hits.\n")
cat("PASS: roflumilast rank 15 and 10/11 core operational hits.\n")
cat("PASS: nine donor omissions per candidate; all 18 rows cross-checked.\n")
cat("No LINCS connectivity queries were rerun.\n")
