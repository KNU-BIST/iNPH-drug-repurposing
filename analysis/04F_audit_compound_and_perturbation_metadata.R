# =============================================================================
# Build: 04F_METADATA_AUDIT_V1_20260907
# Install as analysis/04F_audit_compound_and_perturbation_metadata.R
# =============================================================================
# Read-only audit of local compound annotations and saved LINCS result profiles.
# Uses the completed 04E average-tie summary. Does not change names, combine BRD
# IDs, query LINCS, download metadata, or assert current regulatory approval.
# Multiple source names can be legitimate synonyms; flags request adjudication.
# Dose/time cannot be assigned from a compound-by-cell result without an explicit
# mapping to the executed reference profile. Available signature-file headers are
# inventoried, but compound-only joins are not used to invent that linkage.

local({
  BUILD <- "04F_METADATA_AUDIT_V1_20260907"
  LOCK <- paste0("1e3a2b827a1429e406303871d547a9dec", "07ee4068078c666df69ca82e5428380")
  check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  clean <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | tolower(x) %in% c("", "na", "nan", "null", "unknown", "unidentified", "-666")] <- NA_character_
    x
  }
  values <- function(x) sort(unique(stats::na.omit(clean(x))), method = "radix")
  collapse <- function(x) paste(values(x), collapse = " | ")
  column <- function(d, nm) if (nm %in% names(d)) d[[nm]] else character()
  equal_name <- function(x) tolower(trimws(as.character(x)))
  run <- function() {
    start <- Sys.time()
    for (pkg in c("readr", "signatureSearch"))
      check(requireNamespace(pkg, quietly = TRUE), paste("Required package missing:", pkg))
    selected <- Sys.getenv("INPH_AVERAGE_TIE_SUMMARY_RDS", unset = "")
    files <- if (nzchar(selected)) selected else list.files(
      "results/revision/average_tie_finalization",
      pattern = "^revision_candidate_robustness_average_ties\\.rds$", recursive = TRUE, full.names = TRUE)
    files <- files[file.exists(files) & file.exists(file.path(dirname(files), "_SUCCESS.txt"))]
    check(length(files) > 0L, "Run 04E successfully before the metadata audit.")
    files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
    summary_path <- files[1L]
    revised <- readRDS(summary_path)
    check(identical(revised$status, "PASS") && identical(revised$ranking_method, "average") &&
            identical(revised$manifest_content_sha256, LOCK), "Invalid completed average-tie summary.")
    previous_inputs <- as.data.frame(revised$input_provenance)
    check(all(file.exists(previous_inputs$input)) &&
            identical(unname(tools::md5sum(previous_inputs$input)), as.character(previous_inputs$md5)),
          "An input used by 04E has changed; investigate before continuing.")
    cand <- as.data.frame(revised$candidate_summary)
    check(nrow(cand) == 1376L && !anyDuplicated(cand$pert), "Invalid candidate universe.")
    primary_path <- "data/processed/lincs2_gess_result.rds"
    compound_path <- "data/processed/lincs2_compound_summary.rds"
    sources <- c(GSE70138 = "data/raw/GSE70138_Broad_LINCS_pert_info.txt.gz",
                 compoundinfo_beta = "data/raw/compoundinfo_beta.txt")
    check(all(file.exists(c(primary_path, compound_path))), "Saved primary results/compound summary missing.")
    # A missing annotation source is recorded as a gap, not silently ignored.
    source_present <- file.exists(sources)
    required <- unique(c(summary_path, primary_path, compound_path, unname(sources[source_present])))
    own <- "analysis/04F_audit_compound_and_perturbation_metadata.R"
    if (file.exists(own)) required <- c(required, own)
    md5_before <- unname(tools::md5sum(required))
    check(!anyNA(md5_before), "An input checksum could not be obtained.")
    raw <- list()
    source_inventory <- list()
    source_rows <- list()
    for (src in names(sources)) {
      path <- sources[[src]]
      if (!file.exists(path)) {
        raw[[src]] <- data.frame(pert_id = character())
        source_inventory[[src]] <- data.frame(source = src, file = path, exists = FALSE,
          rows_n = 0L, relevant_rows_n = 0L, columns = "", status = "missing_source")
        next
      }
      parsed <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()),
                                progress = FALSE, show_col_types = FALSE, name_repair = "minimal")
      check(nrow(readr::problems(parsed)) == 0L, paste("Metadata parsing problems:", src))
      d <- as.data.frame(parsed)
      check(!anyDuplicated(names(d)) && "pert_id" %in% names(d), paste("Invalid metadata schema:", src))
      d$pert_id <- trimws(as.character(d$pert_id))
      rows <- which(d$pert_id %in% cand$pert)
      source_inventory[[src]] <- data.frame(source = src, file = path, exists = TRUE,
        rows_n = nrow(d), relevant_rows_n = length(rows), columns = paste(names(d), collapse = " | "), status = "read")
      raw[[src]] <- d[rows, , drop = FALSE]
      # Preserve relevant raw rows, including duplicates and all columns.
      source_rows[[src]] <- cbind(source = rep(src, length(rows)), source_row = rows, d[rows, , drop = FALSE])
    }
    compound <- as.data.frame(readRDS(compound_path))
    check(all(c("pert", "pert_iname") %in% names(compound)) && !anyDuplicated(compound$pert),
          "Saved compound summary must have unique identifiers.")
    matched <- match(cand$pert, compound$pert)
    check(!anyNA(matched), "A ranked compound is absent from the saved compound summary.")
    compound <- compound[matched, , drop = FALSE]
    primary_obj <- readRDS(primary_path)
    check(isS4(primary_obj) && "result" %in% methods::slotNames(primary_obj), "Expected a saved gessResult.")
    profiles <- as.data.frame(methods::slot(primary_obj, "result"))
    check(all(c("pert", "cell", "type", "NCS") %in% names(profiles)), "Primary result schema is incomplete.")
    valid <- !is.na(profiles$pert) & !is.na(profiles$cell) & !is.na(profiles$type)
    thp <- profiles[valid & toupper(profiles$cell) == "THP1" & profiles$type == "trt_cp" &
                      profiles$pert %in% cand$pert, , drop = FALSE]
    check(nrow(thp) == 1376L && !anyDuplicated(thp$pert) && setequal(thp$pert, cand$pert),
          "Primary THP1 profiles do not match the ranked compound universe.")
    thp <- thp[match(cand$pert, thp$pert), ]
    check(identical(as.numeric(thp$NCS), as.numeric(cand$primary_NCS)), "Primary scores differ from 04E.")
    dose_fields <- intersect(c("pert_dose", "dose", "dose_um", "pert_dose_unit"), names(profiles))
    time_fields <- intersect(c("pert_time", "time", "duration", "treatment_time"), names(profiles))
    id_fields <- intersect(c("sig_id", "signature_id", "inst_id", "distil_id"), names(profiles))
    profile_fields <- unique(c("pert", "cell", "type", "NCS", dose_fields, time_fields, id_fields))
    profile_export <- thp[, profile_fields, drop = FALSE]
    profile_export$dose_time_linkage <- if (length(dose_fields) && length(time_fields)) {
      "fields_present_review_values_units_and_lineage"
    } else "not_available_in_executed_result"
    result_rows <- lapply(seq_len(nrow(cand)), function(i) {
      p <- cand$pert[i]
      g <- raw$GSE70138[raw$GSE70138$pert_id == p, , drop = FALSE]
      b <- raw$compoundinfo_beta[raw$compoundinfo_beta$pert_id == p, , drop = FALSE]
      gn <- values(column(g, "pert_iname"))
      bn <- values(c(column(b, "cmap_name"), column(b, "pert_iname")))
      source_names <- values(c(gn, bn))
      keys <- values(c(column(g, "inchi_key"), column(b, "inchi_key")))
      keep_name <- clean(cand$pert_iname[i])
      keep_key <- clean(if ("inchi_key" %in% names(cand)) cand$inchi_key[i] else NA_character_)
      keep_smiles <- clean(if ("canonical_smiles" %in% names(cand)) cand$canonical_smiles[i] else NA_character_)
      d <- profiles[valid & profiles$pert == p & profiles$type == "trt_cp", , drop = FALSE]
      data.frame(
        pert = p, retained_name = keep_name, integrated_average_tie_rank = cand$objective_robustness_rank[i],
        primary_NCS = cand$primary_NCS[i],
        saved_summary_name = clean(compound$pert_iname[i]),
        summary_name_agrees = identical(clean(compound$pert_iname[i]), keep_name),
        GSE70138_rows_n = nrow(g), beta_rows_n = nrow(b),
        GSE70138_names = collapse(gn), beta_names = collapse(bn),
        distinct_source_names_n = length(unique(equal_name(source_names))),
        retained_name_supported_verbatim = !is.na(keep_name) && equal_name(keep_name) %in% equal_name(source_names),
        retained_inchi_key = keep_key, source_inchi_keys = collapse(keys),
        conflicting_source_inchi_keys = length(keys) > 1L,
        retained_key_differs_from_sources = !is.na(keep_key) && length(keys) > 0L && !keep_key %in% keys,
        retained_smiles_available = !is.na(keep_smiles),
        source_phase_or_status = collapse(c(column(b, "clinical_phase"), column(b, "phase"),
                                            column(b, "clinical_status"), column(b, "status"))),
        source_moa = collapse(c(column(b, "moa"), column(g, "moa"))),
        source_targets = collapse(c(column(b, "target"), column(g, "target"))),
        primary_reference_trt_cp_profiles_n = nrow(d),
        primary_reference_cell_lines_n = length(unique(d$cell)),
        duplicate_compound_cell_profiles = anyDuplicated(paste(d$cell, d$type, sep = "::")) > 0L,
        thp1_profiles_n = 1L,
        regulatory_status_verified = FALSE,
        dose_time_linkage = profile_export$dose_time_linkage[i],
        stringsAsFactors = FALSE
      )
    })
    all_compounds <- do.call(rbind, result_rows)
    key_counts <- table(stats::na.omit(all_compounds$retained_inchi_key))
    name_counts <- table(stats::na.omit(equal_name(all_compounds$retained_name)))
    all_compounds$BRD_ids_sharing_retained_key_n <- as.integer(key_counts[all_compounds$retained_inchi_key])
    all_compounds$BRD_ids_sharing_retained_name_n <- as.integer(name_counts[equal_name(all_compounds$retained_name)])
    all_compounds$identity_review_needed <- with(all_compounds,
      !summary_name_agrees | !retained_name_supported_verbatim | conflicting_source_inchi_keys |
      retained_key_differs_from_sources | distinct_source_names_n > 1L |
      (!is.na(BRD_ids_sharing_retained_key_n) & BRD_ids_sharing_retained_key_n > 1L) |
      (!is.na(BRD_ids_sharing_retained_name_n) & BRD_ids_sharing_retained_name_n > 1L))
    top20 <- all_compounds$integrated_average_tie_rank <= 20L
    priority <- top20 |
      all_compounds$pert %in% c("BRD-K82941592", "BRD-K03194791")
    review_queue <- all_compounds[priority | all_compounds$identity_review_needed, ]
    for (nm in c("identity_evidence_url", "regulatory_jurisdiction", "regulatory_status",
                  "regulatory_evidence_url", "verified_on", "dose_time_link_evidence", "resolution_note"))
      review_queue[[nm]] <- NA_character_
    # Header-only inventory: do not equate a raw signature with the executed exemplar.
    possible <- list.files("data/raw", pattern = "sig[_-]?info|inst[_-]?info", recursive = TRUE,
                            full.names = TRUE, ignore.case = TRUE)
    possible <- possible[grepl("\\.(txt|tsv|csv)(\\.gz)?$", possible, ignore.case = TRUE)]
    signature_inventory <- if (length(possible)) do.call(rbind, lapply(possible, function(p) {
      con <- if (grepl("\\.gz$", p, ignore.case = TRUE)) gzfile(p, "rt") else file(p, "rt")
      on.exit(close(con), add = TRUE)
      header <- readLines(con, n = 1L, warn = FALSE)
      data.frame(file = p, size_bytes = file.info(p)$size, header = paste(header, collapse = ""),
                   linkage_status = "header_inventoried_mapping_to_executed_profiles_not_established")
    })) else data.frame(file = character(), size_bytes = numeric(), header = character(), linkage_status = character())
    inventory <- do.call(rbind, source_inventory)
    check(identical(md5_before, unname(tools::md5sum(required))), "An input changed during the metadata audit.")
    base_out <- "results/revision/compound_metadata_audit"
    dir.create(base_out, recursive = TRUE, showWarnings = FALSE)
    out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = base_out)
    check(dir.create(out), "Cannot create a fresh metadata audit output folder.")
    csv <- function(d, stem) readr::write_csv(d, file.path(out, paste0(stem, ".csv")), na = "")
    csv(all_compounds, "all_compound_metadata_audit")
    csv(all_compounds[top20, ], "top20_compound_metadata_audit")
    csv(review_queue, "compound_metadata_review_queue")
    csv(inventory, "annotation_source_inventory")
    csv(signature_inventory, "signature_metadata_header_inventory")
    csv(profile_export, "executed_thp1_profile_metadata")
    for (src in names(source_rows)) csv(source_rows[[src]], paste0("raw_relevant_rows_", src))
    notes <- c(
      BUILD, "LOCAL METADATA AUDIT COMPLETE; EXTERNAL IDENTITY/REGULATORY REVIEW REMAINS PENDING.",
      "No names, structures, rankings, or result profiles have been changed; no LINCS query was executed.",
      "All relevant annotation rows are retained, including duplicate BRD IDs previously reduced to a first row.",
      "Different names can be legitimate synonyms. Distinct InChIKeys or shared IDs require source review, not automatic merging.",
      "A retained name absent verbatim from local sources may reflect a documented manual curation; check the evidence.",
      "Clinical phase/status fields from Broad are descriptive source metadata, not current regulatory-approval verification.",
      "All regulatory_status_verified fields are FALSE until dated authoritative evidence is recorded by jurisdiction.",
      "Dose/time fields in the executed result, if present, still require units and provenance review.",
      "Raw signature metadata cannot establish exemplar dose/time by a compound-only join.",
      "An absent dose/time linkage must be reported as unavailable, not filled from an arbitrary matching exposure.",
      "The editable review queue prioritizes the top 20 and every flagged identity record; flags are not proven errors.",
      paste("Average-tie summary used:", summary_path), paste("Locked manifest:", LOCK)
    )
    writeLines(notes, file.path(out, "metadata_audit_notes.txt"))
    result <- list(analysis_id = "revision_compound_metadata_audit_v1", build_id = BUILD,
      status = "PASS_LOCAL_AUDIT_EXTERNAL_REVIEW_PENDING", manifest_content_sha256 = LOCK,
      source_summary = summary_path, all_compounds = all_compounds, review_queue = review_queue,
      source_inventory = inventory, signature_metadata_inventory = signature_inventory,
      executed_thp1_profiles = profile_export,
      input_provenance = data.frame(input = required, md5 = md5_before), notes = notes)
    path <- file.path(out, "revision_compound_metadata_audit.rds")
    saveRDS(result, path)
    check(identical(readRDS(path), result), "Metadata audit read-back failed.")
    exported <- list.files(out, full.names = TRUE)
    check(all(file.info(exported)$size > 0), "An audit export is empty.")
    csv(data.frame(file = basename(exported), md5 = unname(tools::md5sum(exported))), "output_checksums")
    writeLines(c("PASS_LOCAL_AUDIT_EXTERNAL_REVIEW_PENDING", BUILD), file.path(out, "_SUCCESS.txt"))
    cat("\nREVISION ANALYSIS 04F COMPLETE\n")
    cat("Compounds audited:", nrow(all_compounds), "\n")
    cat("Identity records flagged for review:", sum(all_compounds$identity_review_needed), "\n")
    cat("Top-20 identity records flagged:", sum(all_compounds$identity_review_needed[top20]), "\n")
    cat("Annotation sources available:", sum(source_present), "/", length(sources), "\n")
    cat("Dose fields in executed result:", paste(dose_fields, collapse = ", "), "\n")
    cat("Time fields in executed result:", paste(time_fields, collapse = ", "), "\n")
    cat("Current regulatory approval independently verified: NO (review queue exported)\n")
    cat("Inputs unchanged: PASS\nOutput directory:", out, "\n")
    cat("Elapsed minutes:", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), "\n")
    invisible(result)
  }
  run()
})
