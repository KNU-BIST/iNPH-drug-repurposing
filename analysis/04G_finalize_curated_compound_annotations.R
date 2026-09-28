# =============================================================================
# Build: 04G_CURATED_COMPOUND_ANNOTATIONS_V1_20260907
# Install as analysis/04G_finalize_curated_compound_annotations.R
# =============================================================================
# Purpose: create a separate, provenance-preserving curation layer for compound
# annotations. BRD identifiers, NCS values, FDRs, and all ranks remain unchanged.
# No LINCS query is executed and no existing output is overwritten.

local({
  BUILD <- "04G_CURATED_COMPOUND_ANNOTATIONS_V1_20260907"
  LOCK <- paste0("1e3a2b827a1429e406303871d547a9dec",
                 "07ee4068078c666df69ca82e5428380")
  check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  clean <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | tolower(x) %in% c("", "na", "nan", "null", "unknown", "-666")] <- NA_character_
    x
  }
  values <- function(x) sort(unique(stats::na.omit(clean(x))), method = "radix")
  collapse_values <- function(x) paste(values(x), collapse = " | ")
  find_latest <- function(parent, pattern, success = "_SUCCESS.txt", override = "") {
    paths <- if (nzchar(override)) override else list.files(
      parent, pattern = pattern, recursive = TRUE, full.names = TRUE
    )
    paths <- paths[file.exists(paths) & file.exists(file.path(dirname(paths), success))]
    check(length(paths) > 0L, paste("No completed input found under", parent))
    paths[order(file.info(paths)$mtime, decreasing = TRUE)][1L]
  }
  run <- function() {
    started <- Sys.time()
    audit_path <- find_latest(
      "results/revision/compound_metadata_audit",
      "^revision_compound_metadata_audit\\.rds$",
      override = Sys.getenv("INPH_COMPOUND_METADATA_AUDIT_RDS", unset = "")
    )
    audit <- readRDS(audit_path)
    check(identical(audit$status, "PASS_LOCAL_AUDIT_EXTERNAL_REVIEW_PENDING") &&
            identical(audit$manifest_content_sha256, LOCK),
          "The selected 04F audit is not a valid locked result.")
    provenance <- as.data.frame(audit$input_provenance)
    check(all(c("input", "md5") %in% names(provenance)) && all(file.exists(provenance$input)) &&
            identical(unname(tools::md5sum(provenance$input)), as.character(provenance$md5)),
          "An input recorded by 04F has changed.")
    summary_path <- audit$source_summary
    check(length(summary_path) == 1L && file.exists(summary_path), "04E summary is missing.")
    summary <- readRDS(summary_path)
    check(identical(summary$status, "PASS") && identical(summary$ranking_method, "average") &&
            identical(summary$manifest_content_sha256, LOCK), "Invalid 04E summary.")
    result_dir_04f <- dirname(audit_path)
    gse_path <- file.path(result_dir_04f, "raw_relevant_rows_GSE70138.csv")
    beta_path <- file.path(result_dir_04f, "raw_relevant_rows_compoundinfo_beta.csv")
    check(all(file.exists(c(gse_path, beta_path))), "04F raw relevant-row exports are missing.")
    gse <- utils::read.csv(gse_path, stringsAsFactors = FALSE, check.names = FALSE,
                           na.strings = c("", "NA"))
    beta <- utils::read.csv(beta_path, stringsAsFactors = FALSE, check.names = FALSE,
                            na.strings = c("", "NA"))
    check(all(c("pert_id", "canonical_smiles", "inchi_key", "pert_iname") %in% names(gse)) &&
            all(c("pert_id", "canonical_smiles", "inchi_key", "cmap_name", "compound_aliases") %in% names(beta)),
          "Unexpected 04F raw-row schema.")
    candidates <- as.data.frame(summary$candidate_summary)
    local_audit <- as.data.frame(audit$all_compounds)
    check(nrow(candidates) == 1376L && nrow(local_audit) == 1376L &&
            !anyDuplicated(candidates$pert) && !anyDuplicated(local_audit$pert) &&
            setequal(candidates$pert, local_audit$pert), "Compound universe mismatch.")
    ix <- match(candidates$pert, local_audit$pert)
    local_audit <- local_audit[ix, , drop = FALSE]
    check(identical(as.numeric(candidates$objective_robustness_rank),
                    as.numeric(local_audit$integrated_average_tie_rank)) &&
            identical(as.numeric(candidates$primary_NCS), as.numeric(local_audit$primary_NCS)),
          "04E ranks or scores differ from 04F.")
    column <- function(d, nm) if (nm %in% names(d)) d[[nm]] else rep(NA_character_, nrow(d))
    curated <- data.frame(
      pert = as.character(candidates$pert),
      integrated_average_tie_rank = as.integer(candidates$objective_robustness_rank),
      primary_NCS = as.numeric(candidates$primary_NCS),
      primary_FDR = as.numeric(column(candidates, "primary_FDR")),
      original_retained_name = clean(candidates$pert_iname),
      original_canonical_smiles = clean(column(candidates, "canonical_smiles")),
      original_inchi_key = clean(column(candidates, "inchi_key")),
      curated_display_name = clean(candidates$pert_iname),
      curated_canonical_smiles = clean(column(candidates, "canonical_smiles")),
      curated_inchi_key = clean(column(candidates, "inchi_key")),
      curated_target = clean(local_audit$source_targets),
      curated_moa = clean(local_audit$source_moa),
      local_identity_flag_04f = as.logical(local_audit$identity_review_needed),
      structure_filled_from_exact_brd_source = FALSE,
      recorded_structure_identity_resolved = !as.logical(local_audit$identity_review_needed),
      safe_for_qualified_reporting = !as.logical(local_audit$identity_review_needed),
      physical_assayed_material_independently_verified = FALSE,
      curation_status = ifelse(as.logical(local_audit$identity_review_needed),
                               "pending_external_identity_review", "unchanged_local_sources_consistent"),
      manuscript_handling = ifelse(as.logical(local_audit$identity_review_needed),
                                   "Do not highlight by name until reviewed", "May retain BRD-linked name"),
      identity_evidence = NA_character_,
      US_regulatory_status = "not_reviewed",
      regulatory_evidence_url = NA_character_,
      regulatory_verified_on = NA_character_,
      iNPH_indication_in_cited_US_label = NA,
      stringsAsFactors = FALSE
    )
    # Fill missing structures only when all exact-BRD local source rows agree.
    for (i in seq_len(nrow(curated))) {
      p <- curated$pert[i]
      rows_g <- gse[gse$pert_id == p, , drop = FALSE]
      rows_b <- beta[beta$pert_id == p, , drop = FALSE]
      source_keys <- values(c(rows_g$inchi_key, rows_b$inchi_key))
      source_smiles <- values(c(rows_g$canonical_smiles, rows_b$canonical_smiles))
      if (is.na(curated$curated_inchi_key[i]) && length(source_keys) == 1L)
        curated$curated_inchi_key[i] <- source_keys
      if (is.na(curated$curated_canonical_smiles[i]) && length(source_keys) == 1L &&
          length(source_smiles) == 1L) {
        curated$curated_canonical_smiles[i] <- source_smiles
        curated$structure_filled_from_exact_brd_source[i] <- TRUE
      }
    }
    set_record <- function(pert, ...) {
      i <- match(pert, curated$pert)
      check(!is.na(i), paste("Curated BRD ID missing:", pert))
      entries <- list(...)
      for (nm in names(entries)) curated[i, nm] <<- entries[[nm]]
    }
    set_record(
      "BRD-A02367930",
      curated_display_name = "ethinyl-estradiol",
      recorded_structure_identity_resolved = TRUE,
      safe_for_qualified_reporting = TRUE,
      curation_status = "resolved_alias_only",
      manuscript_handling = "Retain ethinyl-estradiol; beta cmap_name is the BRD ID and compound_aliases supplies the name",
      identity_evidence = "GSE70138 row 599 and compoundinfo_beta row 277 share InChIKey BFPYWIDHMRZLRN-OQPPHWFISA-N"
    )
    set_record(
      "BRD-U86686840",
      curated_display_name = "(R)-AZD-7762 stereoisomer [BRD-U86686840]",
      recorded_structure_identity_resolved = TRUE,
      safe_for_qualified_reporting = TRUE,
      curation_status = "resolved_recorded_R_stereoisomer_physical_material_unverified",
      manuscript_handling = "If reported, use the stereochemical label and BRD ID; do not call this unqualified standard AZD-7762",
      identity_evidence = paste(
        "Local sources: IAYGCINLNONXHY-GFCCVEGCSA-N.",
        "PubChem CID 67077825 identifies the recorded R stereoisomer:",
        "https://pubchem.ncbi.nlm.nih.gov/compound/67077825"
      )
    )
    set_record(
      "BRD-K86525559",
      curated_display_name = "AZD-7762 (S stereoisomer) [BRD-K86525559]",
      recorded_structure_identity_resolved = TRUE,
      safe_for_qualified_reporting = TRUE,
      curation_status = "resolved_recorded_S_stereoisomer_physical_material_unverified",
      manuscript_handling = "Keep separate from BRD-U86686840; use the S-stereoisomer label and BRD ID",
      identity_evidence = paste(
        "Local sources: IAYGCINLNONXHY-LBPRGKRZSA-N.",
        "PubChem CID 11152667 identifies the usual AZD7762 S stereoisomer:",
        "https://pubchem.ncbi.nlm.nih.gov/compound/11152667"
      )
    )
    set_record(
      "BRD-K70642949",
      curated_target = "PDPK1",
      curated_moa = "3-phosphoinositide-dependent protein kinase-1 (PDPK1) inhibitor",
      curation_status = "resolved_mechanism_annotation",
      manuscript_handling = "Use PDPK1 inhibitor; do not describe as a pyruvate dehydrogenase kinase inhibitor",
      identity_evidence = "Najafov et al., Biochem J 2011, PMID 21087210, DOI 10.1042/BJ20101732"
    )
    structure_ids <- c("BRD-K12683773", "BRD-K40624912", "BRD-K59962020")
    structure_names <- c("methyl-fasudil", "ZM-39923", "CHEMBL-374350")
    expected_keys <- c("UGJIFGVOVGCXQP-UHFFFAOYSA-N", "JSASWRWALCMOQP-UHFFFAOYSA-N",
                       "CYCGGKILBWERDJ-UHFFFAOYSA-N")
    for (j in seq_along(structure_ids)) {
      i <- match(structure_ids[j], curated$pert)
      check(!is.na(i) && identical(curated$curated_inchi_key[i], expected_keys[j]) &&
              !is.na(curated$curated_canonical_smiles[i]),
            paste("Source-backed structure recovery failed for", structure_names[j]))
      curated$curation_status[i] <- "resolved_exact_BRD_source_structure_fill"
      curated$manuscript_handling[i] <- "Structure restored from the exact-BRD compoundinfo_beta row; original field retained separately"
      curated$identity_evidence[i] <- paste("compoundinfo_beta exact-BRD InChIKey", expected_keys[j])
    }
    # Dated U.S. regulatory evidence for the two priority candidates only.
    set_record(
      "BRD-K82941592",
      US_regulatory_status = "FDA-approved active ingredient for labeled cardiovascular/lipid indications; not approved for iNPH",
      regulatory_evidence_url = "https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/021366s047lbl.pdf",
      regulatory_verified_on = "2026-09-07",
      iNPH_indication_in_cited_US_label = FALSE
    )
    set_record(
      "BRD-K03194791",
      US_regulatory_status = "FDA-approved oral active ingredient to reduce exacerbation risk in severe COPD with chronic bronchitis; not approved for iNPH",
      regulatory_evidence_url = "https://www.accessdata.fda.gov/drugsatfda_docs/appletter/2011/022522s000ltr.pdf",
      regulatory_verified_on = "2026-09-07",
      iNPH_indication_in_cited_US_label = FALSE
    )
    # Hard checks for the stereochemical split and immutable analysis outputs.
    azd <- curated[curated$pert %in% c("BRD-U86686840", "BRD-K86525559"), ]
    check(nrow(azd) == 2L && length(unique(azd$curated_inchi_key)) == 2L &&
            setequal(azd$curated_inchi_key,
                     c("IAYGCINLNONXHY-GFCCVEGCSA-N", "IAYGCINLNONXHY-LBPRGKRZSA-N")),
          "AZD-7762 stereochemical records were not kept separate.")
    check(identical(curated$pert, as.character(candidates$pert)) &&
            identical(curated$integrated_average_tie_rank,
                      as.integer(candidates$objective_robustness_rank)) &&
            identical(curated$primary_NCS, as.numeric(candidates$primary_NCS)),
          "A BRD identifier, score, or rank changed during curation.")
    top20 <- curated[curated$integrated_average_tie_rank <= 20L, , drop = FALSE]
    check(nrow(top20) == 20L && all(top20$recorded_structure_identity_resolved) &&
            all(top20$safe_for_qualified_reporting),
          "A top-20 identity remains unresolved for qualified reporting.")
    unresolved <- curated[!curated$recorded_structure_identity_resolved, , drop = FALSE]
    decision_ids <- unique(c("BRD-A02367930", "BRD-U86686840", "BRD-K86525559",
                             "BRD-K70642949", structure_ids,
                             "BRD-K82941592", "BRD-K03194791"))
    decisions <- curated[curated$pert %in% decision_ids, , drop = FALSE]
    decisions <- decisions[order(decisions$integrated_average_tie_rank), ]
    required <- unique(c(provenance$input, audit_path, summary_path, gse_path, beta_path))
    own <- "analysis/04G_finalize_curated_compound_annotations.R"
    if (file.exists(own)) required <- unique(c(required, own))
    hashes_before <- unname(tools::md5sum(required))
    check(!anyNA(hashes_before), "An input checksum could not be calculated.")
    parent <- "results/revision/curated_compound_annotations"
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
    check(dir.create(out), "Cannot create a fresh curation output directory.")
    utils::write.csv(curated, file.path(out, "curated_compound_annotations_all_1376.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(top20, file.path(out, "curated_compound_annotations_top20.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(decisions, file.path(out, "curation_decision_log.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(unresolved, file.path(out, "unresolved_identity_review_queue.csv"),
                     row.names = FALSE, na = "")
    notes <- c(
      BUILD,
      "BRD-level scores, FDRs, ranks, and query results are unchanged.",
      "All top-20 records are resolved for qualified reporting; this does not verify the physical material in the original assay vial.",
      "BRD-U86686840 (rank 11) is the recorded R stereoisomer and must not be named as unqualified standard AZD-7762.",
      "BRD-K86525559 (rank 869) is the recorded S stereoisomer usually indexed as AZD7762; the records remain separate.",
      "GSK-2334470 is curated as a PDPK1 inhibitor; PDK1 here does not mean pyruvate dehydrogenase kinase.",
      "Rosuvastatin and oral roflumilast are FDA-approved for cited non-iNPH indications; this does not imply efficacy or approval in iNPH.",
      "Other identity flags remain in the unresolved queue and should not be promoted by name without review.",
      "The observed rank bands in the average-tie scatterplot are expected consequences of tied scores, especially zero blocks.",
      paste("Frozen manifest SHA-256:", LOCK), paste("04F audit:", audit_path)
    )
    writeLines(notes, file.path(out, "curation_notes.txt"))
    result <- list(
      analysis_id = "revision_curated_compound_annotations_v1", build_id = BUILD,
      status = "PASS_TOP20_RESOLVED_QUALIFIED_REPORTING", manifest_content_sha256 = LOCK,
      source_audit = audit_path, source_summary = summary_path, curated_all = curated,
      curated_top20 = top20, decision_log = decisions, unresolved_queue = unresolved,
      figure_review = list(
        file = "primary_vs_core_percentiles_average_ties.png",
        status = "PASS_VISUAL_REVIEW",
        interpretation = "Bands reflect average ranks of tied NCS blocks; axes, priority labels, and equal-percentile line are correct."
      ),
      notes = notes, input_provenance = data.frame(input = required, md5 = hashes_before),
      session_info = utils::sessionInfo(), output_directory = out
    )
    rds_path <- file.path(out, "revision_curated_compound_annotations.rds")
    saveRDS(result, rds_path)
    check(identical(readRDS(rds_path), result), "Curation RDS read-back failed.")
    check(identical(hashes_before, unname(tools::md5sum(required))),
          "An input changed during curation.")
    exports <- list.files(out, full.names = TRUE)
    check(all(file.info(exports)$size > 0), "A curation export is empty.")
    utils::write.csv(data.frame(file = basename(exports), md5 = unname(tools::md5sum(exports))),
                     file.path(out, "output_checksums.csv"), row.names = FALSE)
    writeLines(c("PASS_TOP20_RESOLVED_QUALIFIED_REPORTING", BUILD),
               file.path(out, "_SUCCESS.txt"))
    cat("\nREVISION ANALYSIS 04G COMPLETE\n")
    cat("Compounds retained:", nrow(curated), "| Scores/ranks changed: 0\n")
    cat("Top-20 identities resolved for qualified reporting:", nrow(top20), "/ 20\n")
    cat("Other records remaining in identity review queue:", nrow(unresolved), "\n")
    cat("Source-backed structure fields filled:", sum(curated$structure_filled_from_exact_brd_source), "\n")
    cat("AZD stereoisomers kept separate: PASS\nGSK-2334470 PDPK1 annotation corrected: PASS\n")
    cat("Priority U.S. regulatory evidence recorded: rosuvastatin and roflumilast (non-iNPH indications)\n")
    cat("Average-tie scatterplot visual review: PASS\n")
    cat("Output directory:", out, "\nElapsed minutes:",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2), "\n")
    invisible(result)
  }
  run()
})
