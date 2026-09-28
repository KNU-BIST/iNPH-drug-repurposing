# =============================================================================
# Build: 04H_LINCS_EXPOSURE_METADATA_AUDIT_V1_20260907
# Install as analysis/04H_audit_lincs_exposure_metadata.R
# =============================================================================
# Purpose: determine whether the executed THP1 LINCS results can be linked
# unambiguously to dose and exposure duration. This is a read-only provenance
# audit. It does not run a LINCS query, alter a score or rank, or infer exposure
# metadata by an unsafe compound-only join.

local({
  BUILD <- "04H_LINCS_EXPOSURE_METADATA_AUDIT_V1_20260907"
  LOCK <- paste0(
    "1e3a2b827a1429e406303871d547a9dec",
    "07ee4068078c666df69ca82e5428380"
  )

  check <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }

  find_latest <- function(parent, pattern, override = "") {
    paths <- if (nzchar(override)) {
      override
    } else {
      list.files(
        parent,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE
      )
    }

    paths <- paths[
      file.exists(paths) &
        file.exists(file.path(dirname(paths), "_SUCCESS.txt"))
    ]

    check(length(paths) > 0L, paste("No completed input found under", parent))
    paths[order(file.info(paths)$mtime, decreasing = TRUE)][1L]
  }

  existing_fields <- function(column_names, candidates) {
    column_names[tolower(column_names) %in% tolower(candidates)]
  }

  collapse_fields <- function(x) {
    if (length(x)) paste(x, collapse = " | ") else ""
  }

  run <- function() {
    started <- Sys.time()

    curated_path <- find_latest(
      "results/revision/curated_compound_annotations",
      "^revision_curated_compound_annotations\\.rds$",
      override = Sys.getenv(
        "INPH_CURATED_COMPOUND_ANNOTATIONS_RDS",
        unset = ""
      )
    )

    curated_result <- readRDS(curated_path)

    check(
      identical(
        curated_result$status,
        "PASS_TOP20_RESOLVED_QUALIFIED_REPORTING"
      ) &&
        identical(curated_result$manifest_content_sha256, LOCK),
      "The selected 04G result is not a valid locked curation result."
    )

    curated_provenance <- as.data.frame(curated_result$input_provenance)
    check(
      all(c("input", "md5") %in% names(curated_provenance)) &&
        all(file.exists(curated_provenance$input)) &&
        identical(
          unname(tools::md5sum(curated_provenance$input)),
          as.character(curated_provenance$md5)
        ),
      "An input recorded by 04G has changed."
    )

    audit_path <- curated_result$source_audit
    check(
      length(audit_path) == 1L && file.exists(audit_path),
      "The source 04F audit is missing."
    )

    audit_result <- readRDS(audit_path)
    check(
      identical(
        audit_result$status,
        "PASS_LOCAL_AUDIT_EXTERNAL_REVIEW_PENDING"
      ) &&
        identical(audit_result$manifest_content_sha256, LOCK),
      "The source 04F audit is invalid."
    )

    audit_provenance <- as.data.frame(audit_result$input_provenance)
    check(
      all(c("input", "md5") %in% names(audit_provenance)) &&
        all(file.exists(audit_provenance$input)) &&
        identical(
          unname(tools::md5sum(audit_provenance$input)),
          as.character(audit_provenance$md5)
        ),
      "An input recorded by 04F has changed."
    )

    audit_directory <- dirname(audit_path)
    executed_export_path <- file.path(
      audit_directory,
      "executed_thp1_profile_metadata.csv"
    )
    inventory_path <- file.path(
      audit_directory,
      "signature_metadata_header_inventory.csv"
    )
    primary_result_path <- "data/processed/lincs2_gess_result.rds"

    check(
      all(file.exists(c(
        executed_export_path,
        inventory_path,
        primary_result_path
      ))),
      "One or more exposure-audit inputs are missing."
    )

    executed <- utils::read.csv(
      executed_export_path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA")
    )

    inventory <- utils::read.csv(
      inventory_path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA")
    )

    check(
      identical(
        names(executed),
        c("pert", "cell", "type", "NCS", "dose_time_linkage")
      ),
      "Unexpected executed THP1 metadata export schema."
    )
    check(
      identical(
        names(inventory),
        c("file", "size_bytes", "header", "linkage_status")
      ),
      "Unexpected signature metadata inventory schema."
    )

    curated <- as.data.frame(curated_result$curated_all)
    required_curated_columns <- c(
      "pert",
      "integrated_average_tie_rank",
      "primary_NCS",
      "curated_display_name"
    )

    check(
      all(required_curated_columns %in% names(curated)) &&
        nrow(curated) == 1376L &&
        !anyDuplicated(curated$pert),
      "The curated compound universe is invalid."
    )

    check(
      nrow(executed) == 1376L &&
        !anyDuplicated(executed$pert) &&
        setequal(executed$pert, curated$pert) &&
        all(toupper(executed$cell) == "THP1") &&
        all(executed$type == "trt_cp") &&
        !anyNA(executed$NCS),
      "The executed THP1 export does not define the expected 1,376 profiles."
    )

    executed <- executed[match(curated$pert, executed$pert), , drop = FALSE]
    check(
      identical(as.character(executed$pert), as.character(curated$pert)) &&
        identical(as.numeric(executed$NCS), as.numeric(curated$primary_NCS)),
      "Executed THP1 scores differ from the curated score universe."
    )

    primary_object <- readRDS(primary_result_path)
    check(
      isS4(primary_object) &&
        "result" %in% methods::slotNames(primary_object),
      "Expected the saved primary result to be a gessResult-like S4 object."
    )

    primary_table <- as.data.frame(methods::slot(primary_object, "result"))
    check(
      all(c("pert", "cell", "type", "NCS") %in% names(primary_table)),
      "The saved primary result is missing a required result column."
    )

    valid_primary <-
      !is.na(primary_table$pert) &
      !is.na(primary_table$cell) &
      !is.na(primary_table$type)

    primary_thp1 <- primary_table[
      valid_primary &
        toupper(primary_table$cell) == "THP1" &
        primary_table$type == "trt_cp" &
        primary_table$pert %in% curated$pert,
      ,
      drop = FALSE
    ]

    check(
      nrow(primary_thp1) == 1376L &&
        !anyDuplicated(primary_thp1$pert) &&
        setequal(primary_thp1$pert, curated$pert),
      "The direct primary result does not contain one THP1 row per compound."
    )

    primary_thp1 <- primary_thp1[
      match(curated$pert, primary_thp1$pert),
      ,
      drop = FALSE
    ]

    check(
      identical(
        as.numeric(primary_thp1$NCS),
        as.numeric(curated$primary_NCS)
      ),
      "Direct primary-result scores differ from the curated score universe."
    )

    default_rownames <- identical(
      rownames(primary_table),
      as.character(seq_len(nrow(primary_table)))
    )

    check(
      default_rownames,
      paste(
        "The primary result has non-default row names.",
        "Review whether they encode signature identifiers before continuing."
      )
    )

    dose_candidates <- c(
      "pert_dose",
      "dose",
      "dose_um",
      "concentration",
      "pert_dose_unit",
      "dose_unit"
    )
    duration_candidates <- c(
      "pert_time",
      "time",
      "duration",
      "treatment_time",
      "pert_time_unit",
      "time_unit"
    )
    signature_id_candidates <- c(
      "sig_id",
      "signature_id",
      "inst_id",
      "distil_id",
      "exemplar_id"
    )

    direct_dose_fields <- existing_fields(
      names(primary_table),
      dose_candidates
    )
    direct_duration_fields <- existing_fields(
      names(primary_table),
      duration_candidates
    )
    direct_signature_fields <- existing_fields(
      names(primary_table),
      signature_id_candidates
    )
    export_dose_fields <- existing_fields(names(executed), dose_candidates)
    export_duration_fields <- existing_fields(
      names(executed),
      duration_candidates
    )
    export_signature_fields <- existing_fields(
      names(executed),
      signature_id_candidates
    )

    check(
      length(direct_dose_fields) == 0L &&
        length(direct_duration_fields) == 0L &&
        length(direct_signature_fields) == 0L &&
        length(export_dose_fields) == 0L &&
        length(export_duration_fields) == 0L &&
        length(export_signature_fields) == 0L,
      paste(
        "An exposure or signature-level field is present.",
        "Review it directly instead of applying the unavailable-metadata conclusion."
      )
    )

    check(
      all(
        executed$dose_time_linkage ==
          "not_available_in_executed_result"
      ),
      "The prior audit did not classify every exposure linkage as unavailable."
    )

    schema_audit <- data.frame(
      source = c(
        "saved_primary_result",
        "saved_primary_result",
        "saved_primary_result",
        "04F_executed_THP1_export",
        "04F_executed_THP1_export",
        "04F_executed_THP1_export",
        "local_data_raw_inventory"
      ),
      field_family = c(
        "signature_or_instance_identifier",
        "dose_or_concentration",
        "exposure_duration",
        "signature_or_instance_identifier",
        "dose_or_concentration",
        "exposure_duration",
        "signature_or_instance_metadata_files"
      ),
      candidate_fields = c(
        collapse_fields(signature_id_candidates),
        collapse_fields(dose_candidates),
        collapse_fields(duration_candidates),
        collapse_fields(signature_id_candidates),
        collapse_fields(dose_candidates),
        collapse_fields(duration_candidates),
        "sig_info | inst_info"
      ),
      present_fields_or_files_n = c(
        length(direct_signature_fields),
        length(direct_dose_fields),
        length(direct_duration_fields),
        length(export_signature_fields),
        length(export_dose_fields),
        length(export_duration_fields),
        nrow(inventory)
      ),
      present_names = c(
        collapse_fields(direct_signature_fields),
        collapse_fields(direct_dose_fields),
        collapse_fields(direct_duration_fields),
        collapse_fields(export_signature_fields),
        collapse_fields(export_dose_fields),
        collapse_fields(export_duration_fields),
        collapse_fields(inventory$file)
      ),
      linkage_status = c(
        rep("not_available", 6L),
        if (nrow(inventory)) {
          "files_present_but_no_executed_profile_linkage"
        } else {
          "no_local_signature_metadata_files_found"
        }
      ),
      stringsAsFactors = FALSE
    )

    exposure_availability <- data.frame(
      pert = as.character(curated$pert),
      curated_display_name = as.character(curated$curated_display_name),
      integrated_average_tie_rank = as.integer(
        curated$integrated_average_tie_rank
      ),
      primary_NCS = as.numeric(curated$primary_NCS),
      executed_cell = as.character(executed$cell),
      executed_perturbation_type = as.character(executed$type),
      executed_signature_identifier_available = FALSE,
      dose_available = FALSE,
      duration_available = FALSE,
      dose = NA_real_,
      dose_unit = NA_character_,
      duration = NA_real_,
      duration_unit = NA_character_,
      linkage_status = "UNAVAILABLE_NO_SIGNATURE_LEVEL_LINKAGE",
      reporting_status = paste0(
        "Do not infer from a compound-only join; report exposure metadata ",
        "as unavailable in the saved processed result"
      ),
      stringsAsFactors = FALSE
    )

    top20_exposure <- exposure_availability[
      exposure_availability$integrated_average_tie_rank <= 20L,
      ,
      drop = FALSE
    ]

    check(
      nrow(top20_exposure) == 20L &&
        all(!top20_exposure$dose_available) &&
        all(!top20_exposure$duration_available),
      "Unexpected top-20 exposure-availability result."
    )

    methods_text <- paste(
      "Connectivity scores were obtained from the processed LINCS2",
      "reference accessed through signatureSearch and signatureSearchData.",
      "The saved search result contained one THP1 trt_cp result per BRD",
      "perturbagen, but did not retain a signature- or instance-level",
      "identifier, concentration, or exposure duration. Compound-specific",
      "dose and duration therefore could not be reconstructed without an",
      "unverifiable compound-only linkage to lower-level profiles and were",
      "not reported."
    )

    limitations_text <- paste(
      "The processed connectivity results do not retain dose or exposure-",
      "duration metadata. Consequently, these predictions identify",
      "perturbagens whose reference transcriptional programs oppose the",
      "iNPH signature, but they do not define an experimentally validated",
      "dose or treatment duration."
    )

    reviewer_text <- paste(
      "We audited the saved THP1 connectivity output and the available local",
      "metadata. All 1,376 THP1 trt_cp results were represented once at the",
      "BRD perturbagen level, but neither the saved result nor its local",
      "metadata retained concentration, exposure duration, or a signature-",
      "or instance-level identifier that could establish the executed",
      "exposure. Because joining by compound alone could attach an arbitrary",
      "LINCS exposure, we did not infer these fields. We now report this",
      "limitation explicitly and provide the schema and provenance audit.",
      "This clarification changes no NCS value, FDR, or compound rank."
    )

    manuscript_and_reviewer_text <- c(
      "METHODS ADDITION",
      methods_text,
      "",
      "LIMITATIONS ADDITION",
      limitations_text,
      "",
      "RESPONSE TO REVIEWER",
      reviewer_text
    )

    documentation <- data.frame(
      resource = c(
        "signatureSearch Bioconductor vignette",
        "signatureSearchData Bioconductor package manual",
        "signatureSearchData source repository"
      ),
      url = c(
        paste0(
          "https://www.bioconductor.org/packages/release/bioc/vignettes/",
          "signatureSearch/inst/doc/signatureSearch.html"
        ),
        paste0(
          "https://www.bioconductor.org/packages/release/data/experiment/",
          "manuals/signatureSearchData/man/signatureSearchData.pdf"
        ),
        "https://github.com/girke-lab/signatureSearchData"
      ),
      relevance = c(
        "Documents connectivity searches against processed reference signatures",
        "Documents the processed experiment-data resources used by signatureSearch",
        "Describes the packaged LINCS reference resources"
      ),
      stringsAsFactors = FALSE
    )

    required_inputs <- unique(c(
      curated_provenance$input,
      audit_provenance$input,
      curated_path,
      audit_path,
      executed_export_path,
      inventory_path,
      primary_result_path
    ))
    own_script <- "analysis/04H_audit_lincs_exposure_metadata.R"
    if (file.exists(own_script)) {
      required_inputs <- unique(c(required_inputs, own_script))
    }

    input_hashes_before <- unname(tools::md5sum(required_inputs))
    check(
      !anyNA(input_hashes_before),
      "An exposure-audit input checksum could not be calculated."
    )

    result <- list(
      analysis_id = "revision_lincs_exposure_metadata_audit_v1",
      build_id = BUILD,
      status = "PASS_EXPOSURE_METADATA_UNAVAILABLE_NO_UNSUPPORTED_LINKAGE",
      manifest_content_sha256 = LOCK,
      source_curated_result = curated_path,
      source_metadata_audit = audit_path,
      executed_result_granularity = paste(
        "One saved THP1 trt_cp result per BRD perturbagen;",
        "no signature- or instance-level identifier retained"
      ),
      primary_result_columns = names(primary_table),
      primary_result_default_rownames = default_rownames,
      exposure_availability = exposure_availability,
      top20_exposure_availability = top20_exposure,
      schema_audit = schema_audit,
      local_signature_metadata_inventory = inventory,
      documentation = documentation,
      manuscript_text = list(
        methods = methods_text,
        limitations = limitations_text,
        response_to_reviewer = reviewer_text
      ),
      score_and_rank_changes_n = 0L,
      new_lincs_queries_n = 0L,
      input_provenance = data.frame(
        input = required_inputs,
        md5 = input_hashes_before,
        stringsAsFactors = FALSE
      ),
      session_info = utils::sessionInfo()
    )

    check(
      identical(
        input_hashes_before,
        unname(tools::md5sum(required_inputs))
      ),
      "An input changed during the exposure-metadata audit."
    )

    output_parent <- "results/revision/lincs_exposure_metadata_audit"
    dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
    output_directory <- tempfile(
      paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
      tmpdir = output_parent
    )
    check(
      dir.create(output_directory),
      "Cannot create a fresh exposure-metadata audit directory."
    )

    result$output_directory <- output_directory

    utils::write.csv(
      exposure_availability,
      file.path(
        output_directory,
        "executed_thp1_exposure_availability_all_1376.csv"
      ),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      top20_exposure,
      file.path(
        output_directory,
        "executed_thp1_exposure_availability_top20.csv"
      ),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      schema_audit,
      file.path(output_directory, "executed_result_schema_audit.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      inventory,
      file.path(output_directory, "local_signature_metadata_inventory.csv"),
      row.names = FALSE,
      na = ""
    )
    utils::write.csv(
      documentation,
      file.path(output_directory, "documentation_sources.csv"),
      row.names = FALSE,
      na = ""
    )
    writeLines(
      manuscript_and_reviewer_text,
      file.path(output_directory, "manuscript_and_reviewer_text.txt")
    )

    notes <- c(
      BUILD,
      "PASS: EXPOSURE METADATA UNAVAILABLE; NO UNSUPPORTED LINKAGE PERFORMED.",
      "No LINCS query was run. No NCS, FDR, BRD identifier, annotation, or rank was changed.",
      paste("Executed THP1 compounds audited:", nrow(exposure_availability)),
      "Saved signature/instance identifiers available: 0",
      "Saved dose/concentration fields available: 0",
      "Saved exposure-duration fields available: 0",
      paste("Local sig_info/inst_info files inventoried:", nrow(inventory)),
      "A compound-only join is insufficient because a perturbagen can have multiple doses, durations, cell lines, and instances.",
      "The absence of exposure metadata limits experimental interpretation but does not invalidate the saved perturbagen-level connectivity score.",
      paste("Frozen manifest SHA-256:", LOCK),
      paste("04G curated result:", curated_path),
      paste("04F metadata audit:", audit_path)
    )
    writeLines(
      notes,
      file.path(output_directory, "exposure_metadata_audit_notes.txt")
    )

    result_path <- file.path(
      output_directory,
      "revision_lincs_exposure_metadata_audit.rds"
    )
    saveRDS(result, result_path)
    check(
      identical(readRDS(result_path), result),
      "Exposure-metadata audit RDS read-back failed."
    )

    exports <- list.files(output_directory, full.names = TRUE)
    check(
      all(file.info(exports)$size > 0),
      "An exposure-metadata audit export is empty."
    )

    utils::write.csv(
      data.frame(
        file = basename(exports),
        md5 = unname(tools::md5sum(exports)),
        stringsAsFactors = FALSE
      ),
      file.path(output_directory, "output_checksums.csv"),
      row.names = FALSE
    )

    writeLines(
      c(
        "PASS_EXPOSURE_METADATA_UNAVAILABLE_NO_UNSUPPORTED_LINKAGE",
        BUILD
      ),
      file.path(output_directory, "_SUCCESS.txt")
    )

    check(
      identical(
        input_hashes_before,
        unname(tools::md5sum(required_inputs))
      ),
      "An input changed while exposure-audit outputs were written."
    )

    cat("\nREVISION ANALYSIS 04H COMPLETE\n")
    cat("Status: PASS_EXPOSURE_METADATA_UNAVAILABLE_NO_UNSUPPORTED_LINKAGE\n")
    cat("Executed THP1 compounds matched:", nrow(exposure_availability), "/ 1376\n")
    cat("Signature/instance identifiers retained: 0\n")
    cat("Dose/concentration fields retained: 0\n")
    cat("Exposure-duration fields retained: 0\n")
    cat("Local sig_info/inst_info files found:", nrow(inventory), "\n")
    cat("Unsupported compound-only joins performed: 0\n")
    cat("Scores or ranks changed: 0\n")
    cat("New LINCS queries run: 0\n")
    cat("Manuscript/reviewer wording exported: PASS\n")
    cat("Output directory:", output_directory, "\n")
    cat(
      "Elapsed minutes:",
      round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
      "\n"
    )

    invisible(result)
  }

  run()
})
