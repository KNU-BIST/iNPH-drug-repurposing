# LINCS2 exemplar provenance audit. Run from the project root.
# Source: the signatureSearchData vignette, section 5.2.1.
# This script reads saved scores; it does not run connectivity queries.
# Output: results/revision/lincs2_reference_provenance/run_*/
# The source metadata download is approximately 465 MB (once, then reused).

local({
  needed <- c("data.table", "digest", "jsonlite", "signatureSearch")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Required installed packages missing: ",
                           paste(absent, collapse = ", "))
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  primary_path <- file.path(root, "data/processed/lincs2_gess_result.rds")
  lock_path <- file.path(root, "renv.lock")
  if (!file.exists(primary_path) || !file.exists(lock_path)) {
    stop("Run from the inph-drug-repurposing project root with the saved ",
         "primary result and renv.lock present.")
  }
  source_url <- paste0("https://s3.amazonaws.com/macchiato.clue.io/",
                       "builds/LINCS2020/siginfo_beta.txt")
  expected_sha256 <- "1a38d7ea2a804be79804af4a27aff9f2537af8f13d3c8af1fdc1fef780f40201"
  expected_bytes <- 465242319
  metadata_path <- Sys.getenv("INPH_LINCS_SIGINFO_FILE", unset =
    file.path(root, "data/raw/LINCS2020/siginfo_beta.txt"))
  metadata_path <- path.expand(metadata_path)
  sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
  write_csv <- function(x, path) {
    utils::write.csv(as.data.frame(x), path, row.names = FALSE, na = "",
                     fileEncoding = "UTF-8")
  }
  if (!file.exists(metadata_path)) {
    dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
    partial <- tempfile("siginfo_download_", tmpdir = dirname(metadata_path))
    old_timeout <- getOption("timeout")
    options(timeout = max(1800, old_timeout))
    cat("Downloading the official LINCS2020 signature metadata (about 465 MB)...\n")
    tryCatch({
      status <- utils::download.file(source_url, partial, mode = "wb",
                                    method = "libcurl", quiet = FALSE)
      if (status != 0L || file.info(partial)$size != expected_bytes ||
          !identical(sha256(partial), expected_sha256)) {
        stop("Source size/checksum differs from the verified LINCS2020 file. ",
             "Inspect the source before using a different build.")
      }
      if (!file.rename(partial, metadata_path)) stop("Cannot save metadata file.")
    }, finally = {
      options(timeout = old_timeout)
      if (file.exists(partial)) unlink(partial)
    })
  }
  if (file.info(metadata_path)$size != expected_bytes ||
      !identical(sha256(metadata_path), expected_sha256)) {
    stop("Existing siginfo metadata does not match the verified source checksum. ",
         "It was not overwritten. Check its provenance before proceeding.")
  }
  cat("Reading verified metadata and the saved primary result...\n")
  meta <- data.table::fread(metadata_path, sep = "\t", quote = "",
    colClasses = "character", na.strings = c("", "NA"),
    showProgress = interactive(), data.table = FALSE)
  required <- c("pert_id", "cell_iname", "pert_type", "sig_id", "cmap_name",
    "is_exemplar_sig", "pert_dose", "pert_dose_unit", "pert_time",
    "pert_time_unit", "nsample", "distil_ids", "tas", "cc_q75",
    "ss_ngene", "is_hiq", "qc_pass")
  if (!all(required %in% names(meta))) stop("Required source columns missing.")
  exemplar <- meta[which(meta$pert_type == "trt_cp" &
                           meta$is_exemplar_sig == "1"), , drop = FALSE]
  exemplar$reference_column_id <- paste(exemplar$pert_id, exemplar$cell_iname,
                                       exemplar$pert_type, sep = "__")
  if (nrow(exemplar) != 136460L || anyDuplicated(exemplar$reference_column_id)) {
    stop("Source exemplar membership is not the expected unique 136,460 entries.")
  }
  saved <- readRDS(primary_path)
  primary <- if (is.data.frame(saved)) saved else
    as.data.frame(signatureSearch::result(saved))
  if (!all(c("pert", "cell", "type", "NCS", "WTCS", "WTCS_FDR") %in%
           names(primary))) stop("Saved primary result has unexpected fields.")
  primary$key <- paste(primary$pert, primary$cell, primary$type, sep = "__")
  membership_matches <- nrow(primary) == 136460L &&
    !anyDuplicated(primary$key) && setequal(primary$key, exemplar$reference_column_id)

  parent <- file.path(root, "results/revision/lincs2_reference_provenance")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
                  tmpdir = parent)
  if (!dir.create(out)) stop("Cannot create an audit output directory.")
  if (!membership_matches) {
    write_csv(data.frame(reference_column_id =
      setdiff(primary$key, exemplar$reference_column_id)),
      file.path(out, "saved_keys_without_source_exemplar.csv"))
    write_csv(data.frame(reference_column_id =
      setdiff(exemplar$reference_column_id, primary$key)),
      file.path(out, "source_exemplars_without_saved_key.csv"))
    stop("Full reference membership differs. Mismatch tables saved to: ", out,
         ". No automatic remapping or query rerun was performed.")
  }
  thp <- primary[which(primary$cell == "THP1" & primary$type == "trt_cp"),
                 , drop = FALSE]
  if (nrow(thp) != 1376L) stop("Unexpected THP1 reference size.")
  mapped <- exemplar[match(thp$key, exemplar$reference_column_id), , drop = FALSE]
  mapped$primary_NCS <- thp$NCS
  mapped$primary_WTCS <- thp$WTCS
  mapped$primary_WTCS_FDR <- thp$WTCS_FDR
  mapped$distil_ids_count <- vapply(mapped$distil_ids, function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_integer_)
    length(unique(strsplit(x, "|", fixed = TRUE)[[1]]))
  }, integer(1))
  mapped$replicate_count_agrees <- mapped$distil_ids_count ==
    suppressWarnings(as.integer(mapped$nsample))
  mapped$mapping_basis <- "Source exemplar flag and exact perturbagen/cell/type match"
  mapped$resource_id <- "EH7297"
  ids <- c("BRD-K82941592", "BRD-K03194791", "BRD-K53581288")
  target <- mapped[match(ids, mapped$pert_id), , drop = FALSE]
  if (anyNA(target$pert_id) || !all(target$replicate_count_agrees)) {
    stop("Candidate mapping or replicate identifiers require inspection.")
  }
  conditions <- meta[which(meta$pert_id %in% ids & meta$cell_iname == "THP1" &
                              meta$pert_type == "trt_cp"), , drop = FALSE]
  # Preserve the annotated exemplar even when another condition has higher TAS.
  conditions$maximum_TAS_for_compound_THP1 <- ave(
    as.numeric(conditions$tas), conditions$pert_id,
    FUN = function(x) max(x, na.rm = TRUE))
  conditions$maximum_TAS_record <- as.numeric(conditions$tas) ==
    conditions$maximum_TAS_for_compound_THP1
  replicate_rows <- do.call(rbind, lapply(seq_len(nrow(target)), function(i) {
    reps <- strsplit(target$distil_ids[i], "|", fixed = TRUE)[[1]]
    data.frame(compound = target$cmap_name[i], pert_id = target$pert_id[i],
      sig_id = target$sig_id[i], replicate_index = seq_along(reps),
      distil_id = reps, stringsAsFactors = FALSE)
  }))
  write_csv(target, file.path(out, "lincs2_thp1_candidate_exemplars.csv"))
  write_csv(conditions, file.path(out, "lincs2_thp1_candidate_all_conditions.csv"))
  write_csv(replicate_rows, file.path(out, "lincs2_thp1_candidate_replicate_ids.csv"))
  write_csv(mapped, file.path(out, "lincs2_all_thp1_exemplar_metadata.csv"))
  write_csv(exemplar[, c("reference_column_id", "sig_id")],
            file.path(out, "lincs2_full_reference_exemplar_map.csv"))
  write_csv(as.data.frame(table(qc_pass = mapped$qc_pass, is_hiq = mapped$is_hiq)),
            file.path(out, "lincs2_thp1_quality_flag_counts.csv"))
  requested_fields <- setdiff(required, c("pert_id", "cell_iname", "pert_type"))
  availability <- do.call(rbind, lapply(seq_len(nrow(target)), function(i) {
    values <- unlist(target[i, requested_fields], use.names = FALSE)
    data.frame(compound = target$cmap_name[i], field = requested_fields,
      value = values, available = !is.na(values) & nzchar(values),
      stringsAsFactors = FALSE)
  }))
  write_csv(availability, file.path(out, "lincs2_candidate_metadata_availability.csv"))
  dictionary_url <- paste0("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/",
    "LINCS2020%20Release%20Metadata%20Field%20Definitions.xlsx")
  write_csv(data.frame(field = c("nsample", "cc_q75", "tas", "qc_pass", "is_hiq"),
    meaning = c(
      "Count of Level-4 replicate profiles combined into a Level-5 signature.",
      "75th percentile of pairwise Spearman correlations between replicate landmark profiles.",
      "Normalized transcriptional activity measure combining signature strength and replicate concordance.",
      "Equals 1 when at least half the constituent replicates pass source processing QC.",
      "Equals 1 when qc_pass is 1 and at least one median recall-rank measure (Spearman or WTCS-50) is at most 5."),
    source_url = dictionary_url), file.path(out, "source_quality_field_definitions.csv"))

  lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
  packages <- c("signatureSearch", "signatureSearchData", "ExperimentHub",
                "AnnotationHub", "rhdf5", "data.table", "digest", "jsonlite")
  version_rows <- lapply(packages, function(p) {
    locked <- lock$Packages[[p]]$Version
    if (is.null(locked)) locked <- NA_character_
    installed <- if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else NA_character_
    data.frame(component = p, renv_lock_version = locked,
      current_audit_version = installed, stringsAsFactors = FALSE)
  })
  versions <- rbind(data.frame(component = "R", renv_lock_version = lock$R$Version,
    current_audit_version = as.character(getRversion())), do.call(rbind, version_rows))
  versions$record_scope <- "Declared lockfile versus current audit session; not proof of historical execution versions"
  write_csv(versions, file.path(out, "software_environment_records.csv"))
  capture.output(utils::sessionInfo(), file = file.path(out, "metadata_audit_session_info.txt"))
  write_csv(data.frame(source_url = source_url, resource_id = "EH7297",
    source_release = "LINCS2020 beta", source_bytes = file.info(metadata_path)$size,
    source_sha256 = expected_sha256, audit_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    signature_rows = nrow(meta), chemical_exemplars = nrow(exemplar),
    saved_reference_membership_matches = membership_matches,
    thp1_profiles = nrow(mapped), primary_result_sha256 = sha256(primary_path),
    lockfile_sha256 = sha256(lock_path)), file.path(out, "metadata_mapping_audit.csv"))
  writeLines(c(
    "Reference: LINCS2020 beta exemplar collection exposed as ExperimentHub EH7297.",
    "Replicate aggregation within a condition is distinct from selection across conditions.",
    "Mapping used is_exemplar_sig == 1; it did not reselect the current maximum-TAS record.",
    "All source quality fields are retained without filtering or recoding.",
    "Report is_hiq, qc_pass, TAS and replicate information; computational consistency does not establish reference-experiment quality.",
    "Quality-field meanings are recorded from the official LINCS2020 metadata dictionary.",
    "Signature qc_pass requires at least half its replicate profiles to pass processing QC; is_hiq additionally requires the stated recall criterion.",
    "Specific reasons for individual replicate QC failures are not recovered by this audit.",
    "Exact reference membership was checked against the saved full primary result.",
    "This is an identifier-based provenance mapping, not a comparison of Level-5 expression vectors.",
    "The current audit session and renv.lock are recorded separately; no historical session is invented.",
    "No connectivity queries, primary scores, candidate rankings, or manuscript figures were changed."
  ), file.path(out, "metadata_audit_notes.txt"))
  generated <- list.files(out, full.names = TRUE)
  write_csv(data.frame(file = basename(generated),
    sha256 = vapply(generated, sha256, character(1))), file.path(out, "output_checksums.csv"))
  writeLines("PASS_REFERENCE_MEMBERSHIP_AUDIT; SOURCE_QUALITY_FLAGS_RETAINED_FOR_REPORTING",
             file.path(out, "_SUCCESS.txt"))
  cat("\nMetadata audit completed. Output folder:\n", out, "\n", sep = "")
  print(target[, c("cmap_name", "sig_id", "pert_dose", "pert_dose_unit",
                   "pert_time", "pert_time_unit", "nsample", "tas", "is_hiq", "qc_pass")],
        row.names = FALSE)
})
