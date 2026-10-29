# Self-reported diagnosis dates from the separate corrected exports.
# Match Participant ID, Instance and Array; retain unknown-year flags.

# 1. Packages and local paths ----
library(data.table)
library(survival)
source_dir <- Sys.getenv('UKB_SOURCE_DIR', unset = file.path(getwd(), 'data', 'raw'))
results_dir <- Sys.getenv('UKB_RESULTS_DIR', unset = file.path(getwd(), 'results'))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, 'checks'), showWarnings = FALSE)
options(stringsAsFactors = FALSE)
data.table::setDTthreads(4L)

# 2. Definitions used by this analysis ----
self_report_labels <- list(BPH = c("1516", "bph / benign prostatic hypertrophy"), depress = c(
  "1286",
  "depression"
), PCa = c("1044", "prostate cancer"))

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "04"
stage_dependencies <- c("01")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "04_self_report_dates.R"), mustWork = TRUE)
stage_code_hash <- tools::md5sum(stage_code_file)
manifest_file <- file.path(results_dir, "Private_Run_Manifest.rds")
manifest <- if (file.exists(manifest_file)) readRDS(manifest_file) else list()
remaining_dependencies <- stage_dependencies
checked_dependencies <- character()
while (length(remaining_dependencies)) {
  dependency <- remaining_dependencies[1]
  remaining_dependencies <- remaining_dependencies[-1]
  if (dependency %in% checked_dependencies) next
  record <- manifest[[dependency]]
  if (is.null(record)) stop("Run upstream part ", dependency, " with this version before part ", stage_name)
  recorded_hashes <- c(record$code, record$inputs, record$outputs)
  if (!length(recorded_hashes) || any(file.exists(as.character(record$absent_inputs))) ||
      any(!file.exists(names(recorded_hashes))) ||
      !identical(unname(tools::md5sum(names(recorded_hashes))), unname(recorded_hashes))) {
    stop("Stale or changed files for upstream part ", dependency, "; regenerate that part and its dependents")
  }
  remaining_dependencies <- c(remaining_dependencies, record$dependencies)
  checked_dependencies <- c(checked_dependencies, dependency)
}
stage_absent_input_files <- character()
stage_input_files <- c(file.path(results_dir, c("01_Baseline_First_GP.rds")),
                         file.path(source_dir, c("self-non-code_participant.csv", "self-non-date-1_participant.csv", "self-can-1_participant.csv")))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

cohort <- readRDS(file.path(results_dir, "01_Baseline_First_GP.rds"))

out <- data.table(Participant_ID = cohort$Participant_ID)
for (condition in c("BPH", "depress", "PCa")) {
  set(out, j = paste0("self_", condition, "_baseline_reported"), value = FALSE)
  set(out, j = paste0("self_", condition, "_first_reported_instance"), value = NA_integer_)
}


qc <- list()

for (disease in c("BPH", "depress", "PCa")) {
  set(out, j = paste0("self_", disease, "_date"), value = as.IDate(rep(NA_character_, nrow(out))))
  set(out, j = paste0("self_", disease, "_unknown"), value = FALSE)
  set(out, j = paste0("self_", disease, "_reported"), value = FALSE)
}

for (domain in c("non_cancer", "cancer")) {
  if (domain == "non_cancer") {
    code <- fread(file.path(source_dir, "self-non-code_participant.csv"),
      colClasses = "character",
      na.strings = c("", "NA")
    )
    year <- fread(file.path(source_dir, "self-non-date-1_participant.csv"),
      colClasses = "character",
      na.strings = c("", "NA")
    )
    code_prefix <- "Non-cancer illness code, self-reported"
    year_prefix <- "Interpolated Year when non-cancer illness first diagnosed"
    diseases <- c("BPH", "depress")
  } else {
    code <- fread(file.path(source_dir, "self-can-1_participant.csv"),
      colClasses = "character",
      na.strings = c("", "NA")
    )
    year <- copy(code)
    code_prefix <- "Cancer code, self-reported"
    year_prefix <- "Interpolated Year when cancer first diagnosed"
    diseases <- "PCa"
  }
  setnames(code, "Participant ID", "Participant_ID")
  setnames(year, "Participant ID", "Participant_ID")
  stopifnot(!anyDuplicated(code$Participant_ID), !anyDuplicated(year$Participant_ID))
  code <- code[match(cohort$Participant_ID, code$Participant_ID)]
  year <- year[match(cohort$Participant_ID, year$Participant_ID)]
  stopifnot(identical(code$Participant_ID, cohort$Participant_ID), identical(
    year$Participant_ID,
    cohort$Participant_ID
  ))
  code_columns <- names(code)[startsWith(names(code), paste0(code_prefix, " | Instance "))]
  for (disease in diseases) {
    earliest <- rep(Inf, nrow(code))
    unknown <- rep(FALSE, nrow(code))
    reported <- rep(FALSE, nrow(code))
    baseline_reported <- rep(FALSE, nrow(code))
    first_instance <- rep(Inf, nrow(code))
    for (code_column in code_columns) {
      # Preserve the complete Instance and Array suffix when pairing the two exports.
      suffix <- substring(code_column, nchar(code_prefix) + 1L)
      year_column <- paste0(year_prefix, suffix)
      stopifnot(year_column %in% names(year))
      match_disease <- tolower(trimws(code[[code_column]])) %in% self_report_labels[[disease]]
      rows <- which(match_disease)
      if (length(rows) == 0L)
        next
      reported[rows] <- TRUE
      instance_text <- sub(".*Instance ", "", code_column)
      instance <- as.integer(strsplit(instance_text, " | ", fixed = TRUE)[[1]][1])
      stopifnot(!is.na(instance))
      if (instance == 0L) baseline_reported[rows] <- TRUE
      first_instance[rows] <- pmin(first_instance[rows], instance)
      raw_year <- year[[year_column]][rows]
      diagnosis_year <- suppressWarnings(as.numeric(raw_year))
      # A positive target report with an empty or unparseable year remains unresolved.
      invalid_year <- is.na(diagnosis_year) | diagnosis_year < 1900 | diagnosis_year >= 2100
      unknown[rows[invalid_year]] <- TRUE
      unknown[rows[diagnosis_year %in% c(-1, -3)]] <- TRUE
      # Unknown/declined years remain flags and are not candidates for the minimum.
      valid <- !is.na(diagnosis_year) & diagnosis_year >= 1900 & diagnosis_year < 2100
      valid_rows <- rows[valid]
      if (length(valid_rows) == 0L)
        next
      y <- diagnosis_year[valid]
      integer_year <- floor(y)
      first_day <- as.IDate(paste0(integer_year, "-01-01"))
      next_year <- as.IDate(paste0(integer_year + 1L, "-01-01"))
      days_in_year <- as.numeric(next_year - first_day)
      approximate_date <- as.numeric(first_day) + round((y - integer_year) * days_in_year)
      earliest[valid_rows] <- pmin(earliest[valid_rows], approximate_date, na.rm = TRUE)
    }
    earliest[!is.finite(earliest)] <- NA_real_
    set(out, j = paste0("self_", disease, "_date"), value = as.IDate(earliest, origin = "1970-01-01"))
    set(out, j = paste0("self_", disease, "_unknown"), value = unknown)
    set(out, j = paste0("self_", disease, "_reported"), value = reported)
    set(out, j = paste0("self_", disease, "_baseline_reported"), value = baseline_reported)
    first_instance[!is.finite(first_instance)] <- NA_real_
    set(out, j = paste0("self_", disease, "_first_reported_instance"), value = as.integer(first_instance))
    set(out, j = paste0("self_", disease, "_reported_without_valid_date"),
        value = reported & is.na(earliest))
    set(out, j = paste0("self_", disease, "_unknown_or_invalid_year"), value = unknown)

    qc[[disease]] <- data.table(
      disease = disease, reported_n = sum(reported), valid_date_n = sum(!is.na(earliest)),
      unknown_year_flag_n = sum(unknown), reported_without_date_n = sum(reported & is.na(earliest))
    )
  }
}

saveRDS(out, file.path(results_dir, "04_Self_Report_Dates.rds"))

fwrite(rbindlist(qc), file.path(results_dir, "checks", "04_Self_Report_Counts.csv"))

cat("Step 04 complete. Self-report sources correctly paired.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^04_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^04_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
