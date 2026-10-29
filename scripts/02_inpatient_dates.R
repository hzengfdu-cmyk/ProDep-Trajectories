# Hospital diagnosis dates matched to their original array positions.
# Depression includes ICD10 F32/F33 and the declared ICD9 depression definition.

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
inpatient_icd9 <- list(BPH = "6009", depress = "3119", PCa = "1859")

inpatient_icd10 <- list(BPH = "N40", depress = c("F32", "F33"), PCa = "C61")

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "02"
stage_dependencies <- c("01")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "02_inpatient_dates.R"), mustWork = TRUE)
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
                         file.path(source_dir, c("outcomes_impatient+death+self.csv")))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

cohort <- readRDS(file.path(results_dir, "01_Baseline_First_GP.rds"))

file <- file.path(source_dir, "outcomes_impatient+death+self.csv")

header <- names(fread(file, nrows = 0L))

columns <- c("Participant ID", "Diagnoses - ICD10", "Diagnoses - ICD9", grep("^Date of first in-patient diagnosis - ICD(10|9) \\| Array ",
  header,
  value = TRUE
), "Date of death | Instance 0", "Date of death | Instance 1")

hospital <- fread(file, select = columns, colClasses = "character", na.strings = c("", "NA"))

setnames(hospital, "Participant ID", "Participant_ID")

stopifnot(!anyDuplicated(hospital$Participant_ID))

hospital <- hospital[match(cohort$Participant_ID, hospital$Participant_ID)]

stopifnot(identical(hospital$Participant_ID, cohort$Participant_ID))

out <- data.table(Participant_ID = cohort$Participant_ID)

qc <- list()

for (system in c("10", "9")) {
  diagnosis_column <- paste0("Diagnoses - ICD", system)
  # Keep the original diagnosis positions so missing slots cannot shift Array matches.
  pieces <- strsplit(hospital[[diagnosis_column]], "|", fixed = TRUE)
  codes <- lapply(pieces, function(x) gsub(".", "", sub("\\s.*$", "", trimws(x)), fixed = TRUE))
  for (disease in c("BPH", "depress", "PCa")) {
    target <- if (system == "10")
      inpatient_icd10[[disease]]
    else inpatient_icd9[[disease]]
    if (system == "10") {
      positions <- lapply(codes, function(x) which(substr(x, 1, 3) %in% target))
    } else {
      positions <- lapply(codes, function(x) which(x %in% target))
    }
    hits <- data.table(row_number = rep(seq_len(nrow(hospital)), lengths(positions)), array_number = unlist(positions,
      use.names = FALSE
    ) - 1L)
    date_value <- rep(NA_real_, nrow(hospital))
    if (nrow(hits) > 0L) {
      hits[, `:=`(diagnosis_date, as.IDate(NA_character_))]
      for (array in unique(hits$array_number)) {
        date_column <- paste0(
          "Date of first in-patient diagnosis - ICD", system, " | Array ",
          array
        )
        stopifnot(date_column %in% names(hospital))
        rows <- hits[array_number == array, row_number]
        hits[array_number == array, `:=`(diagnosis_date, as.IDate(hospital[[date_column]][rows]))]
      }
      # Multiple matching diagnoses can have different dates: take the minimum valid date.
      earliest <- hits[!is.na(diagnosis_date), .(date = min(diagnosis_date)), by = row_number]
      date_value[earliest$row_number] <- as.numeric(earliest$date)
    }
    column <- paste0("inpatient_", disease, "_ICD", system)
    set(out, j = column, value = as.IDate(date_value, origin = "1970-01-01"))
    matched <- lengths(positions) > 0L
    unknown <- matched & is.na(date_value)
    set(out, j = paste0(column, "_matched"), value = matched)
    set(out, j = paste0(column, "_unknown"), value = unknown)
    qc[[column]] <- data.table(source = column, matched_records = nrow(hits),
      matched_participants = sum(matched), dated_participants = sum(!is.na(date_value)),
      matched_without_valid_date = sum(unknown))
  }
}

d0 <- as.numeric(as.IDate(hospital[["Date of death | Instance 0"]]))

d1 <- as.numeric(as.IDate(hospital[["Date of death | Instance 1"]]))

early_death <- pmin(d0, d1, na.rm = TRUE)

early_death[!is.finite(early_death)] <- NA_real_

out[, `:=`(inpatient_death, as.IDate(early_death, origin = "1970-01-01"))]

saveRDS(out, file.path(results_dir, "02_Inpatient_Dates.rds"))

fwrite(rbindlist(qc), file.path(results_dir, "checks", "02_Inpatient_Counts.csv"))

cat("Step 02 complete. Aligned male rows:", nrow(out), "\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^02_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^02_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
