# Male baseline, first-occurrence records, and screened GP dates.
# Read the existing GP summary, or run 00_gp_screening.R locally to prepare it.

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
admin_end_dates <- data.table::data.table(region = c("England", "Scotland", "Wales"), admin_end_date = as.Date(c(
  "2023-03-31",
  "2022-08-31", "2022-05-31"
)))

centre_region <- data.table::data.table(centre = c(
  "Barts", "Birmingham", "Bristol", "Bury",
  "Cardiff", "Croydon", "Edinburgh", "Glasgow", "Hounslow", "Leeds", "Liverpool", "Manchester",
  "Middlesbrough", "Middlesborough", "Newcastle", "Nottingham", "Oxford", "Reading", "Sheffield",
  "Stockport (pilot)", "Stoke", "Swansea", "Wrexham"
), region = c(
  "England", "England", "England",
  "England", "Wales", "England", "Scotland", "Scotland", "England", "England", "England",
  "England", "England", "England", "England", "England", "England", "England", "England",
  "England", "England", "Wales", "Wales"
))

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Start from the supplied male cohort; external sources never define the cohort size.
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "01"
stage_dependencies <- character()
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "01_baseline_first_occurrence_GP.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("00_GP_Dates.csv")),
                         file.path(source_dir, c("base date + censor region.csv", "first+dea_participant.csv")))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

baseline <- fread(file.path(source_dir, "base date + censor region.csv"),
  colClasses = "character",
  na.strings = c("", "NA")
)

setnames(baseline, names(baseline)[1], "Participant_ID")

stopifnot(!anyNA(baseline$Participant_ID), !anyDuplicated(baseline$Participant_ID))

baseline[, `:=`(baseline_date, as.IDate(`Date of consenting to join UK Biobank`))]

stopifnot(!anyNA(baseline$baseline_date))

baseline[, `:=`(centre, `UK Biobank assessment centre | Instance 0`)]

baseline[, `:=`(region, centre_region$region[match(centre, centre_region$centre)])]

stopifnot(!anyNA(baseline$region))

baseline[, `:=`(admin_end_date, admin_end_dates$admin_end_date[match(region, admin_end_dates$region)])]

baseline <- baseline[, .(Participant_ID, baseline_date, centre, region, admin_end_date)]

# First-occurrence depression uses F32 only; F33 is not read from this source.
first <- fread(file.path(source_dir, "first+dea_participant.csv"),
  select = c(
    "Participant ID",
    "Date N40 first reported (hyperplasia of prostate)", "Date F32 first reported (depressive episode)",
    "Date of death | Instance 0", "Date of death | Instance 1"
  ), colClasses = "character",
  na.strings = c("", "NA")
)

setnames(first, names(first)[1], "Participant_ID")

stopifnot(!anyDuplicated(first$Participant_ID))

first <- first[match(baseline$Participant_ID, first$Participant_ID)]

stopifnot(identical(first$Participant_ID, baseline$Participant_ID))

# UKB Coding 819 applies to first-occurrence fields only.
first_columns <- c("Date N40 first reported (hyperplasia of prostate)",
  "Date F32 first reported (depressive episode)")
first_conditions <- c("BPH", "depress")
for (condition in c("BPH", "depress")) {
  set(first, j = paste0("first_", condition), value = as.IDate(rep(NA_character_, nrow(first))))
  set(first, j = paste0("first_", condition, "_baseline"), value = FALSE)
  set(first, j = paste0("first_", condition, "_unknown"), value = FALSE)
}
date_quality <- list()
for (i in seq_along(first_columns)) {
  raw_date <- first[[first_columns[i]]]
  parsed_date <- as.IDate(raw_date)
  baseline_evidence <- raw_date %in% c("1902-02-02", "1903-03-03")
  unknown_evidence <- raw_date %in% c("1900-01-01", "1901-01-01", "1909-09-09", "2037-07-07") |
    (!is.na(raw_date) & is.na(parsed_date))
  parsed_date[baseline_evidence | unknown_evidence] <- as.IDate(NA_character_)
  condition <- first_conditions[i]
  old_date <- as.numeric(first[[paste0("first_", condition)]])
  minimum <- pmin(old_date, as.numeric(parsed_date), na.rm = TRUE)
  minimum[!is.finite(minimum)] <- NA_real_
  set(first, j = paste0("first_", condition), value = as.IDate(minimum, origin = "1970-01-01"))
  set(first, j = paste0("first_", condition, "_baseline"),
      value = first[[paste0("first_", condition, "_baseline")]] | baseline_evidence)
  set(first, j = paste0("first_", condition, "_unknown"),
      value = first[[paste0("first_", condition, "_unknown")]] | unknown_evidence)
  date_quality[[i]] <- data.table(field = first_columns[i], baseline_mask = sum(baseline_evidence),
    unknown_or_invalid = sum(unknown_evidence), valid_dates = sum(!is.na(parsed_date)))
}
fwrite(rbindlist(date_quality), file.path(results_dir, "checks", "01_First_Occurrence_Date_Quality.csv"))

d0 <- as.numeric(as.IDate(first[["Date of death | Instance 0"]]))

d1 <- as.numeric(as.IDate(first[["Date of death | Instance 1"]]))

death_date <- pmin(d0, d1, na.rm = TRUE)

death_date[!is.finite(death_date)] <- NA_real_

first[, `:=`(first_death, as.IDate(death_date, origin = "1970-01-01"))]

first <- first[, .(Participant_ID, first_BPH, first_depress, first_death,
  first_BPH_baseline, first_BPH_unknown, first_depress_baseline, first_depress_unknown)]

# Read the GP dates and date-quality flags prepared in part 00.
gp <- fread(file.path(results_dir, "00_GP_Dates.csv"), colClasses = "character", na.strings = c(
  "",
  "NA"
))

setnames(gp, c("Participant ID", "BPH-date", "depress-date", "PCa-date"), c(
  "Participant_ID",
  "GP_BPH", "GP_depress", "GP_PCa"
))

stopifnot(!anyDuplicated(gp$Participant_ID))

for (column in c("GP_BPH", "GP_depress", "GP_PCa")) {
  old_value <- gp[[column]]
  gp[[column]] <- as.IDate(old_value)
  stopifnot(!any(!is.na(old_value) & is.na(gp[[column]])))
}

cohort <- merge(baseline, first, by = "Participant_ID", all.x = TRUE, sort = FALSE)

cohort <- merge(cohort, gp, by = "Participant_ID", all.x = TRUE, sort = FALSE)

stopifnot(nrow(cohort) == nrow(baseline), !anyDuplicated(cohort$Participant_ID))



fwrite(data.table(check = c("male_baseline_rows", "first_source_rows", "GP_rows"), value = c(
  nrow(baseline),
  nrow(first), nrow(gp)
)), file.path(results_dir, "checks", "01_Source_Counts.csv"))



for (column in grep("^GP_.*_(baseline|unknown)$", names(cohort), value = TRUE)) {
  flag <- as.integer(cohort[[column]])
  flag[is.na(flag)] <- 0L
  set(cohort, j = column, value = flag)
}
saveRDS(cohort, file.path(results_dir, "01_Baseline_First_GP.rds"))

cat("Step 01 complete. Male baseline:", nrow(cohort), "\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^01_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^01_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
