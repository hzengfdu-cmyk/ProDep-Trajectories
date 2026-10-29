# Earliest diagnosis dates, follow-up endpoints and cohort exclusions.
# Retain source dates, apply censoring, and record the exclusion reason for each participant.

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
exclude_self_report_unknown <- TRUE

exclude_diagnosis_at_censor <- TRUE

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "05"
stage_dependencies <- c("01", "02", "03", "04")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "05_followup_and_exclusions.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("01_Baseline_First_GP.rds", "02_Inpatient_Dates.rds", "03_Cancer_Registry_Dates.rds", "04_Self_Report_Dates.rds")),
                         file.path(source_dir, c("lost_follow_up.csv")))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

cohort <- readRDS(file.path(results_dir, "01_Baseline_First_GP.rds"))

for (file in c("02_Inpatient_Dates.rds", "03_Cancer_Registry_Dates.rds", "04_Self_Report_Dates.rds")) {
  current <- readRDS(file.path(results_dir, file))
  stopifnot(!anyDuplicated(current$Participant_ID))
  cohort <- merge(cohort, current, by = "Participant_ID", all.x = TRUE, sort = FALSE)
}

initial_n <- nrow(cohort)

stopifnot(!anyDuplicated(cohort$Participant_ID))

# Apply the source-specific definitions before taking cross-source minimum dates.
source_columns <- list(BPH = c(
  "first_BPH", "GP_BPH", "inpatient_BPH_ICD10", "inpatient_BPH_ICD9",
  "self_BPH_date"
), depress = c(
  "first_depress", "GP_depress", "inpatient_depress_ICD10",
  "inpatient_depress_ICD9", "self_depress_date"
), PCa = c(
  "GP_PCa", "inpatient_PCa_ICD10", "inpatient_PCa_ICD9", "registry_PCa",
  "self_PCa_date"
))

source_summary <- list()

for (disease in names(source_columns)) {
  columns <- source_columns[[disease]]
  date_matrix <- do.call(cbind, lapply(cohort[, ..columns], as.numeric))
  minimum <- do.call(pmin, c(as.data.frame(date_matrix), na.rm = TRUE))
  minimum[!is.finite(minimum)] <- NA_real_
  set(cohort, j = paste0(disease, "_date_all_sources"), value = as.IDate(minimum, origin = "1970-01-01"))
  set(cohort, j = paste0(disease, "_date"), value = as.IDate(minimum, origin = "1970-01-01"))
  source_summary[[disease]] <- data.table(disease = disease, source = columns, dated_n = colSums(!is.na(date_matrix)))
}

death <- pmin(as.numeric(cohort$first_death), as.numeric(cohort$inpatient_death), na.rm = TRUE)

death[!is.finite(death)] <- NA_real_

cohort[, `:=`(death_date, as.IDate(death, origin = "1970-01-01"))]

cohort[, `:=`(death_source_disagreement, !is.na(first_death) & !is.na(inpatient_death) & first_death !=
  inpatient_death)]

lost <- fread(file.path(source_dir, "lost_follow_up.csv"),
  select = c("Participant ID", "Date lost to follow-up"),
  colClasses = "character", na.strings = c("", "NA")
)

setnames(lost, "Participant ID", "Participant_ID")

stopifnot(!anyDuplicated(lost$Participant_ID))

lost[, `:=`(loss_date, as.IDate(`Date lost to follow-up`))]

cohort <- merge(cohort, lost[, .(Participant_ID, loss_date)],
  by = "Participant_ID", all.x = TRUE,
  sort = FALSE
)

cohort[, `:=`(censor_date, as.IDate(pmin(as.numeric(admin_end_date), as.numeric(loss_date),
  na.rm = TRUE
), origin = "1970-01-01"))]

stopifnot(nrow(cohort) == initial_n, !anyNA(cohort$censor_date))

cohort[, `:=`(exclusion_reason, "")]

# Identify baseline disease before deleting events beyond the observation endpoint.
raw_dates <- c("BPH_date", "depress_date", "PCa_date", "death_date")

raw_time_matrix <- do.call(cbind, lapply(cohort[, ..raw_dates], function(x) as.numeric(x -
  cohort$baseline_date)))

baseline_flag <- rowSums(raw_time_matrix <= 0, na.rm = TRUE) > 0

gp_baseline_columns <- paste0("GP_", c("BPH", "depress", "PCa"), "_baseline")

baseline_flag <- baseline_flag | rowSums(as.matrix(cohort[, ..gp_baseline_columns]), na.rm = TRUE) >
  0

baseline_columns <- c("first_BPH_baseline", "first_depress_baseline",
  paste0("self_", c("BPH", "depress", "PCa"), "_baseline_reported"))
baseline_flag <- baseline_flag | rowSums(as.matrix(cohort[, ..baseline_columns]), na.rm = TRUE) > 0
cohort[, qc_baseline_event := baseline_flag]
cohort[baseline_flag, `:=`(exclusion_reason, "baseline_or_negative_event")]

invalid_followup <- cohort$censor_date <= cohort$baseline_date
cohort[, qc_nonpositive_followup := invalid_followup]

cohort[exclusion_reason == "" & invalid_followup, `:=`(exclusion_reason, "nonpositive_followup")]

# Preserve inpatient/registry evidence independently of the current exclusion policy.
# Author-confirmed: these source-only unknown dates add no exclusion; retain their QC flags.
hospital_registry_flags <- lapply(c("BPH", "depress", "PCa"), function(disease) {
  columns <- paste0("inpatient_", disease, "_ICD", c("10", "9"), "_unknown")
  if (disease == "PCa") columns <- c(columns, "registry_PCa_unknown")
  stopifnot(all(columns %in% names(cohort)))
  unresolved <- rowSums(as.matrix(cohort[, ..columns]), na.rm = TRUE) > 0 &
    is.na(cohort[[paste0(disease, "_date_all_sources")]])
  set(cohort, j = paste0("qc_", disease, "_hospital_registry_unknown"), value = unresolved)
  data.table(condition = disease, unmatched_date_all_sources = sum(unresolved),
    exclusion_policy = "Audit only; no additional exclusion for hospital/registry-only unknown dates")
})
fwrite(rbindlist(hospital_registry_flags), file.path(results_dir, "checks", "05_Hospital_Registry_Unknown.csv"))

# Apply the existing first-occurrence, GP and self-report unknown-date exclusion.
# Hospital/registry-only unknown flags do not trigger this exclusion; no diagnosis day is invented.
unknown_flag <- rep(FALSE, nrow(cohort))

for (disease in c("BPH", "depress", "PCa")) {
  self_unknown <- cohort[[paste0("self_", disease, "_unknown")]]
  gp_unknown <- as.logical(cohort[[paste0("GP_", disease, "_unknown")]])
  no_valid_date <- is.na(cohort[[paste0(disease, "_date_all_sources")]])
  first_unknown <- if (disease == "PCa") rep(FALSE, nrow(cohort)) else cohort[[paste0("first_", disease, "_unknown")]]
  unknown_flag <- unknown_flag | ((self_unknown | gp_unknown | first_unknown) & no_valid_date)
}

cohort[, qc_unknown_target_date := unknown_flag & exclude_self_report_unknown]

if (exclude_self_report_unknown) {
  cohort[exclusion_reason == "" & unknown_flag, `:=`(exclusion_reason, "target_condition_unknown_date")]
}

after_censor <- list()

for (column in raw_dates) {
  too_late <- !is.na(cohort[[column]]) & cohort[[column]] > cohort$censor_date
  after_censor[[column]] <- data.table(event = column, dates_after_censor = sum(too_late))
  set(cohort, i = which(too_late), j = column, value = as.IDate(NA_character_))
}

disease_matrix <- as.matrix(cohort[, .(as.numeric(BPH_date), as.numeric(depress_date), as.numeric(PCa_date))])

all_matrix <- cbind(disease_matrix, as.numeric(cohort$death_date))

same_day <- rep(FALSE, nrow(cohort))

for (a in 1:3) for (b in (a + 1):4) {
  same_day <- same_day | (!is.na(all_matrix[, a]) & !is.na(all_matrix[, b]) & all_matrix[
    ,
    a
  ] == all_matrix[, b])
}

cohort[, qc_same_day := same_day]
cohort[exclusion_reason == "" & same_day, `:=`(exclusion_reason, "same_day_unresolved")]

at_censor <- rowSums(disease_matrix == as.numeric(cohort$censor_date), na.rm = TRUE) > 0

cohort[, qc_at_censor := at_censor & exclude_diagnosis_at_censor]

if (exclude_diagnosis_at_censor) {
  cohort[exclusion_reason == "" & at_censor, `:=`(exclusion_reason, "diagnosis_at_censor")]
}

after_death <- !is.na(cohort$death_date) & rowSums(disease_matrix > as.numeric(cohort$death_date),
  na.rm = TRUE
) > 0

cohort[, qc_after_death := after_death]
cohort[, pre_tie_eligible := !qc_baseline_event & !qc_nonpositive_followup & !qc_unknown_target_date & !qc_after_death]
cohort[exclusion_reason == "" & after_death, `:=`(exclusion_reason, "diagnosis_after_death")]

for (disease in c("BPH", "depress", "PCa", "death")) {
  set(cohort, j = paste0(disease, "_time"), value = as.numeric(cohort[[paste0(disease, "_date")]] -
    cohort$baseline_date) / 365.25)
}

cohort[, `:=`(censor_time, as.numeric(censor_date - baseline_date) / 365.25)]

kept <- cohort[exclusion_reason == ""]

stopifnot(!anyDuplicated(kept$Participant_ID), all(kept$censor_time > 0))

saveRDS(cohort, file.path(results_dir, "05_All_Outcome_Dates.rds"))

saveRDS(kept, file.path(results_dir, "05_Eligible_Outcome_Cohort.rds"))

fwrite(cohort[, .(Participant_ID, exclusion_reason)], file.path(results_dir, "05_Private_Exclusion_Ledger.csv"))

flow <- cohort[, .(n = .N), by = exclusion_reason]

flow[exclusion_reason == "", `:=`(exclusion_reason, "retained")]

fwrite(flow, file.path(results_dir, "checks", "05_Exclusion_Counts.csv"))

fwrite(rbindlist(source_summary), file.path(results_dir, "checks", "05_Source_Date_Counts.csv"))

fwrite(rbindlist(after_censor), file.path(results_dir, "checks", "05_After_Censor_Counts.csv"))

fwrite(
  data.table(
    initial_n = initial_n, retained_n = nrow(kept), loss_before_admin_n = sum(cohort$loss_date <
      cohort$admin_end_date, na.rm = TRUE), death_source_disagreement_n = sum(cohort$death_source_disagreement)
  ),
  file.path(results_dir, "checks", "05_Run_Summary.csv")
)

cat("Step 05 complete. Initial:", initial_n, "retained:", nrow(kept), "\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^05_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^05_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
