# Prostate cancer dates from the cancer registry.
# Match the diagnosis code and date within the same Instance.

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
# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "03"
stage_dependencies <- c("01")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "03_cancer_registry_dates.R"), mustWork = TRUE)
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
                         file.path(source_dir, c("cancer.csv")))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

cohort <- readRDS(file.path(results_dir, "01_Baseline_First_GP.rds"))

cancer <- fread(file.path(source_dir, "cancer.csv"), colClasses = "character", na.strings = c(
  "",
  "NA"
))

setnames(cancer, "Participant ID", "Participant_ID")

stopifnot(!anyDuplicated(cancer$Participant_ID))

cancer <- cancer[match(cohort$Participant_ID, cancer$Participant_ID)]

stopifnot(identical(cancer$Participant_ID, cohort$Participant_ID))

registry_date <- rep(Inf, nrow(cancer))

matched_instances <- 0L
matched_participant <- rep(FALSE, nrow(cancer))

# Match each code to its own Instance date before calculating any minimum.
for (system in c("10", "9")) {
  code_columns <- grep(paste0("^Type of cancer: ICD", system, " \\| Instance "), names(cancer),
    value = TRUE
  )
  for (code_column in code_columns) {
    instance <- sub(".*Instance ", "", code_column)
    date_column <- paste0("Date of cancer diagnosis | Instance ", instance)
    stopifnot(date_column %in% names(cancer))
    code <- gsub(".", "", sub("\\s.*$", "", trimws(cancer[[code_column]])), fixed = TRUE)
    match_pca <- if (system == "10")
      substr(code, 1, 3) == "C61"
    else code %in% c("185", "1859")
    rows <- which(!is.na(match_pca) & match_pca)
    matched_instances <- matched_instances + length(rows)
    matched_participant[rows] <- TRUE
    current_date <- as.numeric(as.IDate(cancer[[date_column]][rows]))
    registry_date[rows] <- pmin(registry_date[rows], current_date, na.rm = TRUE)
  }
}

registry_date[!is.finite(registry_date)] <- NA_real_

out <- data.table(Participant_ID = cohort$Participant_ID, registry_PCa = as.IDate(registry_date,
  origin = "1970-01-01"
))

out[, registry_PCa_matched := matched_participant]
out[, registry_PCa_unknown := matched_participant & is.na(registry_PCa)]

saveRDS(out, file.path(results_dir, "03_Cancer_Registry_Dates.rds"))

fwrite(
  data.table(matched_instances = matched_instances, matched_participants = sum(matched_participant),
    dated_participants = sum(!is.na(registry_date)), matched_without_valid_date = sum(out$registry_PCa_unknown)),
  file.path(results_dir, "checks", "03_Cancer_Counts.csv")
)

cat("Step 03 complete. PCa registry dates:", sum(!is.na(registry_date)), "\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^03_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^03_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
