# Standalone score categorical Cox analysis.
# All paths and algorithms required by this analysis are defined here.
# Inspect long_data, model_data, exposure_columns, formula, fit and results.

# Packages and local paths ----
library(data.table)

library(survival)

source_dir <- Sys.getenv("UKB_SOURCE_DIR", unset = file.path(getwd(), "data", "raw"))

results_dir <- Sys.getenv("UKB_RESULTS_DIR", unset = file.path(getwd(), "results"))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(results_dir, "checks"), showWarnings = FALSE)

options(stringsAsFactors = FALSE)

data.table::setDTthreads(4L)

checked_coxph <- function(formula, data, ties = "efron", ...) {
  warnings <- character()
  fit <- withCallingHandlers(survival::coxph(formula,
    data = data, ties = ties, na.action = stats::na.fail,
    model = TRUE, x = TRUE, y = TRUE, ...
  ), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  fit$v1_warnings <- warnings
  if (any(grepl("infinite|converg|iterations", warnings, ignore.case = TRUE)))
    stop("Unstable Cox fit: ", paste(warnings, collapse = "; "))
  if (any(!is.finite(stats::coef(fit))))
    stop("Non-estimable Cox coefficients; sparse/constant covariates or singular model")
  fit
}

coefficient_results <- function(fit, map) {
  cf <- stats::coef(fit)
  V <- stats::vcov(fit)
  ans <- map
  idx <- match(ans$term, names(cf))
  ans$beta <- cf[idx]
  ans$SE <- sqrt(diag(V))[idx]
  ans$HR <- exp(ans$beta)
  ans$lower <- exp(ans$beta - 1.96 * ans$SE)
  ans$upper <- exp(ans$beta + 1.96 * ans$SE)
  ans$P_value <- 2 * stats::pnorm(-abs(ans$beta / ans$SE))
  ans
}

# Inspect this saved risk set before building the model.
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "08b"
stage_dependencies <- c("07")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "08b_categorical_score_cox.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("07_Private_Risk_Sets.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

# Read the risk set and prepare the model ----
long_data <- readRDS(file.path(results_dir, "07_Private_Risk_Sets.rds"))

model_data <- as.data.frame(long_data)

{
  score_group <- ifelse(model_data$lifestyle <= 1, "0-1", ifelse(model_data$lifestyle >=
    5, "5-6", as.character(model_data$lifestyle)))
  exposure_columns <- c("score_2", "score_3", "score_4", "score_5_6")
  for (j in seq_along(exposure_columns)) {
    model_data[[exposure_columns[j]]] <- as.integer(score_group == c("2", "3", "4", "5-6")[j])
  }
}

stopifnot(!anyNA(model_data[, c(exposure_columns, "Townsend", "region")]))

model_data$region <- factor(model_data$region, levels = c("England", "Scotland", "Wales"))

stopifnot(!anyNA(model_data$region))

model_data$age5 <- floor(model_data$Tstart_age / 5) * 5

# Each exposure has a separate coefficient for each transition.
model_terms <- character()

coefficient_map <- list()

for (exposure in exposure_columns) {
  value <- as.numeric(as.character(model_data[[exposure]]))
  stopifnot(!anyNA(value))
  for (transition in 1:10) {
    term <- sprintf("%s_T%02d", exposure, transition)
    model_data[[term]] <- value * as.integer(model_data$trans == transition)
    model_terms <- c(model_terms, term)
    coefficient_map[[term]] <- data.frame(term = term, exposure = exposure, trans = transition)
  }
}

# Score models use a common TDI effect and do not adjust for their component factors.
model_terms <- c(model_terms, "Townsend")

# Attained-age intervals, joint strata, participant clustering, and Efron ties.
formula <- as.formula(paste("Surv(Tstart_age, Tstop_age, status) ~", paste(c(
  model_terms, "strata(trans, age5, region)",
  "cluster(Participant_ID)"
), collapse = " + ")))

# Fit the specified model and save results ----
fit <- checked_coxph(formula, model_data, ties = "efron")

# Use the fitted robust covariance for coefficient intervals.
results <- coefficient_results(fit, rbindlist(coefficient_map))

# Event counts use transition events in the exact fitted risk set.
event_counts <- list()
for (exposure in "lifestyle_group") {
  count_data <- as.data.table(model_data)
  count_data[, group := score_group]
  counts <- count_data[, .(risk_n = uniqueN(Participant_ID), events = sum(status)), by = .(trans, group)]
  counts[, exposure := exposure]
  event_counts[[exposure]] <- counts
}
fwrite(rbindlist(event_counts), file.path(results_dir, "08_score_categorical_Event_Counts.csv"))
results$total_events <- vapply(results$trans, function(transition) sum(model_data$status[model_data$trans == transition]), numeric(1))
results$exposed_events <- vapply(seq_len(nrow(results)), function(i) {
  rows <- model_data$trans == results$trans[i] & model_data[[results$exposure[i]]] == 1
  sum(model_data$status[rows])
}, numeric(1))
results$reference_events <- vapply(results$trans, function(transition) sum(model_data$status[model_data$trans == transition & score_group == "0-1"]), numeric(1))

fwrite(results, file.path(results_dir, "08_score_categorical_HR.csv"))

saveRDS(fit, file.path(results_dir, "08_Private_score_categorical_Model.rds"))

cat("Completed score_categorical model.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^08_(score_categorical|Private_score_categorical)", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^08b_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
