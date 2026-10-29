# State occupancy probabilities by baseline age and risk-score strata.
# Prediction starts in BPH at baseline time zero; this is not an individual diagnosis-time reset.

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
STATE_NAMES <- c("Healthy", "BPH", "Depression", "PCa", "Death")

TRANSITIONS <- data.frame(trans = 1:10, from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4), to = c(
  2, 3,
  4, 5, 3, 4, 5, 4, 5, 5
))

TIME_COLUMNS <- c("BPH_time", "depress_time", "PCa_time", "death_time", "censor_time")

# 3. Reused algorithms needed repeatedly in this analysis ----
# Algorithm: build_paths ----
build_paths <- function(d, start_time = 0) {
  stopifnot(!anyNA(d$Participant_ID), !anyDuplicated(d$Participant_ID), all(nzchar(as.character(d$Participant_ID))))
  stopifnot(all(c(TIME_COLUMNS, "Age_base") %in% names(d)))
  for (v in TIME_COLUMNS) if (is.logical(d[[v]]) && all(is.na(d[[v]])))
    d[[v]] <- as.numeric(d[[v]])
  for (v in c(TIME_COLUMNS, "Age_base")) if (!is.numeric(d[[v]]))
    stop("Times and Age_base must be numeric: ", v)
  out <- list()
  off <- list()
  elig <- list()
  ev <- list()
  k <- 0L
  g <- 0L
  h <- 0L
  for (i in seq_len(nrow(d))) {
    z <- as.numeric(d[i, TIME_COLUMNS])
    names(z) <- TIME_COLUMNS
    if (!is.finite(z["censor_time"]))
      stop("Finite censor_time is required to construct risk sets")
    stop_at <- if (!is.na(z["death_time"]) && z["death_time"] <= z["censor_time"])
      z["death_time"]
    else z["censor_time"]
    raw <- z[1:4]
    reason <- ""
    if (!is.finite(stop_at) || stop_at <= start_time)
      reason <- "invalid_followup_end"
    if (any(!is.na(raw) & !is.finite(raw)))
      reason <- "nonfinite_time"
    if (any(!is.na(raw) & raw <= 0))
      reason <- "baseline_or_negative_event"
    if (!is.na(z["death_time"]) && any(raw[1:3] > z["death_time"], na.rm = TRUE))
      reason <- "diagnosis_after_death"
    if (!is.na(stop_at) && any(raw[1:3] == stop_at, na.rm = TRUE))
      reason <- "diagnosis_at_endpoint"
    use <- which(!is.na(raw) & is.finite(raw) & raw <= stop_at)
    if (anyDuplicated(raw[use]))
      reason <- "same_day_unresolved"
    if (start_time > 0 && any(raw[use] <= start_time))
      reason <- "event_before_landmark"
    elig[[i]] <- data.frame(Participant_ID = d$Participant_ID[i], eligible = reason ==
      "", reason = reason)
    if (nzchar(reason))
      next
    current <- 1L
    entry <- start_time
    for (j in use[order(raw[use])]) {
      dest <- j + 1L
      tt <- TRANSITIONS$trans[TRANSITIONS$from == current & TRANSITIONS$to == dest]
      if (!length(tt)) {
        g <- g + 1L
        off[[g]] <- data.frame(
          Participant_ID = d$Participant_ID[i], state = current,
          diagnosis = STATE_NAMES[dest], time = raw[j]
        )
        next
      }
      k <- k + 1L
      out[[k]] <- data.frame(
        Participant_ID = d$Participant_ID[i], from = current, Tstart = entry,
        Tstop = raw[j], actual_trans = tt, Age_base = d$Age_base[i]
      )
      h <- h + 1L
      ev[[h]] <- data.frame(Participant_ID = d$Participant_ID[i], trans = tt, time = raw[j])
      current <- dest
      entry <- raw[j]
      if (current == 5L)
        break
    }
    if (current != 5L && stop_at > entry) {
      k <- k + 1L
      out[[k]] <- data.frame(
        Participant_ID = d$Participant_ID[i], from = current, Tstart = entry,
        Tstop = stop_at, actual_trans = NA_integer_, Age_base = d$Age_base[i]
      )
    }
  }
  list(
    episodes = as.data.frame(data.table::rbindlist(out, use.names = TRUE, fill = TRUE)),
    events = as.data.frame(data.table::rbindlist(ev, use.names = TRUE, fill = TRUE)), off_graph = as.data.frame(data.table::rbindlist(off,
      use.names = TRUE, fill = TRUE
    )), eligibility = as.data.frame(data.table::rbindlist(elig,
      use.names = TRUE, fill = TRUE
    ))
  )
}

# Algorithm: expand_risk_sets ----
expand_risk_sets <- function(episodes, covariates = NULL) {
  stopifnot(all(c("Participant_ID", "from", "Tstart", "Tstop", "actual_trans", "Age_base") %in%
    names(episodes)))
  if (any(!is.finite(episodes$Tstart) | !is.finite(episodes$Tstop) | episodes$Tstop <= episodes$Tstart))
    stop("Invalid or zero risk interval")
  rows <- lapply(seq_len(nrow(episodes)), function(i) {
    e <- episodes[i, , drop = FALSE]
    trans <- TRANSITIONS[TRANSITIONS$from == e$from, , drop = FALSE]
    z <- e[rep(1, nrow(trans)), , drop = FALSE]
    z$trans <- trans$trans
    z$to <- trans$to
    z$status <- as.integer(!is.na(e$actual_trans) & trans$trans == e$actual_trans)
    z
  })
  long <- as.data.frame(data.table::rbindlist(rows, use.names = TRUE, fill = TRUE))
  if (!is.null(covariates)) {
    stopifnot(
      !anyNA(covariates$Participant_ID), !anyDuplicated(covariates$Participant_ID),
      all(nzchar(as.character(covariates$Participant_ID)))
    )
    idx <- match(long$Participant_ID, covariates$Participant_ID)
    if (anyNA(idx))
      stop("Covariates do not cover risk set IDs")
    vars <- setdiff(names(covariates), names(long))
    long <- cbind(long, covariates[idx, vars, drop = FALSE])
  }
  long$time <- long$Tstop - long$Tstart
  long$Tstart_age <- long$Age_base + long$Tstart
  long$Tstop_age <- long$Age_base + long$Tstop
  long$age5 <- floor(long$Tstart_age / 5) * 5
  attr(long, "trans") <- transition_matrix()
  class(long) <- c("msdata", "data.frame")
  long
}

# Algorithm: checked_coxph ----
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

# Algorithm: transition_matrix ----
transition_matrix <- function() {
  m <- matrix(NA_integer_, 5, 5, dimnames = list(STATE_NAMES, STATE_NAMES))
  m[cbind(TRANSITIONS$from, TRANSITIONS$to)] <- TRANSITIONS$trans
  m
}

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "10"
stage_dependencies <- c("06")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "10_state_occupancy.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("06_Analysis_Data.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

data <- readRDS(file.path(results_dir, "06_Analysis_Data.rds"))

young <- data$Age_base < 55

low_score <- data$lifestyle < 3

subgroups <- list(age_lt55 = young, age_ge55 = !young, young_low = young & low_score, older_low = !young &
  low_score, young_high = young & !low_score, older_high = !young & !low_score)

probabilities <- list()

for (stratum in names(subgroups)) {
  subgroup <- data[subgroups[[stratum]], , drop = FALSE]
  paths <- build_paths(subgroup)
  long <- expand_risk_sets(paths$episodes)
  fit <- checked_coxph(Surv(Tstart, Tstop, status) ~ strata(trans), long, ties = "breslow")
  hazards <- mstate::msfit(fit, trans = transition_matrix(), variance = FALSE)
  probability <- mstate::probtrans(hazards, predt = 0, variance = FALSE)[[2]]
  probability$stratum <- stratum
  probability$initial_state <- "BPH"
  probabilities[[stratum]] <- probability
}

fwrite(rbindlist(probabilities), file.path(results_dir, "10_BPH_Initial_State_Probabilities.csv"))

cat("Step 10 complete. Time origin remains baseline, not individual BPH diagnosis.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^10_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^10_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
