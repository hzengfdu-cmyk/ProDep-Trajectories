# Event-timing sensitivity analyses.
# Use scored pre-tie records, not an attempt to recreate excluded patients from the final cohort.

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

RISK_FACTORS <- c(
  "smoke", "sport", "sleep", "social_support_score", "new_food_score_score",
  "sexpartners_score"
)

TIME_COLUMNS <- c("BPH_time", "depress_time", "PCa_time", "death_time", "censor_time")

# 3. Reused algorithms needed repeatedly in this analysis ----
# Algorithm: run_timing_sensitivity ----
run_timing_sensitivity <- function(d, pre_tie_data = NULL) {
  reg <- timing_registry()
  outputs <- list()
  skips <- list()
  eligibility_audit <- list()
  risk_audit <- list()
  for (i in seq_len(nrow(reg))) {
    s <- reg[i, ]
    z <- d
    start <- 0
    only <- 1:10
    if (s$type == "same_day") {
      if (is.null(pre_tie_data)) {
        skips[[s$spec_id]] <- data.frame(spec_id = s$spec_id, status = "NOT_RUN", message = "Pre-tie cohort required; excluded cases cannot be reconstructed from final cohort")
        next
      }
      hard_flags <- c("qc_baseline_event", "qc_nonpositive_followup", "qc_unknown_target_date", "qc_after_death")
      stopifnot(all(hard_flags %in% names(pre_tie_data)))
      stopifnot(!anyNA(pre_tie_data[, hard_flags, drop = FALSE]),
        !any(as.matrix(pre_tie_data[, hard_flags, drop = FALSE])))
      z <- order_same_day(pre_tie_data, if (s$spec_id == "SD_UPSTREAM_FIRST")
        "upstream_first"
      else "endpoint_priority")
      filtered <- filter_cohort(z)
      filter_counts <- as.data.table(filtered$exclusions)[, .(n = .N), by = reason]
      filter_counts[, `:=`(spec_id = s$spec_id, eligible = reason == "retained")]
      eligibility_audit[[paste0(s$spec_id, "_filter")]] <- filter_counts
      z <- filtered$data
    }
    if (s$type == "landmark")
      start <- s$lag_months / 12
    if (s$type == "reclassification") {
      short <- !is.na(z$depress_time) & !is.na(z$PCa_time) & z$PCa_time >= z$depress_time &
        z$PCa_time - z$depress_time < s$lag_months / 12
      z$depress_time[short] <- NA_real_
    }
    p <- build_paths(z, start_time = start)
    eligibility <- as.data.table(p$eligibility)
    eligibility[, spec_id := s$spec_id]
    eligibility_audit[[s$spec_id]] <- eligibility[, .(n = .N), by = .(spec_id, eligible, reason)]

    if (s$type == "delayed_entry") {
      long <- depression_delayed_risk(p$episodes, z, s$lag_months / 12)
      only <- 8L
    } else long <- expand_risk_sets(p$episodes, z)
    if (!nrow(long)) {
      skips[[s$spec_id]] <- data.frame(
        spec_id = s$spec_id, status = "NOT_ESTIMABLE",
        message = "No eligible risk intervals"
      )
      next
    }
    outputs[[s$spec_id]] <- fit_timing_models(long, s$spec_id, only_trans = only)
    outputs[[s$spec_id]]$specification_type <- s$type
    outputs[[s$spec_id]]$lag_months <- s$lag_months
    risk_table <- as.data.table(long)[, .(risk_n = uniqueN(Participant_ID),
      events = sum(status), person_years = sum(Tstop - Tstart)), by = trans]
    risk_table[, `:=`(spec_id = s$spec_id, lag_months = s$lag_months,
                     time_origin = "baseline; delayed entry retains original time")]
    risk_audit[[s$spec_id]] <- risk_table
  }
  list(registry = reg, results = as.data.frame(data.table::rbindlist(outputs,
    use.names = TRUE,
    fill = TRUE
  )), eligibility = as.data.frame(data.table::rbindlist(eligibility_audit, fill = TRUE)),
    transition_risk = as.data.frame(data.table::rbindlist(risk_audit, fill = TRUE)),
    not_run = as.data.frame(data.table::rbindlist(skips,
    use.names = TRUE,
    fill = TRUE
  )))
}

# Algorithm: timing_registry ----
timing_registry <- function() {
  rbind(data.frame(spec_id = c(
    "BASE", "SD_UPSTREAM_FIRST", "SD_ENDPOINT_PRIORITY", "LANDMARK_1Y",
    "LANDMARK_2Y"
  ), type = c("base", "same_day", "same_day", "landmark", "landmark"), lag_months = c(
    NA,
    NA, NA, 12, 24
  )), data.frame(
    spec_id = paste0("DP_RECLASS_", c(3, 6, 12, 24), "M"),
    type = "reclassification", lag_months = c(3, 6, 12, 24)
  ), data.frame(spec_id = paste0(
    "DP_DELAYED_",
    c(3, 6, 12, 24), "M"
  ), type = "delayed_entry", lag_months = c(3, 6, 12, 24)))
}

# Algorithm: order_same_day ----
order_same_day <- function(d, rule = c("upstream_first", "endpoint_priority")) {
  rule <- match.arg(rule)
  z <- d
  for (i in seq_len(nrow(z))) {
    times <- as.numeric(z[i, TIME_COLUMNS[1:4]])
    groups <- split(which(!is.na(times)), times[!is.na(times)])
    for (g in groups) if (length(g) > 1) {
      t <- times[g[1]]
      if (t <= 0)
        next
      if (rule == "endpoint_priority") {
        winner <- max(g)
        times[setdiff(g, winner)] <- NA_real_
      } else {
        prior <- max(c(0, times[!is.na(times) & times < t]))
        delta <- min(0.25 / 365.25, (t - prior) / (length(g) + 1))
        g <- sort(g)
        times[g] <- t - rev(seq.int(0, length(g) - 1)) * delta
      }
    }
    boundary <- min(c(z$censor_time[i], times[4]), na.rm = TRUE)
    at_boundary <- which(!is.na(times[1:3]) & times[1:3] == boundary)
    if (length(at_boundary)) {
      if (rule == "endpoint_priority")
        times[at_boundary] <- NA_real_
      else {
        prior <- max(c(0, times[!is.na(times) & times < boundary]))
        delta <- min(0.25 / 365.25, (boundary - prior) / (length(at_boundary) + 1))
        times[at_boundary] <- boundary - rev(seq_along(at_boundary)) * delta
      }
    }
    z[i, TIME_COLUMNS[1:4]] <- as.list(times)
  }
  z
}

# Algorithm: filter_cohort ----
filter_cohort <- function(d) {
  stopifnot(!anyNA(d$Participant_ID), !anyDuplicated(d$Participant_ID), all(nzchar(as.character(d$Participant_ID))))
  stopifnot(all(TIME_COLUMNS %in% names(d)))
  for (v in TIME_COLUMNS) if (is.logical(d[[v]]) && all(is.na(d[[v]])))
    d[[v]] <- as.numeric(d[[v]])
  for (v in TIME_COLUMNS) if (!is.numeric(d[[v]]))
    stop("Time columns must be numeric years: ", v)
  tm <- as.matrix(d[, TIME_COLUMNS[1:4], drop = FALSE])
  diag <- tm[, 1:3, drop = FALSE]
  reason <- rep("", nrow(d))
  mark <- function(flag, label) {
    idx <- which(flag & reason == "")
    reason[idx] <<- label
  }
  mark(apply(tm, 1, function(x) any(x <= 0, na.rm = TRUE)), "baseline_or_negative_time")
  # Read immutable source QC; scenario-removed dates are not missing source dates.
  hard_flags <- c("qc_baseline_event", "qc_nonpositive_followup", "qc_unknown_target_date", "qc_after_death")
  stopifnot(all(hard_flags %in% names(d)), !anyNA(d[, hard_flags, drop = FALSE]))
  for (flag in hard_flags) mark(d[[flag]], paste0("source_", flag))
  mark(apply(tm, 1, function(x) anyDuplicated(x[!is.na(x)]) > 0), "same_day_unresolved")
  mark(rowSums(diag == d$censor_time, na.rm = TRUE) > 0, "diagnosis_at_censor")
  mark(rowSums(diag > d$death_time, na.rm = TRUE) > 0 & !is.na(d$death_time), "diagnosis_after_death")
  list(data = d[reason == "", , drop = FALSE], exclusions = data.frame(
    Participant_ID = d$Participant_ID,
    reason = ifelse(reason == "", "retained", reason)
  ))
}

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

# Algorithm: depression_delayed_risk ----
depression_delayed_risk <- function(episodes, d, lag_years) {
  ep <- episodes[episodes$from == 3 & episodes$Tstop - episodes$Tstart > lag_years, , drop = FALSE]
  if (!nrow(ep))
    return(data.frame())
  ep$Tstart <- ep$Tstart + lag_years
  long <- expand_risk_sets(ep, d)
  as.data.frame(long[long$trans == 8, , drop = FALSE])
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

# Algorithm: fit_timing_models ----
fit_timing_models <- function(long, spec_id, only_trans = 1:10) {
  rows <- list()
  k <- 0L
  # Each specification/transition requires one mutually adjusted fit and one score fit.
  for (transition in only_trans) {
    z <- as.data.frame(long[long$trans == transition, , drop = FALSE])
    for (model_type in c("individual", "score")) {
      covariates <- if (model_type == "individual") RISK_FACTORS else "lifestyle"
      failure <- ""
      fit <- tryCatch({
        stopifnot(all(c(covariates, "Age_base", "region") %in% names(z)))
        stopifnot(!anyNA(z[, c(covariates, "Age_base", "region"), drop = FALSE]))
        if (!nrow(z) || sum(z$status) == 0) stop("No events/risk set")
        z$region <- factor(z$region)
        formula <- as.formula(paste("Surv(Tstart,Tstop,status) ~",
          paste(c(covariates, "Age_base", "strata(region)", "cluster(Participant_ID)"), collapse = " + ")))
        checked_coxph(formula, z, ties = "efron")
      }, error = function(e) { failure <<- conditionMessage(e); NULL })
      if (is.null(fit)) {
        answer <- data.frame(term = covariates, exposure = covariates, trans = transition,
          beta = NA_real_, SE = NA_real_, HR = NA_real_, lower = NA_real_, upper = NA_real_, P_value = NA_real_)
      } else {
        answer <- coefficient_results(fit, data.frame(term = covariates,
          exposure = covariates, trans = transition))
      }
      answer$spec_id <- spec_id
      answer$status <- if (is.null(fit)) "NOT_ESTIMABLE" else "FITTED"
      answer$message <- failure
      answer$events <- if (nrow(z)) sum(z$status) else 0L
      answer$risk_n <- length(unique(z$Participant_ID))
      answer$person_years <- sum(z$Tstop - z$Tstart)
      answer$time_origin <- "baseline; no clock reset"
      k <- k + 1L
      rows[[k]] <- answer
    }
  }
  answer <- as.data.frame(data.table::rbindlist(rows, use.names = TRUE, fill = TRUE))
  answer$BH_FDR <- NA_real_
  for (index in split(seq_len(nrow(answer)), interaction(answer$spec_id, answer$exposure, drop = TRUE))) {
    answer$BH_FDR[index] <- p.adjust(answer$P_value[index], method = "BH")
  }
  answer$star <- ifelse(is.na(answer$P_value), "", ifelse(answer$P_value < .001, "***",
    ifelse(answer$P_value < .01, "**", ifelse(answer$P_value < .05, "*", ""))))
  answer
}

# Algorithm: transition_matrix ----
transition_matrix <- function() {
  m <- matrix(NA_integer_, 5, 5, dimnames = list(STATE_NAMES, STATE_NAMES))
  m[cbind(TRANSITIONS$from, TRANSITIONS$to)] <- TRANSITIONS$trans
  m
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

# Algorithm: coefficient_results ----
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

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "12"
stage_dependencies <- c("06")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "12_timing_sensitivity.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("06_Analysis_Data.rds", "06_Pre_Tie_Scored.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

data <- readRDS(file.path(results_dir, "06_Analysis_Data.rds"))

pre_tie <- readRDS(file.path(results_dir, "06_Pre_Tie_Scored.rds"))

# Keep sparse estimates explicit; BH is grouped by specification and exposure only.
results <- run_timing_sensitivity(data, pre_tie_data = pre_tie)

fwrite(results$results, file.path(results_dir, "12_S2_Timing_Results.csv"))
fwrite(results$eligibility, file.path(results_dir, "12_S2_Eligibility_Counts.csv"))
fwrite(results$transition_risk, file.path(results_dir, "12_S2_Transition_Risk.csv"))


if (nrow(results$not_run) == 0L) {
  status_table <- data.table(status = "ALL_SPECIFICATIONS_EXECUTED",
    detail = "Inspect coefficient-level estimability in the main result table")
} else {
  status_table <- as.data.table(results$not_run)
}
fwrite(status_table, file.path(results_dir, "12_S2_Not_Run.csv"))

saveRDS(results, file.path(results_dir, "12_Private_Timing_Results.rds"))

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^12_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^12_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
