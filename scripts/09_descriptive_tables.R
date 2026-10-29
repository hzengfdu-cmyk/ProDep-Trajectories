# Baseline tables, incidence rates, state durations and annual observed stock.
# Transition rates and the Methods-specific hierarchical endpoint rates have different denominators.

# 1. Packages and local paths ----
library(data.table)
library(survival)
source_dir <- Sys.getenv('UKB_SOURCE_DIR', unset = file.path(getwd(), 'data', 'raw'))
results_dir <- Sys.getenv('UKB_RESULTS_DIR', unset = file.path(getwd(), 'results'))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, 'checks'), showWarnings = FALSE)
options(stringsAsFactors = FALSE)
data.table::setDTthreads(4L)

# S8 describes participants already recruited and still under observation at the year boundary.

# 2. Definitions used by this analysis ----
STATE_NAMES <- c("Healthy", "BPH", "Depression", "PCa", "Death")

TRANSITIONS <- data.frame(trans = 1:10, from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4), to = c(
  2, 3,
  4, 5, 3, 4, 5, 4, 5, 5
))

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "09"
stage_dependencies <- c("06", "07")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "09_descriptive_tables.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("06_Analysis_Data.rds", "07_Private_Paths.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

data <- as.data.table(readRDS(file.path(results_dir, "06_Analysis_Data.rds")))

paths <- readRDS(file.path(results_dir, "07_Private_Paths.rds"))

episodes <- as.data.table(paths$episodes)

continuous <- c("Age_base", "Townsend", "BMI", "new_food_score", "sexpartners")

categorical <- c(
  "Race", "region", "Education", "sport", "sleep", "smoke", "social_support_score",
  "sexpartners_score", "lifestyle"
)

table1 <- list()

normality <- list()

for (variable in continuous) {
  value <- data[[variable]]
  value <- value[!is.na(value)]
  if (length(value) == 0L) {
    table1[[variable]] <- data.table(
      variable = variable, level = "", n = NA_integer_,
      value = "Unavailable in supplied source"
    )
    next
  }
  set.seed(20261029L)
  sample_value <- if (length(value) > 5000L)
    sample(value, 5000L)
  else value
  sd_value <- sd(value)
  shapiro_p <- if (length(sample_value) >= 3L && sd(sample_value) > 0)
    shapiro.test(sample_value)$p.value
  else NA_real_
  skewness <- if (sd_value > 0)
    mean((value - mean(value))^3) / sd_value^3
  else NA_real_
  kurtosis <- if (sd_value > 0)
    mean((value - mean(value))^4) / sd_value^4 - 3
  else NA_real_
  use_mean <- variable == "Age_base" || (!is.na(shapiro_p) && shapiro_p >= 0.05 && abs(skewness) <=
    1 && abs(kurtosis) <= 1)
  summary_value <- if (use_mean)
    sprintf("%.2f (%.2f)", mean(value), sd_value)
  else {
    q <- quantile(value, c(0.25, 0.5, 0.75))
    sprintf("%.2f (%.2f, %.2f)", q[2], q[1], q[3])
  }
  table1[[variable]] <- data.table(variable = variable, level = "", n = length(value), value = summary_value)
  normality[[variable]] <- data.table(
    variable = variable, shapiro_p = shapiro_p, skewness = skewness,
    excess_kurtosis = kurtosis, method = if (use_mean)
      "mean_SD"
    else "median_IQR"
  )
}

for (variable in categorical) {
  value <- as.character(data[[variable]])
  if (variable == "smoke") {
    value[value == "0"] <- "Never smoked"
    value[value == "1"] <- "Ever smoked"
  }
  denominator <- sum(!is.na(value))
  if (denominator == 0L) {
    table1[[variable]] <- data.table(
      variable = variable, level = "", n = NA_integer_,
      value = "Unavailable in supplied source"
    )
    next
  }
  counts <- as.data.table(table(value, useNA = "no"))
  setnames(counts, c("level", "n"))
  counts[, `:=`(variable = variable, value = sprintf("%d (%.1f%%)", n, n / denominator * 100))]
  table1[[variable]] <- counts[, .(variable, level, n, value)]
}

fwrite(rbindlist(table1, fill = TRUE), file.path(results_dir, "09_Table1.csv"))

fwrite(rbindlist(normality, fill = TRUE), file.path(results_dir, "checks", "09_Normality.csv"))

missing_covariates <- c("BMI", "Race", "Education")[vapply(
  data[, .(BMI, Race, Education)],
  function(x) all(is.na(x)), logical(1)
)]

fwrite(data.table(
  status = if (length(missing_covariates)) "INCOMPLETE_SOURCE" else "COMPLETE",
  missing_fields = paste(missing_covariates, collapse = ";")
), file.path(
  results_dir, "checks",
  "09_Table1_Status.csv"
))

# Transition rates use actual time at risk in the originating state.
tir <- list()

for (transition in 1:10) {
  risk <- episodes[from == TRANSITIONS$from[transition]]
  events <- sum(risk$actual_trans == transition, na.rm = TRUE)
  person_years <- sum(risk$Tstop - risk$Tstart)
  estimate <- lower <- upper <- NA_real_
  if (person_years > 0) {
    confidence <- poisson.test(events, T = person_years)$conf.int * 1000
    estimate <- events / person_years * 1000
    lower <- confidence[1]
    upper <- confidence[2]
  }
  tir[[transition]] <- data.table(
    trans = transition, events = events, person_years = person_years,
    rate_per_1000 = estimate, lower = lower, upper = upper
  )
}

fwrite(rbindlist(tir), file.path(results_dir, "09_S5_Transition_Rates.csv"))

# Endpoint rates follow the separate hierarchical person-time definition in the Methods.
base_time <- episodes[from == 1, .(Participant_ID, baseline_stop = Tstop)]

ir <- list()

for (target in 2:5) {
  observed <- episodes[!is.na(actual_trans) & from %in% seq_len(target - 1L), .(last_transition_time = max(Tstop)),
    by = Participant_ID
  ]
  contribution <- merge(base_time, observed, by = "Participant_ID", all.x = TRUE)
  contribution[is.na(last_transition_time), `:=`(last_transition_time, baseline_stop)]
  person_years <- sum(contribution$last_transition_time)
  events <- sum(episodes$actual_trans %in% TRANSITIONS$trans[TRANSITIONS$to == target])
  confidence <- poisson.test(events, T = person_years)$conf.int * 1000
  ir[[target]] <- data.table(
    endpoint = STATE_NAMES[target], events = events, person_years = person_years,
    rate_per_1000 = events / person_years * 1000, lower = confidence[1], upper = confidence[2]
  )
}

fwrite(rbindlist(ir), file.path(results_dir, "09_S6_Endpoint_Rates.csv"))

episodes[, `:=`(entry_age = Age_base + Tstart, exit_age = Age_base + Tstop, duration = Tstop -
  Tstart, transition_age = ifelse(is.na(actual_trans), NA_real_, Age_base + Tstop))]

state_summary <- list()

for (measure in c("entry_age", "exit_age", "duration", "transition_age")) {
  summary <- episodes[!is.na(get(measure)), .(
    n = .N, mean = mean(get(measure)), SD = sd(get(measure)),
    median = median(get(measure)), p25 = quantile(get(measure), 0.25), p75 = quantile(
      get(measure),
      0.75
    )
  ), by = from]
  summary[, `:=`(measure, measure)]
  state_summary[[measure]] <- summary
}

fwrite(rbindlist(state_summary), file.path(results_dir, "09_State_Interval_Summary.csv"))

# Formal S7 uses the same event population for entry age, transition age and duration.
actual_events <- episodes[!is.na(actual_trans)]
S7 <- actual_events[, .(events = .N,
  entry_age_mean = mean(entry_age), entry_age_SD = sd(entry_age),
  transition_age_mean = mean(transition_age), transition_age_SD = sd(transition_age),
  duration_mean = mean(duration), duration_SD = sd(duration)), by = actual_trans]
S7 <- merge(data.table(actual_trans = 1:10), S7, by = "actual_trans", all.x = TRUE, sort = TRUE)
S7[is.na(events), events := 0L]
fwrite(S7, file.path(results_dir, "09_S7_Transition_Summary.csv"))


# S8 uses stock immediately before events at Jan 1 and calendar flows [Jan 1, next Jan 1).
# A participant recruited on Jan 1 enters that year's recruitment flow, not the pre-boundary stock.
years <- seq(as.integer(format(min(data$baseline_date), "%Y")),
             as.integer(format(max(data$censor_date), "%Y")))
annual <- list()
for (year in years) {
  start <- as.numeric(as.Date(paste0(year, "-01-01")) - data$baseline_date) / 365.25
  end <- as.numeric(as.Date(paste0(year + 1L, "-01-01")) - data$baseline_date) / 365.25
  death <- data$death_time
  censor <- data$censor_time
  death_observed <- !is.na(death) & death <= censor
  observation_end <- ifelse(death_observed, death, censor)
  alive_start <- start > 0 & observation_end >= start
  alive_next <- end > 0 & observation_end >= end
  new_entries <- start <= 0 & end > 0
  death_flow <- death_observed & death >= pmax(0, start) & death < end
  censor_flow <- !death_observed & censor >= pmax(0, start) & censor < end
  events <- numeric(4)
  for (target in 2:5) {
    actual <- episodes[!is.na(actual_trans) & actual_trans %in% TRANSITIONS$trans[TRANSITIONS$to == target]]
    index <- match(actual$Participant_ID, data$Participant_ID)
    events[target - 1L] <- sum(actual$Tstop >= pmax(0, start[index]) & actual$Tstop < end[index])
  }
  stopifnot(events[4] == sum(death_flow))
  age <- data$Age_base[alive_start] + start[alive_start]
  annual[[as.character(year)]] <- data.table(
    Year = year, Interval = "[Jan 1, next Jan 1)",
    Stock_timing = "Immediately before Jan 1 recruitment, death and censoring events",
    BPH_events = events[1], Depression_events = events[2], PCa_events = events[3],
    Death_events_calendar = events[4], Alive_at_risk_start = sum(alive_start),
    New_entries = sum(new_entries), Deaths_flow = sum(death_flow), Censorings_flow = sum(censor_flow),
    Alive_at_risk_next = sum(alive_next),
    Flow_residual = sum(alive_start) + sum(new_entries) - sum(death_flow) - sum(censor_flow) - sum(alive_next),
    Mean_age = if (length(age)) mean(age) else NA_real_,
    SD_age = if (length(age) > 1L) sd(age) else NA_real_,
    Median_age = if (length(age)) median(age) else NA_real_,
    P25_age = if (length(age)) quantile(age, .25, names = FALSE) else NA_real_,
    P75_age = if (length(age)) quantile(age, .75, names = FALSE) else NA_real_)
}

annual <- rbindlist(annual)

stopifnot(all(annual$Flow_residual == 0))

fwrite(annual, file.path(results_dir, "09_S8_Annual_Observed_Stock.csv"))
fwrite(data.table(status = "DEFINITION_CONFIRMED",
  population = "Already recruited, alive and still under observation at the year boundary",
  annual_window = "[Jan 1, next Jan 1)",
  stock_timing = "Immediately before Jan 1 events"),
  file.path(results_dir, "checks", "09_S8_Status.csv"))

cat("Step 09 complete. Table1 status and annual stock-flow checks are explicit.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^09_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^09_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
