# Historical individual analysis: fixed median model, cohort distribution and a high-benefit example.
# Individual predictions are deterministic conditional profiles, not causal effects or individual confidence intervals.
# Packages and local paths ----
library(data.table)
library(survival)
source_dir <- Sys.getenv("UKB_SOURCE_DIR", unset = file.path(getwd(), "data", "raw"))
results_dir <- Sys.getenv("UKB_RESULTS_DIR", unset = file.path(getwd(), "results"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "checks"), showWarnings = FALSE)
STATE_NAMES <- c("Healthy", "BPH", "Depression", "PCa", "Death")
RISK_FACTORS <- c("smoke", "sport", "sleep", "social_support_score", "new_food_score_score", "sexpartners_score")
LIFESTYLE_VARS <- RISK_FACTORS
PROFILE_VARS <- c("Age_base", "region", RISK_FACTORS)
TRANSITIONS <- data.frame(trans = 1:10, from = c(1,1,1,1,2,2,2,3,3,4), to = c(2,3,4,5,3,4,5,4,5,5))

# Algorithm: evaluate_scenarios ----
evaluate_scenarios <- function(models, cohort, definitions, grid = seq(0, 15, 0.5)) {
    curves <- lapply(definitions, function(v) standardised_sop(models, cohort, v, grid))
    rmts <- lapply(curves, restricted_mean_times)
    deltas <- lapply(curves, function(x) {
        z <- x
        z[, STATE_NAMES] <- 100 * (x[, STATE_NAMES] - curves$S0[, STATE_NAMES])
        z
    })
    list(curves = curves, rmts = rmts, years_gained = lapply(rmts, function(x) x - rmts$S0), 
        probability_difference_pp = deltas)
}

# Algorithm: restricted_mean_times ----
restricted_mean_times <- function(sop) {
    stopifnot(all(c("time", STATE_NAMES) %in% names(sop)))
    validate_grid(sop$time)
    P <- as.matrix(sop[, STATE_NAMES])
    if (any(P < -1e-10) || any(abs(rowSums(P) - 1) > 1e-08)) 
        stop("Invalid state probabilities")
    metrics <- cbind(healthy = P[, 1], pca_free = rowSums(P[, 1:3, drop = FALSE]), depression_free = rowSums(P[, 
        c(1, 2, 4), drop = FALSE]), alive = rowSums(P[, 1:4, drop = FALSE]))
    colSums((metrics[-1, , drop = FALSE] + metrics[-nrow(metrics), , drop = FALSE])/2 * diff(sop$time))
}

# Algorithm: standardised_sop ----
standardised_sop <- function(models, cohort, replace = character(), grid = seq(0, 15, 0.5)) {
    validate_grid(grid)
    if (!length(models) || length(models) != 10) 
        stop("Ten transition models required")
    if (any(!replace %in% RISK_FACTORS)) 
        stop("Invalid scenario replacement")
    d <- {
        analysis_frame <- cohort
        stopifnot(all(c(RISK_FACTORS, "lifestyle") %in% names(analysis_frame)))
        stopifnot(!anyNA(analysis_frame[, c(RISK_FACTORS, "lifestyle"), drop = FALSE]))
        for (column in c(RISK_FACTORS, "lifestyle")) analysis_frame[[column]] <- as.numeric(as.character(analysis_frame[[column]]))
        stopifnot(all(as.matrix(analysis_frame[, RISK_FACTORS]) %in% c(0, 1)))
        stopifnot(all(analysis_frame$lifestyle == rowSums(analysis_frame[, RISK_FACTORS])))
        analysis_frame
    }
    {
        stopifnot(all(c(RISK_FACTORS, "Age_base", "region") %in% names(d)))
        stopifnot(!anyNA(d[, c(RISK_FACTORS, "Age_base", "region"), drop = FALSE]))
    }
    if (!nrow(d)) 
        stop("Empty standardisation cohort")
    for (v in replace) d[[v]] <- 0L
    d$region <- factor(d$region, levels = c("England", "Scotland", "Wales"))
    dt <- data.table::as.data.table(d[, c(RISK_FACTORS, "Age_base", "region"), drop = FALSE])
    profiles <- as.data.frame(dt[, list(weight = .N), by = c(RISK_FACTORS, "Age_base", "region")])
    weights <- profiles$weight/sum(profiles$weight)
    risk <- sapply(models, function(m) {
        fit <- m$fit
        X <- stats::model.matrix(stats::delete.response(stats::terms(fit)), profiles, contrasts.arg = fit$contrasts)
        X <- X[, names(stats::coef(fit)), drop = FALSE]
        as.numeric(exp(X %*% stats::coef(fit) - m$center))
    })
    risk <- matrix(risk, nrow = nrow(profiles), ncol = 10)
    if (any(!is.finite(risk))) 
        stop("Nonfinite scenario relative hazard")
    cum <- sapply(models, function(m) {
        h <- m$baseline
        ix <- findInterval(grid, h$time)
        c(0, h$hazard)[ix + 1L]
    })
    if (any(abs(cum[1, ]) > 1e-12)) 
        stop("Baseline hazard before time zero is nonzero")
    inc <- matrix(apply(cum, 2, diff), nrow = length(grid) - 1L, ncol = 10)
    probabilities <- matrix(0, nrow(profiles), 5)
    probabilities[, 1] <- 1
    out <- matrix(0, length(grid), 5, dimnames = list(NULL, STATE_NAMES))
    out[1, 1] <- 1
    for (i in seq_len(length(grid) - 1L)) {
        probabilities <- aj_update(probabilities, sweep(risk, 2, inc[i, ], `*`))
        out[i + 1, ] <- as.numeric(crossprod(weights, probabilities))
    }
    data.frame(time = grid, out, check.names = FALSE)
}

# Algorithm: validate_grid ----
validate_grid <- function(grid) {
    if (length(grid) < 2 || grid[1] != 0 || any(!is.finite(grid)) || any(diff(grid) <= 0)) 
        stop("Grid must start at zero and increase strictly")
}

# Algorithm: aj_update ----
aj_update <- function(probabilities, increments) {
    p <- as.matrix(probabilities)
    a <- as.matrix(increments)
    if (ncol(p) != 5 || ncol(a) != 10 || nrow(p) != nrow(a)) 
        stop("AJ matrix dimensions differ")
    if (any(!is.finite(a)) || any(a < 0)) 
        stop("Invalid hazard increments")
    old <- p
    for (tt in 1:10) {
        f <- TRANSITIONS$from[tt]
        to <- TRANSITIONS$to[tt]
        flow <- old[, f] * a[, tt]
        p[, f] <- p[, f] - flow
        p[, to] <- p[, to] + flow
    }
    if (any(!is.finite(p)) || any(p < -1e-10) || any(abs(rowSums(p) - 1) > 1e-08)) 
        stop("Half-year AJ update invalid: aggregated outgoing hazards are too large; do not clip or silently renormalise")
    p[abs(p) < 1e-14] <- 0
    p
}

# Algorithm: derive_all_lifestyle_ranking ----
derive_all_lifestyle_ranking <- function(yg_summary) {
    required <- c("intervention_type", "intervention", "intervention_label", "metric", "median_yg")
    missing <- setdiff(required, names(yg_summary))
    if (length(missing)) {
        stop("YG summary is missing: ", paste(missing, collapse = ", "))
    }
    selected <- yg_summary[yg_summary$intervention_type == "single" & yg_summary$metric == 
        "healthy", required[required %in% names(yg_summary)], drop = FALSE]
    if (nrow(selected) != 6L || anyDuplicated(selected$intervention) || !setequal(selected$intervention, 
        LIFESTYLE_VARS) || any(!is.finite(selected$median_yg))) {
        stop("Healthy YG ranking requires exactly six finite single interventions.")
    }
    selected <- selected[order(-selected$median_yg, selected$intervention), ]
    rownames(selected) <- NULL
    selected$healthy_yg_rank <- seq_len(nrow(selected))
    selected[c("intervention", "intervention_label", "median_yg", "healthy_yg_rank")]
}

# Algorithm: build_patient_candidate_sets ----
build_patient_candidate_sets <- function(cohort, top_three) {
    top_three <- as.character(top_three)
    required <- c("Participant_ID", "Age_base", top_three, "path1", "path5", "path8", "BPH_time", 
        "depress_time", "PCa_time")
    missing <- setdiff(required, names(cohort))
    if (length(missing)) {
        stop("Candidate data are missing: ", paste(missing, collapse = ", "))
    }
    if (length(top_three) != 3L || anyDuplicated(top_three) || any(!top_three %in% LIFESTYLE_VARS)) {
        stop("top_three must contain three unique model lifestyle variables.")
    }
    cohort$Participant_ID <- as.character(cohort$Participant_ID)
    age_eligible <- cohort[is.finite(cohort$Age_base) & cohort$Age_base >= 65, , drop = FALSE]
    top_three_eligible <- age_eligible[rowSums(age_eligible[top_three] == 1L) == 3L, , drop = FALSE]
    partial_flag <- top_three_eligible$path1 == 1L & top_three_eligible$path5 == 1L & is.finite(top_three_eligible$BPH_time) & 
        is.finite(top_three_eligible$depress_time) & top_three_eligible$BPH_time < top_three_eligible$depress_time
    strict_flag <- partial_flag & top_three_eligible$path8 == 1L & is.finite(top_three_eligible$PCa_time) & 
        top_three_eligible$depress_time < top_three_eligible$PCa_time
    list(all = cohort, age_eligible = age_eligible, top_three_eligible = top_three_eligible, 
        partial = top_three_eligible[partial_flag, , drop = FALSE], strict = top_three_eligible[strict_flag, 
            , drop = FALSE], top_three = top_three)
}

# Algorithm: select_typical_patient ----
select_typical_patient <- function(candidate_sets, candidate_healthy_yg) {
    if (nrow(candidate_sets$strict)) {
        candidates <- candidate_sets$strict
        level <- "strict_h_bph_depression_pca"
    }
    else if (nrow(candidate_sets$partial)) {
        candidates <- candidate_sets$partial
        level <- "fallback_h_bph_depression"
    }
    else if (nrow(candidate_sets$top_three_eligible)) {
        candidates <- candidate_sets$top_three_eligible
        level <- "fallback_age_top_three"
    }
    else {
        stop("No participant remains after the documented selection fallbacks.")
    }
    supplied_names <- names(candidate_healthy_yg)
    candidate_healthy_yg <- as.numeric(candidate_healthy_yg)
    names(candidate_healthy_yg) <- supplied_names
    ids <- as.character(candidates$Participant_ID)
    supplied <- candidate_healthy_yg
    if (is.null(names(supplied)) || any(!ids %in% names(supplied))) {
        stop("Healthy YG values must be named for every candidate Participant ID.")
    }
    candidates$healthy_yg_s4 <- as.numeric(supplied[ids])
    if (any(!is.finite(candidates$healthy_yg_s4))) {
        stop("Candidate Healthy YG values must be finite.")
    }
    pca_tie_break <- ifelse(is.finite(candidates$PCa_time), candidates$PCa_time, Inf)
    order_index <- order(-candidates$healthy_yg_s4, pca_tie_break, as.character(candidates$Participant_ID))
    ranked <- candidates[order_index, , drop = FALSE]
    ranked$selection_rank <- seq_len(nrow(ranked))
    patient <- ranked[1, , drop = FALSE]
    patient$selection_level <- level
    list(patient = patient, ranked_candidates = ranked, selection_level = level, fallback_used = !identical(level, 
        "strict_h_bph_depression_pca"))
}

# Algorithm: empirical_percentile ----
empirical_percentile <- function(population_yg, participant_id) {
    participant_id <- as.character(participant_id)
    if (is.null(names(population_yg)) || !participant_id %in% names(population_yg)) {
        stop("The typical patient is missing from the population YG distribution.")
    }
    finite_values <- population_yg[is.finite(population_yg)]
    target <- unname(population_yg[participant_id])
    if (!length(finite_values) || !is.finite(target)) {
        stop("Healthy YG percentile cannot be calculated from non-finite values.")
    }
    mean(finite_values <= target) * 100
}

# Build the fixed model on one common zero-covariate hazard reference.
build_central_model <- function(model_replicates, grid) {
  good <- Filter(Negate(is.null), model_replicates)
  if (!length(good)) stop("No complete transition models available")
  lapply(1:10, function(transition) {
    template <- good[[1]][[transition]]
    coefficient_names <- names(coef(template$fit))
    beta <- vapply(good, function(models) {
      stopifnot(setequal(names(coef(models[[transition]]$fit)), coefficient_names),
                identical(models[[transition]]$center, 0))
      as.numeric(coef(models[[transition]]$fit)[coefficient_names])
    }, numeric(length(coefficient_names)))
    median_beta <- apply(beta, 1L, median)
    names(median_beta) <- coefficient_names
    cumulative <- vapply(good, function(models) {
      baseline <- models[[transition]]$baseline
      position <- findInterval(grid[-1], baseline$time)
      c(0, baseline$hazard)[position + 1L]
    }, numeric(length(grid) - 1L))
    if (any(!is.finite(cumulative)) || any(cumulative < 0)) stop("Invalid common-reference baseline hazards")
    median_hazard <- apply(cumulative, 1L, median)
    stopifnot(all(diff(median_hazard) >= -1e-10))
    median_hazard <- cummax(pmax(median_hazard, 0))
    template$fit$coefficients <- median_beta
    template$baseline <- data.frame(time = grid[-1], hazard = median_hazard)
    template$center <- 0
    template
  })
}

# Integrate all compressed profiles together, retaining only healthy-state areas.
# This is the historical full-profile calculation with lower memory usage.
profile_healthy_time <- function(models, profiles, grid) {
  validate_grid(grid)
  profiles$region <- factor(profiles$region, levels = c("England", "Scotland", "Wales"))
  risk <- sapply(models, function(model) {
    matrix <- model.matrix(delete.response(terms(model$fit)), profiles,
                           contrasts.arg = model$fit$contrasts)
    matrix <- matrix[, names(coef(model$fit)), drop = FALSE]
    as.numeric(exp(matrix %*% coef(model$fit) - model$center))
  })
  risk <- matrix(risk, nrow = nrow(profiles), ncol = 10L)
  cumulative <- sapply(models, function(model) {
    position <- findInterval(grid, model$baseline$time)
    c(0, model$baseline$hazard)[position + 1L]
  })
  increments <- matrix(apply(cumulative, 2L, diff), nrow = length(grid) - 1L, ncol = 10L)
  probability <- matrix(0, nrow(profiles), 5L)
  probability[, 1] <- 1
  area <- numeric(nrow(profiles))
  for (i in seq_len(length(grid) - 1L)) {
    next_probability <- aj_update(probability, sweep(risk, 2L, increments[i, ], `*`))
    area <- area + (probability[, 1] + next_probability[, 1]) * (grid[i+1L] - grid[i]) / 2
    probability <- next_probability
  }
  area
}

# Read a complete registered prediction-model bootstrap, not association-model HRs.
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "11b"
stage_dependencies <- c("11", "07")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "11b_individual_predictions.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("11_Private_Scenario_Bootstrap.rds", "06_Analysis_Data.rds", "07_Private_Paths.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

bootstrap <- readRDS(file.path(results_dir, "11_Private_Scenario_Bootstrap.rds"))
stopifnot(isTRUE(bootstrap$complete), bootstrap$B_requested == 1000L,
  bootstrap$B_valid == 1000L, length(bootstrap$model_replicates) == 1000L,
  all(vapply(bootstrap$model_replicates, function(x) !is.null(x), logical(1))),
  nrow(bootstrap$log) == 1000L, all(bootstrap$log$success),
  identical(sort(bootstrap$log$iteration), seq_len(1000L)))
data <- as.data.frame(readRDS(file.path(results_dir, "06_Analysis_Data.rds")))
stopifnot(!anyNA(data$Participant_ID), !anyDuplicated(data$Participant_ID))
central_model <- build_central_model(bootstrap$model_replicates, bootstrap$grid)
saveRDS(central_model, file.path(results_dir, "11b_Private_Central_Model.rds"))

# The historical eligibility ranking uses median single-factor healthy gains.
single <- as.data.frame(bootstrap$summary$years)
single <- single[startsWith(single$scenario, "I_") & single$metric == "healthy", , drop = FALSE]
historical_ranking_input <- data.frame(intervention_type = "single", metric = "healthy",
  intervention = sub("^I_", "", single$scenario),
  intervention_label = sub("^I_", "", single$scenario), median_yg = single$median)
ranking <- derive_all_lifestyle_ranking(historical_ranking_input)
top_three <- as.character(ranking$intervention[1:3])
fwrite(ranking, file.path(results_dir, "11b_Individual_Eligibility_Ranking.csv"))

# Compress identical observed and all-low-risk profiles; keep participant alignment explicit.
observed_profiles <- data[, PROFILE_VARS, drop = FALSE]
low_profiles <- observed_profiles
for (factor in RISK_FACTORS) low_profiles[[factor]] <- 0L
profiles <- unique(rbind(observed_profiles, low_profiles))
observed_key <- do.call(paste, c(lapply(observed_profiles, as.character), sep = "\r"))
low_key <- do.call(paste, c(lapply(low_profiles, as.character), sep = "\r"))
profile_key <- do.call(paste, c(lapply(profiles, as.character), sep = "\r"))
observed_index <- match(observed_key, profile_key)
low_index <- match(low_key, profile_key)
stopifnot(!anyNA(observed_index), !anyNA(low_index))
healthy_time <- profile_healthy_time(central_model, profiles, bootstrap$grid)
population_gains <- setNames(healthy_time[low_index] - healthy_time[observed_index],
                            as.character(data$Participant_ID))
stopifnot(all(is.finite(population_gains)))
saveRDS(population_gains, file.path(results_dir, "11b_Private_Population_Healthy_Gains.rds"))

# Reconstruct event membership from the common risk-set output, not historical input path columns.
paths <- readRDS(file.path(results_dir, "07_Private_Paths.rds"))
for (transition in c(1L, 5L, 8L)) {
  ids <- paths$events$Participant_ID[paths$events$trans == transition]
  data[[paste0("path", transition)]] <- as.integer(data$Participant_ID %in% ids)
}
candidates <- build_patient_candidate_sets(data, top_three)
selection <- select_typical_patient(candidates, population_gains)
patient <- selection$patient
patient_id <- as.character(patient$Participant_ID[1])
percentile <- empirical_percentile(population_gains, patient_id)
if (selection$fallback_used) warning("Historical fallback eligibility was used; the figure caption must show the actual selection tier")
saveRDS(selection, file.path(results_dir, "11b_Private_Selection.rds"))
fwrite(data.table(display_label = "Illustrative high-benefit profile",
  healthy_S4_gain = unname(population_gains[patient_id]), percentile = percentile,
  selection_level = selection$selection_level, fallback_used = selection$fallback_used,
  eligibility_factors = paste(top_three, collapse = ";"),
  fixed_model = "Component-wise median coefficients and pointwise median uncentered cumulative hazards",
  individual_intervals = "None; deterministic fixed-model prediction"),
  file.path(results_dir, "11b_Illustrative_Profile_Summary.csv"))

# The frozen scenario definitions remain those used by the population analysis.
definitions <- bootstrap$definitions[paste0("S", 0:4)]
individual <- evaluate_scenarios(central_model, patient, definitions, bootstrap$grid)
curves <- list()
years <- list()
exposures <- list()
for (scenario in names(definitions)) {
  curve <- as.data.table(individual$probability_difference_pp[[scenario]])
  curve <- melt(curve, id.vars = "time", measure.vars = STATE_NAMES,
                variable.name = "state", value.name = "difference_pp")
  curve[, scenario := scenario]
  curves[[scenario]] <- curve
  value <- individual$years_gained[[scenario]]
  years[[scenario]] <- data.table(scenario = scenario, metric = names(value), years_gained = unname(value))
  profile <- patient[, RISK_FACTORS, drop = FALSE]
  for (factor in definitions[[scenario]]) profile[[factor]] <- 0L
  exposures[[scenario]] <- data.table(scenario = scenario, factor = RISK_FACTORS,
    higher_risk = as.numeric(unlist(profile[1, RISK_FACTORS], use.names = FALSE)))
}
fwrite(rbindlist(curves), file.path(results_dir, "11b_Illustrative_Profile_Contrasts.csv"))
fwrite(rbindlist(years), file.path(results_dir, "11b_Illustrative_Profile_Years.csv"))
fwrite(rbindlist(exposures), file.path(results_dir, "11b_Illustrative_Profile_Exposures.csv"))
# Histogram counts can be plotted without publishing participant identifiers.
histogram <- hist(unname(population_gains), breaks = 30, plot = FALSE)
fwrite(data.table(left = head(histogram$breaks, -1), right = tail(histogram$breaks, -1),
  midpoint = histogram$mids, count = histogram$counts),
  file.path(results_dir, "11b_Individual_Distribution_Bins.csv"))
cat("Historical individual analysis complete. Public plot tables use an illustrative alias, not a participant ID.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^11b_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^11b_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
