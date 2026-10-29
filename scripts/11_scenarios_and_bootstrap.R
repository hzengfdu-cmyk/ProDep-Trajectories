# Hypothetical exposure scenarios and paired participant bootstrap.
# Use 1000 participant bootstrap samples and a 15-year prediction horizon.

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

bootstrap_replicates <- 1000L

bootstrap_seed <- 20260810L

# 3. Reused algorithms needed repeatedly in this analysis ----
# Store only model elements needed by prediction, all on the same zero-covariate baseline.
# Converting centered hazards before aggregation prevents mixing different bootstrap centers.
compact_transition_models <- function(models, grid) {
  lapply(models, function(model) {
    fit <- model$fit
    position <- findInterval(grid[-1], model$baseline$time)
    hazard <- c(0, model$baseline$hazard)[position + 1L] * exp(-model$center)
    if (any(!is.finite(hazard)) || any(diff(hazard) < -1e-10)) stop("Invalid uncentered cumulative hazard")
    compact_fit <- structure(list(coefficients = coef(fit), terms = terms(fit),
      contrasts = fit$contrasts, xlevels = fit$xlevels), class = "coxph")
    list(fit = compact_fit, baseline = data.frame(time = grid[-1], hazard = hazard),
      center = 0, events = fit$nevent, warnings = fit$v1_warnings)
  })
}

bootstrap_scenarios <- function(d, B = 1000L, seed = 20260810L,
                               grid = seq(0, 15, .5), checkpoint_dir = NULL) {
  stopifnot(!anyNA(d$Participant_ID), !anyDuplicated(d$Participant_ID), B >= 1L, B == as.integer(B))
  original_models <- fit_scenario_models(d)
  ranking <- rank_single_scenarios(original_models, d, grid)
  defs <- scenario_definitions(ranking)
  defs <- c(defs, setNames(lapply(RISK_FACTORS, identity), paste0("I_", RISK_FACTORS)))
  point <- evaluate_scenarios(original_models, d, defs, grid)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  signature_columns <- c("Participant_ID", TIME_COLUMNS, RISK_FACTORS, "lifestyle", "Age_base", "region")
  signature <- list(data = d[, signature_columns, drop = FALSE], B = B, seed = seed,
    seeds = seeds, RNGkind = RNGkind(), grid = grid, definitions = defs,
    R_version = R.version.string, survival_version = as.character(packageVersion("survival")),
    data_table_version = as.character(packageVersion("data.table")),
    states = STATE_NAMES, transitions = TRANSITIONS, formula = deparse(scenario_model_formula()),
    script_md5 = getOption("mstate.scenario_script_md5", NA_character_),
    algorithm = lapply(c("compact_transition_models", "bootstrap_scenarios", "summarise_bootstrap_coefficients", "fit_scenario_models", "rank_single_scenarios", "scenario_definitions", "evaluate_scenarios", "resample_participants", "summarise_scenario_bootstrap", "build_paths", "expand_risk_sets", "checked_coxph", "scenario_model_formula", "restricted_mean_times", "standardised_sop", "transition_matrix", "validate_grid", "aj_update"), function(name) {
      value <- get(name, envir = environment(bootstrap_scenarios))
      list(name = name, arguments = formals(value), body = deparse(body(value)))
    }))
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    signature_file <- file.path(checkpoint_dir, "Input_and_Algorithm.rds")
    if (file.exists(signature_file)) {
      if (!identical(readRDS(signature_file), signature)) stop("Checkpoint input, algorithm or random settings differ; use a new result directory")
    } else {
      if (length(list.files(checkpoint_dir, pattern = "^replicate_[0-9]+[.]rds$"))) {
        stop("Replicate files exist without a registered input signature")
      }
      saveRDS(signature, signature_file)
    }
  }
  reps <- vector("list", B)
  model_replicates <- vector("list", B)
  logs <- vector("list", B)
  for (b in seq_len(B)) {
    file <- if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, sprintf("replicate_%04d.rds", b))
    if (!is.null(file) && file.exists(file)) {
      iteration <- readRDS(file)
      stopifnot(identical(iteration$iteration, b), identical(iteration$seed, seeds[b]),
        length(iteration$indices) == nrow(d), all(iteration$indices >= 1L),
        all(iteration$indices <= nrow(d)))
    } else {
      set.seed(seeds[b])
      indices <- sample.int(nrow(d), nrow(d), replace = TRUE)
      iteration <- tryCatch({
        models <- fit_scenario_models(resample_participants(d, indices))
        compact_models <- compact_transition_models(models, grid)
        predictions <- evaluate_scenarios(compact_models, d, defs, grid)
        list(success = TRUE, result = predictions, models = compact_models, message = "")
      }, error = function(e) list(success = FALSE, message = conditionMessage(e)))
      iteration$iteration <- b
      iteration$seed <- seeds[b]
      iteration$indices <- indices
      if (!is.null(file)) saveRDS(iteration, file)
    }
    if (iteration$success) {
      reps[[b]] <- iteration$result
      model_replicates[[b]] <- iteration$models
    }
    logs[[b]] <- data.frame(iteration = b, seed = seeds[b], success = iteration$success,
                           message = iteration$message)
  }
  log <- as.data.frame(data.table::rbindlist(logs))
  complete <- B == 1000L && nrow(log) == B && !anyDuplicated(log$iteration) &&
    identical(sort(log$iteration), seq_len(B)) && all(log$success)
  list(point = point, replicates = reps, model_replicates = model_replicates,
    original_models = compact_transition_models(original_models, grid), log = log,
    ranking = ranking, definitions = defs, grid = grid, B_requested = B,
    B_valid = sum(log$success), complete = complete,
    summary = summarise_scenario_bootstrap(reps, defs, grid))
}

summarise_bootstrap_coefficients <- function(model_replicates) {
  records <- list()
  k <- 0L
  for (b in seq_along(model_replicates)) {
    if (is.null(model_replicates[[b]])) next
    for (transition in 1:10) {
      beta <- coef(model_replicates[[b]][[transition]]$fit)
      k <- k + 1L
      records[[k]] <- data.frame(iteration = b, trans = transition,
        exposure = names(beta), beta = unname(beta))
    }
  }
  raw <- as.data.frame(data.table::rbindlist(records))
  if (!nrow(raw)) return(list(raw = raw, summary = data.frame()))
  groups <- split(raw, interaction(raw$trans, raw$exposure, drop = TRUE))
  summary <- lapply(groups, function(x) {
    q <- quantile(x$beta, c(.025, .5, .975), names = FALSE)
    lower_tail <- (sum(x$beta <= 0) + 1) / (nrow(x) + 1)
    upper_tail <- (sum(x$beta >= 0) + 1) / (nrow(x) + 1)
    p <- min(1, 2 * min(lower_tail, upper_tail))
    data.frame(trans = x$trans[1], exposure = x$exposure[1], beta_median = q[2],
      HR = exp(q[2]), lower = exp(q[1]), upper = exp(q[3]),
      bootstrap_tail_P = p, n_valid = nrow(x),
      star = ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
  })
  list(raw = raw, summary = as.data.frame(data.table::rbindlist(summary)))
}

# Algorithm: fit_scenario_models ----
fit_scenario_models <- function(d) {
  d <- {
    analysis_frame <- d
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
  d$region <- factor(d$region, levels = c("England", "Scotland", "Wales"))
  p <- build_paths(d)
  if (!all(p$eligibility$eligible))
    stop("Scenario model input must be an eligible cohort; unresolved ties/boundaries cannot be silently dropped")
  long <- expand_risk_sets(p$episodes, d)
  lapply(1:10, function(tt) {
    z <- as.data.frame(long[long$trans == tt, , drop = FALSE])
    if (!nrow(z) || sum(z$status) == 0)
      stop("No observed event support for transition ", tt)
    fit <- checked_coxph(scenario_model_formula(), z, ties = "efron")
    baseline <- survival::basehaz(fit, centered = TRUE)
    list(fit = fit, baseline = baseline, center = sum(fit$means * stats::coef(fit)))
  })
}

# Algorithm: rank_single_scenarios ----
rank_single_scenarios <- function(models, cohort, grid = seq(0, 15, 0.5)) {
  observed <- restricted_mean_times(standardised_sop(models, cohort, grid = grid))["healthy"]
  gains <- setNames(vapply(RISK_FACTORS, function(v) restricted_mean_times(standardised_sop(
    models,
    cohort, v, grid
  ))["healthy"] - observed, numeric(1)), RISK_FACTORS)
  names(gains)[order(-gains, match(names(gains), RISK_FACTORS))]
}

# Algorithm: scenario_definitions ----
scenario_definitions <- function(ranking) {
  if (!setequal(ranking, RISK_FACTORS) || length(ranking) != 6)
    stop("Need each of six factors once in fixed ranking")
  list(S0 = character(), S1 = ranking[1], S2 = ranking[1:2], S3 = ranking[1:3], S4 = RISK_FACTORS)
}

# Algorithm: evaluate_scenarios ----
evaluate_scenarios <- function(models, cohort, definitions, grid = seq(0, 15, 0.5)) {
  curves <- lapply(definitions, function(v) standardised_sop(models, cohort, v, grid))
  rmts <- lapply(curves, restricted_mean_times)
  deltas <- lapply(curves, function(x) {
    z <- x
    z[, STATE_NAMES] <- 100 * (x[, STATE_NAMES] - curves$S0[, STATE_NAMES])
    z
  })
  list(
    curves = curves, rmts = rmts, years_gained = lapply(rmts, function(x) x - rmts$S0),
    probability_difference_pp = deltas
  )
}

# Algorithm: resample_participants ----
resample_participants <- function(d, indices) {
  z <- d[indices, , drop = FALSE]
  z$Participant_ID <- paste0("bootstrap_", seq_len(nrow(z)))
  rownames(z) <- NULL
  z
}

# Algorithm: summarise_scenario_bootstrap ----
summarise_scenario_bootstrap <- function(replicates, definitions, grid, probabilities = c(
                                           0.025,
                                           0.5, 0.975
                                         )) {
  good <- Filter(Negate(is.null), replicates)
  if (!length(good))
    return(list(years = data.frame(), probabilities = data.frame()))
  q <- function(v) stats::quantile(v, probabilities, names = FALSE)
  years <- list()
  curves <- list()
  k <- 0L
  j <- 0L
  for (s in names(definitions)) {
    for (metric in names(good[[1]]$years_gained[[s]])) {
      z <- q(vapply(good, function(r) r$years_gained[[s]][metric], numeric(1)))
      k <- k + 1L
      years[[k]] <- data.frame(
        scenario = s, metric = metric, lower = z[1], median = z[2],
        upper = z[3], n_valid = length(good)
      )
    }
    for (state in STATE_NAMES) for (i in seq_along(grid)) {
      z <- q(vapply(good, function(r) r$probability_difference_pp[[s]][[state]][i], numeric(1)))
      j <- j + 1L
      curves[[j]] <- data.frame(
        scenario = s, state = state, time = grid[i], lower = z[1],
        median = z[2], upper = z[3], n_valid = length(good), interval = "pointwise"
      )
    }
  }
  list(
    years = as.data.frame(data.table::rbindlist(years, use.names = TRUE, fill = TRUE)),
    probabilities = as.data.frame(data.table::rbindlist(curves, use.names = TRUE, fill = TRUE))
  )
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

# Algorithm: scenario_model_formula ----
scenario_model_formula <- function() stats::as.formula(paste(
  "Surv(Tstart,Tstop,status) ~",
  paste(c(RISK_FACTORS, "Age_base", "region", "cluster(Participant_ID)"), collapse = " + ")
))

# Algorithm: restricted_mean_times ----
restricted_mean_times <- function(sop) {
  stopifnot(all(c("time", STATE_NAMES) %in% names(sop)))
  validate_grid(sop$time)
  P <- as.matrix(sop[, STATE_NAMES])
  if (any(P < -1e-10) || any(abs(rowSums(P) - 1) > 1e-08))
    stop("Invalid state probabilities")
  metrics <- cbind(healthy = P[, 1], pca_free = rowSums(P[, 1:3, drop = FALSE]), depression_free = rowSums(P[,
    c(1, 2, 4),
    drop = FALSE
  ]), alive = rowSums(P[, 1:4, drop = FALSE]))
  colSums((metrics[-1, , drop = FALSE] + metrics[-nrow(metrics), , drop = FALSE]) / 2 * diff(sop$time))
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
  weights <- profiles$weight / sum(profiles$weight)
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

# Algorithm: transition_matrix ----
transition_matrix <- function() {
  m <- matrix(NA_integer_, 5, 5, dimnames = list(STATE_NAMES, STATE_NAMES))
  m[cbind(TRANSITIONS$from, TRANSITIONS$to)] <- TRANSITIONS$trans
  m
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

# 4. Run or resume the registered participant bootstrap. ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "11"
stage_dependencies <- c("06")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "11_scenarios_and_bootstrap.R"), mustWork = TRUE)
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
options(mstate.scenario_script_md5 = unname(stage_code_hash))

data <- readRDS(file.path(results_dir, "06_Analysis_Data.rds"))
bootstrap <- bootstrap_scenarios(data, B = bootstrap_replicates, seed = bootstrap_seed,
  grid = seq(0, 15, .5), checkpoint_dir = file.path(results_dir, "11_Private_Checkpoints"))
saveRDS(bootstrap, file.path(results_dir, "11_Private_Scenario_Bootstrap.rds"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "checks", "11_SessionInfo.txt"))
fwrite(bootstrap$log, file.path(results_dir, "11_Bootstrap_Log.csv"))
fwrite(data.table(requested = bootstrap$B_requested, valid = bootstrap$B_valid,
  complete_1000 = bootstrap$complete), file.path(results_dir, "checks", "11_Bootstrap_Status.csv"))
fwrite(data.table(rank = seq_along(bootstrap$ranking), factor = bootstrap$ranking),
  file.path(results_dir, "11_Single_Factor_Ranking.csv"))
scenario_table <- data.table(scenario = names(bootstrap$definitions),
  replaced_factors = vapply(bootstrap$definitions, paste, character(1), collapse = ";"))
fwrite(scenario_table, file.path(results_dir, "11_Scenario_Definitions.csv"))
if (!bootstrap$complete) {
  stop("Formal outputs require all 1000 registered samples to succeed. Checkpoint failures are retained and are not replaced with new samples")
}
coefficient_output <- summarise_bootstrap_coefficients(bootstrap$model_replicates)
fwrite(coefficient_output$raw, file.path(results_dir, "11_Bootstrap_Coefficients.csv"))
fwrite(coefficient_output$summary, file.path(results_dir, "11_Figure3A_Bootstrap_HR.csv"))
fwrite(bootstrap$summary$years, file.path(results_dir, "11_Years_Gained.csv"))
fwrite(bootstrap$summary$probabilities, file.path(results_dir, "11_Probability_Differences.csv"))
cat("All registered bootstrap samples passed. Formal scenario outputs saved.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^11_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^11_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
manifest[[stage_name]]$settings <- list(seed = bootstrap_seed, RNGkind = RNGkind(),
  definitions = bootstrap$definitions, grid = bootstrap$grid, requested = bootstrap$B_requested,
  valid = bootstrap$B_valid, complete = bootstrap$complete)
saveRDS(manifest, manifest_file)
