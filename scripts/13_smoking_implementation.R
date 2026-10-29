# Partial implementation of the smoking-risk-indicator replacement.
# Mix paired observed and smoking-only trajectories with fixed weights 0.85 and 0.15.

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

# 3. Reused algorithms needed repeatedly in this analysis ----
# Algorithm: smoking_implementation ----
smoking_implementation <- function(bootstrap, p = 0.15) {
  if (length(p) != 1 || !is.finite(p) || p < 0 || p > 1)
    stop("Implementation probability must be in [0,1]")
  good <- Filter(Negate(is.null), bootstrap$replicates)
  if (!length(good))
    stop("No valid paired bootstrap trajectories")
  if (any(!vapply(good, function(r) all(c("S0", "I_smoke") %in% names(r$curves)), logical(1))))
    stop("Paired observed and smoking-only trajectories required")
  years <- lapply(good, function(r) {
    sop <- r$curves$S0
    sop[, STATE_NAMES] <- (1 - p) * r$curves$S0[, STATE_NAMES] + p * r$curves$I_smoke[
      ,
      STATE_NAMES
    ]
    restricted_mean_times(sop) - restricted_mean_times(r$curves$S0)
  })
  year_rows <- lapply(names(years[[1]]), function(metric) {
    q <- stats::quantile(vapply(years, function(r) r[metric], numeric(1)), c(
      0.025, 0.5,
      0.975
    ), names = FALSE)
    data.frame(metric = metric, lower = q[1], median = q[2], upper = q[3])
  })
  curves <- list()
  k <- 0L
  for (state in STATE_NAMES) for (i in seq_along(bootstrap$grid)) {
    delta <- vapply(
      good, function(r) 100 * p * (r$curves$I_smoke[[state]][i] - r$curves$S0[[state]][i]),
      numeric(1)
    )
    q <- stats::quantile(delta, c(0.025, 0.5, 0.975), names = FALSE)
    k <- k + 1L
    curves[[k]] <- data.frame(
      state = state, time = bootstrap$grid[i], lower = q[1], median = q[2],
      upper = q[3], interval = "pointwise"
    )
  }
  list(
    years = as.data.frame(data.table::rbindlist(year_rows, use.names = TRUE, fill = TRUE)),
    curves = as.data.frame(data.table::rbindlist(curves, use.names = TRUE, fill = TRUE)),
    probability = p, n_valid = length(good)
  )
}

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

# Algorithm: validate_grid ----
validate_grid <- function(grid) {
  if (length(grid) < 2 || grid[1] != 0 || any(!is.finite(grid)) || any(diff(grid) <= 0))
    stop("Grid must start at zero and increase strictly")
}

# 4. Use the same complete, paired predictions at both implementation levels. ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "13"
stage_dependencies <- c("11")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "13_smoking_implementation.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("11_Private_Scenario_Bootstrap.rds")),
                         file.path(source_dir, character()))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

bootstrap <- readRDS(file.path(results_dir, "11_Private_Scenario_Bootstrap.rds"))
stopifnot(isTRUE(bootstrap$complete), bootstrap$B_requested == 1000L,
  bootstrap$B_valid == 1000L, length(bootstrap$replicates) == 1000L,
  nrow(bootstrap$log) == 1000L, all(bootstrap$log$success), !anyDuplicated(bootstrap$log$iteration),
  identical(sort(bootstrap$log$iteration), seq_len(1000L)),
  all(vapply(bootstrap$replicates, function(x) !is.null(x), logical(1))))
outputs <- list()
years <- list()
for (probability in c(.15, 1)) {
  implementation <- smoking_implementation(bootstrap, p = probability)
  stopifnot(implementation$n_valid == 1000L)
  implementation$curves$implementation_probability <- probability
  implementation$years$implementation_probability <- probability
  outputs[[as.character(probability)]] <- implementation$curves
  years[[as.character(probability)]] <- implementation$years
}
fwrite(rbindlist(outputs), file.path(results_dir, "13_S4_Probability_Differences.csv"))
fwrite(rbindlist(years), file.path(results_dir, "13_S4_Years_All_Implementation.csv"))
fwrite(years[["0.15"]], file.path(results_dir, "13_S4_Smoking_Years.csv"))
fwrite(outputs[["0.15"]], file.path(results_dir, "13_S4_Smoking_Probability_Differences.csv"))

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^13_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^13_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
