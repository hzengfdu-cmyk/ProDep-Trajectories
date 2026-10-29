# Optional plotting only; no statistical calculation depends on this script.
# Plot saved statistical outputs. No models are fitted in this part.
# Tables and model objects remain the source for checking estimates and sparse results.
# Packages and local paths ----
library(data.table)
library(ggplot2)
results_dir <- Sys.getenv('UKB_RESULTS_DIR', unset = file.path(getwd(), 'results'))
figure_dir <- file.path(results_dir, 'figures')
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "14"
stage_dependencies <- c("08a", "08b", "08c", "08d", "10", "11b", "12", "13")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "14_plot_figures.R"), mustWork = TRUE)
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
stage_input_files <- file.path(results_dir, c(
  "08_individual_HR.csv",
  "08_score_categorical_HR.csv",
  "08_score_continuous_HR.csv",
  "08_components_HR.csv",
  "10_BPH_Initial_State_Probabilities.csv",
  "checks/11_Bootstrap_Status.csv",
  "11_Figure3A_Bootstrap_HR.csv",
  "11_Years_Gained.csv",
  "11_Probability_Differences.csv",
  "12_S2_Timing_Results.csv",
  "13_S4_Probability_Differences.csv",
  "11b_Individual_Distribution_Bins.csv",
  "11b_Illustrative_Profile_Summary.csv",
  "11b_Illustrative_Profile_Exposures.csv",
  "11b_Illustrative_Profile_Contrasts.csv"
))
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

# Primary sequence displays; complete transition estimates remain in the CSV files.
for (model in c('individual','score_categorical','score_continuous','components')) {
  estimates <- fread(file.path(results_dir, paste0('08_', model, '_HR.csv')))
  estimates <- estimates[trans %in% c(1,5,8,10)]
  stopifnot(nrow(estimates) > 0L, all(is.finite(estimates$HR)),
            all(is.finite(estimates$lower)), all(is.finite(estimates$upper)))
  if (model == "components") {
    estimates[, panel := paste0(exposure, "\n", contrast)]
  } else estimates[, panel := exposure]
  plot <- ggplot(estimates, aes(x = HR, y = factor(trans))) +
    geom_vline(xintercept = 1, linetype = 2, colour = 'grey50') +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = 'y', width = 0.15) +
    geom_point(size = 1.8) + scale_x_log10() + facet_wrap(~panel, scales = 'free_x') +
    labs(x = 'Hazard ratio (95% CI)', y = 'Transition', title = model) + theme_classic()
  if (model == "components") {
    plot <- plot + geom_text(aes(label = paste0("events=", display_events)), vjust = 1.7, size = 2.4) +
      labs(caption = "Binary factors: exposed events; weighted continuous scores: total events")
  }
  if (model != "components") {
    plot <- plot + geom_text(aes(label = ifelse(is.na(exposed_events),
      paste0("events=", total_events), paste0("events=", exposed_events, "; ref=", reference_events))),
      vjust = 1.7, size = 2.4)
  }
  ggsave(file.path(figure_dir, paste0('HR_', model, '.pdf')), plot, width = 11, height = 7)
}

# Probabilities start in BPH at baseline, as specified for this analysis.
probabilities <- fread(file.path(results_dir, '10_BPH_Initial_State_Probabilities.csv'))
for (column in paste0('pstate',1:5)) {
  set(probabilities, j = column, value = as.numeric(probabilities[[column]]))
}
probabilities <- melt(probabilities, id.vars = c('time','stratum'),
                      measure.vars = paste0('pstate',1:5), variable.name = 'state', value.name = 'probability')
probabilities[, state := factor(state, levels = paste0('pstate',1:5),
                                labels = c('Healthy','BPH','Depression','PCa','Death'))]
plot <- ggplot(probabilities, aes(time, probability, colour = state)) + geom_line() +
  facet_wrap(~stratum) + labs(x = 'Years since baseline', y = 'State occupancy probability',
                            colour = 'State', subtitle = 'Initial prediction state: BPH at baseline time zero') + theme_classic()
ggsave(file.path(figure_dir,'State_Occupancy.pdf'), plot, width = 11, height = 7)

# Formal scenario panels require the complete registered bootstrap, not only a nonempty table.
status <- fread(file.path(results_dir, "checks", "11_Bootstrap_Status.csv"))
stopifnot(nrow(status) == 1L, isTRUE(status$complete_1000[1]),
          status$requested[1] == 1000L, status$valid[1] == 1000L)

# Figure 3A uses prediction-model bootstrap coefficients, not association-model HRs.
coefficient_summary <- fread(file.path(results_dir, "11_Figure3A_Bootstrap_HR.csv"))
coefficient_summary <- coefficient_summary[exposure %in% c("smoke", "sport", "sleep", "social_support_score", "new_food_score_score", "sexpartners_score")]
stopifnot(all(coefficient_summary$n_valid == 1000L), nrow(coefficient_summary) == 60L)
plot <- ggplot(coefficient_summary, aes(factor(trans), exposure, fill = log(HR))) +
  geom_tile(colour = "white") + geom_text(aes(label = paste0(sprintf("%.2f", HR), star)), size = 3) +
  scale_fill_gradient2(low = "#3978A8", mid = "white", high = "#B54B3A", midpoint = 0) +
  labs(x = "Transition", y = "Baseline risk factor", fill = "log HR",
       caption = "HR = exp(median bootstrap coefficient); stars use nominal two-sided bootstrap tails") + theme_classic()
ggsave(file.path(figure_dir, "Figure3A_Prediction_Model_HR.pdf"), plot, width = 11, height = 5)

# Figure 3B: all six single-factor scenarios and all four time measures.
single <- fread(file.path(results_dir, "11_Years_Gained.csv"))[startsWith(scenario, "I_")]
stopifnot(nrow(single) == 24L, all(single$n_valid == 1000L))
plot <- ggplot(single, aes(scenario, median)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .15) + geom_point() +
  facet_wrap(~metric, scales = "free_y") + labs(x = "Single-factor scenario", y = "Years gained versus S0") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figure_dir, "Figure3B_Single_Factor_Years.pdf"), plot, width = 12, height = 8)

# Joint scenario time contrasts are paired within the bootstrap.
years <- fread(file.path(results_dir,'11_Years_Gained.csv'))[scenario %in% paste0('S',0:4)]
plot <- ggplot(years, aes(scenario, median)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) + geom_point() +
  facet_wrap(~metric, scales = 'free_y') + labs(x = 'Hypothetical exposure scenario', y = 'Years gained versus S0') + theme_classic()
ggsave(file.path(figure_dir,'Scenario_Years_Gained.pdf'), plot, width = 10, height = 7)

# These are state-occupancy differences, not differences in cancer incidence risk.
contrasts <- fread(file.path(results_dir,'11_Probability_Differences.csv'))[scenario %in% paste0('S',1:4)]
plot <- ggplot(contrasts, aes(time, median, colour = scenario, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) + geom_line() +
  facet_wrap(~state, scales = 'free_y') + labs(x = 'Years since baseline',
      y = 'State occupancy difference (percentage points)', caption = 'Pointwise bootstrap intervals') + theme_classic()
ggsave(file.path(figure_dir,'Scenario_Probability_Differences.pdf'), plot, width = 11, height = 7)

# Timing stars use unadjusted P. Non-estimable rows remain explicit in the source table.
timing <- fread(file.path(results_dir,'12_S2_Timing_Results.csv'))
timing <- timing[exposure == 'lifestyle' & trans %in% c(1,5,8,10) & status == 'FITTED']
stopifnot(nrow(timing) > 0L)
plot <- ggplot(timing, aes(HR, spec_id)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = 'y', width = 0.15) +
  geom_point() + geom_text(aes(label = star), vjust = -0.6, size = 3) +
  scale_x_log10() + facet_wrap(~trans, scales = 'free_y') +
  labs(x = 'Hazard ratio (95% CI)', y = 'Timing specification',
       caption = 'Non-estimable rows are recorded in the CSV and are not plotted as zero effects') + theme_classic()
ggsave(file.path(figure_dir,'Timing_Sensitivity.pdf'), plot, width = 11, height = 8)

# Supplementary S4: paired pointwise occupancy differences at 15% and full implementation.
smoking <- fread(file.path(results_dir, "13_S4_Probability_Differences.csv"))
stopifnot(setequal(unique(smoking$implementation_probability), c(.15, 1)))
smoking[, implementation := factor(implementation_probability, levels = c(.15, 1), labels = c("15%", "100%"))]
plot <- ggplot(smoking, aes(time, median, colour = implementation, fill = implementation)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .15, colour = NA) + geom_line() +
  facet_wrap(~state, scales = "free_y") +
  scale_colour_manual(values = c("15%" = "#D88933", "100%" = "#3978A8")) +
  scale_fill_manual(values = c("15%" = "#D88933", "100%" = "#3978A8")) +
  labs(x = "Years since baseline", y = "State occupancy difference (percentage points)",
    colour = "Implementation", fill = "Implementation", caption = "Pointwise intervals from 1000 paired predictions") + theme_classic()
ggsave(file.path(figure_dir, "S4_Smoking_15_vs_100_Probability.pdf"), plot, width = 11, height = 7)


# Figure 3E-G use the approved historical fixed-median-model individual analysis.
individual_files <- c("11b_Individual_Distribution_Bins.csv", "11b_Illustrative_Profile_Summary.csv",
  "11b_Illustrative_Profile_Exposures.csv", "11b_Illustrative_Profile_Contrasts.csv")
if (!all(file.exists(file.path(results_dir, individual_files)))) {
  stop("Run part 11b to create the approved individual-analysis outputs before plotting Figure 3E-G")
}
distribution <- fread(file.path(results_dir, "11b_Individual_Distribution_Bins.csv"))
example <- fread(file.path(results_dir, "11b_Illustrative_Profile_Summary.csv"))
stopifnot(nrow(example) == 1L)
selection_caption <- paste0("Eligibility risks: ", example$eligibility_factors,
  "; selection tier: ", example$selection_level, "; fallback used: ", example$fallback_used)

plot <- ggplot(distribution) +
  geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = count),
            fill = "#B5C9D6", colour = "white") +
  geom_vline(xintercept = example$healthy_S4_gain, colour = "#A94C3B", linewidth = .7) +
  labs(x = "Individual Healthy years gained under S4 versus S0", y = "Participants",
    subtitle = sprintf("Illustrative high-benefit profile: %.3f years; %.1fth empirical percentile",
                       example$healthy_S4_gain, example$percentile),
    caption = paste("Fixed median model; selected example is not representative of the cohort", selection_caption, sep = "\n")) + theme_classic()
ggsave(file.path(figure_dir, "Figure3E_Individual_Distribution.pdf"), plot, width = 10, height = 6)

profiles <- fread(file.path(results_dir, "11b_Illustrative_Profile_Exposures.csv"))
profiles[, scenario := factor(scenario, levels = paste0("S", 0:4))]
profiles[, risk_group := factor(higher_risk, levels = c(0, 1), labels = c("Lower risk", "Higher risk"))]
plot <- ggplot(profiles, aes(scenario, factor, fill = risk_group)) + geom_tile(colour = "white") +
  scale_fill_manual(values = c("Lower risk" = "#68A486", "Higher risk" = "#B65B50")) +
  labs(x = "Scenario", y = "Baseline risk factor", fill = "Profile",
       subtitle = "Illustrative profile; no participant identifier is displayed", caption = selection_caption) + theme_classic()
ggsave(file.path(figure_dir, "Figure3F_Illustrative_Exposures.pdf"), plot, width = 9, height = 5)

individual <- fread(file.path(results_dir, "11b_Illustrative_Profile_Contrasts.csv"))
individual <- individual[scenario != "S0"]
plot <- ggplot(individual, aes(time, difference_pp, colour = scenario)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line() + facet_wrap(~state, scales = "free_y") +
  labs(x = "Years since baseline", y = "Individual state occupancy difference (percentage points)",
       caption = paste("Deterministic fixed-median-model prediction; no individual uncertainty intervals", selection_caption, sep = "\n")) + theme_classic()
ggsave(file.path(figure_dir, "Figure3G_Illustrative_Contrasts.pdf"), plot, width = 11, height = 7)
cat("Saved Figure3A-G and supplementary result-based templates.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(file.path(results_dir, "figures"), pattern = "[.]pdf$", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^14_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
