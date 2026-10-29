# Baseline covariates, first chained-equation completion and six risk indicators.
# Raw variables are imputed before derived risks. Entirely unavailable covariates stay unavailable.

# 1. Packages and local paths ----
library(data.table)
library(survival)
source_dir <- Sys.getenv('UKB_SOURCE_DIR', unset = file.path(getwd(), 'data', 'raw'))
results_dir <- Sys.getenv('UKB_RESULTS_DIR', unset = file.path(getwd(), 'results'))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, 'checks'), showWarnings = FALSE)
options(stringsAsFactors = FALSE)
data.table::setDTthreads(4L)

# Impute weekly bread/cereal amounts first, divide by seven, then classify daily grain components.
# Field 2139 restores explicit skip-rule zeros before imputation.
sexual_history_file <- Sys.getenv("UKB_SEX_HISTORY_FILE",
  unset = file.path(source_dir, "sex_participant.csv"))

# 2. Definitions used by this analysis ----
RISK_FACTORS <- c(
  "smoke", "sport", "sleep", "social_support_score", "new_food_score_score",
  "sexpartners_score"
)

DIET_COMPONENTS <- c(
  "fruit", "vegetable", "fish", "processed_meat", "red_meat", "whole_grains",
  "refined_grains"
)

SOCIAL_COMPONENTS <- c("visits", "confide", "social_activities")

less_than_one_value <- 0.5

RAW_BASELINE_FIELD_MAP <- c(
  Participant_ID = "Participant ID", Age_base = "Age at recruitment",
  Townsend = "Townsend deprivation index at recruitment", BMI = "Body mass index (BMI) | Instance 0",
  Race = "Ethnic background | Instance 0", Education = "Qualifications | Instance 0", EverSmoked = "Ever smoked | Instance 0",
  PhysicalActivity = "At or above moderate/vigorous/walking recommendation | Instance 0",
  SleepDuration = "Sleep duration | Instance 0", SexPartners = "Lifetime number of sexual partners | Instance 0",
  FreshFruit = "Fresh fruit intake | Instance 0", DriedFruit = "Dried fruit intake | Instance 0",
  CookedVeg = "Cooked vegetable intake | Instance 0", RawVeg = "Salad / raw vegetable intake | Instance 0",
  OilyFish = "Oily fish intake | Instance 0", NonOilyFish = "Non-oily fish intake | Instance 0",
  ProcessedMeat = "Processed meat intake | Instance 0", Beef = "Beef intake | Instance 0",
  Lamb = "Lamb/mutton intake | Instance 0", Pork = "Pork intake | Instance 0", BreadIntake = "Bread intake | Instance 0",
  BreadType = "Bread type | Instance 0", CerealIntake = "Cereal intake | Instance 0", CerealType = "Cereal type | Instance 0",
  FriendVisits = "Frequency of friend/family visits | Instance 0", Confide = "Able to confide | Instance 0",
  SocialActivities = "Leisure/social activities | Instance 0"
)

RAW_IMPUTATION_FIELDS <- c(
  "Townsend", "Age_base", "BMI", "Race", "Education", "EverSmoked", "PhysicalActivity",
  "SleepDuration", "SexPartners", "FreshFruit", "DriedFruit", "CookedVeg", "RawVeg", "OilyFish",
  "NonOilyFish", "ProcessedMeat", "Beef", "Lamb", "Pork", "BreadIntake", "BreadType", "CerealIntake",
  "CerealType", "FriendVisits", "Confide", "SocialActivities"
)

# 3. Reused algorithms needed repeatedly in this analysis ----
# Algorithm: map_raw_baseline ----
map_raw_baseline <- function(d, field_map) {
  stopifnot(all(unname(field_map) %in% names(d)))
  z <- d
  for (canonical in names(field_map)) z[[canonical]] <- d[[field_map[[canonical]]]]
  numeric_fields <- intersect(
    c(
      "Age_base", "Townsend", "BMI", "SleepDuration", "SexPartners",
      "FreshFruit", "DriedFruit", "CookedVeg", "RawVeg", "BreadIntake", "CerealIntake"
    ),
    names(field_map)
  )
  for (v in numeric_fields) z[[v]] <- {
    raw_value <- clean_response(z[[v]])
    number_value <- suppressWarnings(as.numeric(as.character(raw_value)))
    if (any(!is.na(raw_value) & !is.finite(number_value)))
      stop("Unexpected nonnumeric value")
    number_value
  }
  if ("SexPartners" %in% names(z))
    z$SexPartners <- pmin(z$SexPartners, 100)
  for (v in intersect(c("EverSmoked", "PhysicalActivity"), names(field_map))) {
    x <- clean_response(z[[v]])
    if (any(!is.na(x) & !x %in% c("No", "Yes")))
      stop("Invalid ", v)
    z[[v]] <- factor(x, levels = c("No", "Yes"))
  }
  if ("SocialActivities" %in% names(field_map)) {
    x <- clean_response(z$SocialActivities)
    z$SocialActivities <- factor(ifelse(is.na(x), NA, as.integer(x == "None of the above")),
      levels = 0:1
    )
  }
  freq <- c(
    "Never", "Less than once a week", "Once a week", "2-4 times a week", "5-6 times a week",
    "Once or more daily"
  )
  for (v in intersect(
    c("OilyFish", "NonOilyFish", "ProcessedMeat", "Beef", "Lamb", "Pork"),
    names(field_map)
  )) {
    x <- clean_response(z[[v]])
    intake_rank(x)
    z[[v]] <- ordered(x, levels = freq)
  }
  for (v in intersect(c("BreadType", "CerealType", "Education", "Race"), names(field_map))) z[[v]] <- factor(clean_response(z[[v]]))
  if ("FriendVisits" %in% names(field_map))
    z$FriendVisits <- ordered(clean_response(z$FriendVisits), levels = c(
      "Almost daily",
      "2-4 times a week", "About once a week", "About once a month", "Once every few months",
      "Never or almost never", "No friends/family outside household"
    ))
  if ("Confide" %in% names(field_map))
    z$Confide <- ordered(clean_response(z$Confide), levels = c(
      "Almost daily", "2-4 times a week",
      "About once a week", "About once a month", "Once every few months", "Never or almost never"
    ))
  for (v in intersect(c("FriendVisits", "Confide"), names(field_map))) if (any(!is.na(clean_response(d[[field_map[[v]]]])) &
    is.na(z[[v]])))
    stop("Unknown categorical level ", v)
  z
}

# Algorithm: encode_lifestyle_components ----
encode_lifestyle_components <- function(d) {
  fields <- c(
    "EverSmoked", "PhysicalActivity", "SleepDuration", "FriendVisits", "Confide",
    "SocialActivities", "FreshFruit", "DriedFruit", "CookedVeg", "RawVeg", "OilyFish",
    "NonOilyFish", "ProcessedMeat", "Beef", "Lamb", "Pork", "BreadIntake", "BreadType",
    "CerealIntake", "CerealType", "SexPartners"
  )
  stopifnot(all(fields %in% names(d)))
  z <- d
  z$smoke <- encode_yes_no(d$EverSmoked)
  z$sport <- encode_yes_no(d$PhysicalActivity, yes = 0, no = 1)
  z$sleep <- encode_sleep(d$SleepDuration)
  visits <- clean_response(d$FriendVisits)
  confide <- clean_response(d$Confide)
  allowed_visits <- c(
    "Almost daily", "2-4 times a week", "About once a week", "About once a month",
    "Once every few months", "Never or almost never", "No friends/family outside household"
  )
  allowed_confide <- c(
    "Almost daily", "2-4 times a week", "About once a week", "About once a month",
    "Once every few months", "Never or almost never"
  )
  if (any(!is.na(visits) & !visits %in% allowed_visits) || any(!is.na(confide) & !confide %in%
    allowed_confide))
    stop("Unexpected social support response")
  z$visits <- ifelse(is.na(visits), NA, ifelse(visits %in% allowed_visits[1:2], 0, ifelse(visits ==
    "About once a week", 1, 2)))
  z$confide <- ifelse(is.na(confide), NA, ifelse(confide == "Almost daily", 0, ifelse(confide %in%
    allowed_confide[2:4], 1, 2)))
  social <- clean_response(d$SocialActivities)
  if (any(!is.na(social) & !social %in% c("0", "1")))
    stop("SocialActivities must be the canonical 0/1 participation-deficit factor; use map_raw_baseline() for raw leisure text")
  z$social_activities <- as.numeric(social)
  z$sexpartners <- {
    raw_value <- clean_response(d$SexPartners)
    number_value <- suppressWarnings(as.numeric(as.character(raw_value)))
    if (any(!is.na(raw_value) & !is.finite(number_value)))
      stop("Unexpected nonnumeric value")
    number_value
  }
  if (any(z$sexpartners < 0, na.rm = TRUE))
    stop("Negative lifetime partner count")
  z$sexpartners <- pmin(z$sexpartners, 100)
  n <- function(v) {
    x <- {
      raw_value <- clean_response(d[[v]])
      number_value <- suppressWarnings(as.numeric(as.character(raw_value)))
      if (any(!is.na(raw_value) & !is.finite(number_value)))
        stop("Unexpected nonnumeric value")
      number_value
    }
    if (any(x < 0, na.rm = TRUE))
      stop("Negative intake needs an explicit source-code mapping: ", v)
    x
  }
  z$fruit <- as.integer(n("FreshFruit") + n("DriedFruit") >= 3)
  z$vegetable <- as.integer(n("CookedVeg") + n("RawVeg") >= 3)
  of <- intake_rank(d$OilyFish)
  nf <- intake_rank(d$NonOilyFish)
  z$fish <- ifelse(is.na(of) | is.na(nf), NA_integer_, as.integer(of >= 4 | nf >= 4 | (of ==
    3 & nf == 3)))
  pm <- intake_rank(d$ProcessedMeat)
  z$processed_meat <- ifelse(is.na(pm), NA_integer_, as.integer(pm <= 3))
  meat <- cbind(intake_rank(d$Beef), intake_rank(d$Lamb), intake_rank(d$Pork))
  z$red_meat <- apply(meat, 1, function(x) if (anyNA(x))
    NA_integer_
  else as.integer(all(x <= 2) || (sum(x == 3) == 1 && sum(x == 1) == 2)))
  bt <- clean_response(d$BreadType)
  ct <- clean_response(d$CerealType)
  bread_types <- c("Wholemeal or wholegrain", "White", "Brown", "Other type of bread")
  cereal_types <- c(
    "Biscuit cereal (e.g. Weetabix)", "Oat cereal (e.g. Ready Brek, porridge)",
    "Muesli", "Bran cereal (e.g. All Bran, Branflakes)", "Other (e.g. Cornflakes, Frosties)"
  )
  if (any(!is.na(bt) & !bt %in% bread_types) || any(!is.na(ct) & !ct %in% cereal_types))
    stop("Unexpected bread/cereal type")
  wb <- ifelse(is.na(bt), NA_real_, ifelse(bt == bread_types[1], n("BreadIntake"), 0))
  wc <- ifelse(is.na(ct), NA_real_, ifelse(ct %in% cereal_types[1:3], n("CerealIntake"),
    0
  ))
  rb <- ifelse(is.na(bt), NA_real_, ifelse(bt %in% bread_types[2:4], n("BreadIntake"), 0))
  rc <- ifelse(is.na(ct), NA_real_, ifelse(ct %in% cereal_types[4:5], n("CerealIntake"),
    0
  ))
  # Convert the original weekly bread/cereal quantities to daily amounts after imputation.
  whole_grains_daily <- (wb + wc) / 7
  refined_grains_daily <- (rb + rc) / 7
  z$whole_grains <- as.integer(whole_grains_daily >= 3)
  z$refined_grains <- as.integer(refined_grains_daily <= 1.5)
  z$food_score <- rowSums(z[, DIET_COMPONENTS, drop = FALSE])
  derive_risk_scores(z)
}

# Algorithm: recode_demographics ----
recode_demographics <- function(d) {
  stopifnot(all(c("Race", "Education") %in% names(d)))
  x <- as.character(d$Race)
  race <- rep(NA_character_, length(x))
  other <- !is.na(x) & grepl("Mixed|Other ethnic group|White and (Black|Asian)", x, ignore.case = TRUE)
  black <- !is.na(x) & !other & grepl("Black|Caribbean|African", x, ignore.case = TRUE)
  asian <- !is.na(x) & !other & !black & grepl("Asian|Indian|Pakistani|Bangladeshi|Chinese",
    x,
    ignore.case = TRUE
  )
  white <- !is.na(x) & !other & !black & !asian & grepl("White|British|Irish", x, ignore.case = TRUE)
  race[other] <- "Other"
  race[black] <- "Black"
  race[asian] <- "Asian"
  race[white] <- "White"
  y <- as.character(d$Education)
  edu <- rep(NA_character_, length(y))
  degree <- !is.na(y) & grepl("College|University degree", y, ignore.case = TRUE)
  vocational <- !is.na(y) & !degree & grepl("NVQ|HND|HNC|Other professional qualification|nursing|teaching",
    y,
    ignore.case = TRUE
  )
  school <- !is.na(y) & !degree & !vocational & grepl("A levels|AS levels|O levels|GCSEs?|CSEs?",
    y,
    ignore.case = TRUE
  )
  none <- !is.na(y) & !degree & !vocational & !school & grepl("None of the above", y, ignore.case = TRUE)
  edu[degree] <- "Degree"
  edu[vocational] <- "Vocational"
  edu[school] <- "School"
  edu[none] <- "None"
  if (any(!is.na(x) & is.na(race)) || any(!is.na(y) & is.na(edu)))
    stop("Unmapped ethnicity/education label; extend the explicit dictionary rather than guess")
  d$Race <- factor(race, levels = c("White", "Black", "Asian", "Other"))
  d$Education <- factor(edu, levels = c("Degree", "Vocational", "School", "None"))
  d
}

# Algorithm: clean_response ----
clean_response <- function(x) {
  if (is.factor(x))
    x <- as.character(x)
  if (is.numeric(x))
    x[x %in% c(-1, -3)] <- NA_real_
  if (is.character(x)) {
    x <- trimws(x)
    x[x %in% c("", "NA", "Prefer not to answer", "Do not know", "-1", "-3")] <- NA_character_
  }
  x
}

# Algorithm: intake_rank ----
intake_rank <- function(x) {
  allowed <- c(
    "Never", "Less than once a week", "Once a week", "2-4 times a week", "5-6 times a week",
    "Once or more daily"
  )
  x <- clean_response(x)
  if (any(!is.na(x) & !x %in% allowed))
    stop("Unexpected food frequency response")
  match(x, allowed)
}

# Algorithm: encode_yes_no ----
encode_yes_no <- function(x, yes = 1, no = 0) {
  x <- clean_response(x)
  if (any(!is.na(x) & !x %in% c("Yes", "No")))
    stop("Expected Yes/No response")
  ifelse(is.na(x), NA_real_, ifelse(x == "Yes", yes, no))
}

# Algorithm: encode_sleep ----
encode_sleep <- function(x) {
  x <- clean_response(x)
  n <- {
    raw_value <- x
    number_value <- suppressWarnings(as.numeric(as.character(raw_value)))
    if (any(!is.na(raw_value) & !is.finite(number_value)))
      stop("Unexpected nonnumeric value")
    number_value
  }
  if (any(!is.na(n) & (n < 0 | n > 24)))
    stop("Sleep hours outside 0-24")
  ifelse(is.na(n), NA_real_, as.integer(n < 6 | n > 8))
}

# Algorithm: derive_risk_scores ----
derive_risk_scores <- function(d) {
  vars <- c(
    "smoke", "sport", "sleep", "vegetable", "processed_meat", "red_meat", SOCIAL_COMPONENTS,
    "sexpartners"
  )
  stopifnot(all(vars %in% names(d)))
  for (v in vars) d[[v]] <- {
    raw_value <- d[[v]]
    number_value <- suppressWarnings(as.numeric(as.character(raw_value)))
    if (any(!is.na(raw_value) & !is.finite(number_value)))
      stop("Unexpected nonnumeric value")
    number_value
  }
  for (v in setdiff(vars, c("visits", "confide", "sexpartners"))) if (any(!is.na(d[[v]]) &
    !d[[v]] %in% c(0, 1)))
    stop("Non-binary component: ", v)
  for (v in c("visits", "confide")) if (any(!is.na(d[[v]]) & !d[[v]] %in% 0:2))
    stop("Invalid social component: ", v)
  d$new_food_score <- d$vegetable + d$processed_meat + d$red_meat
  d$social_support <- d$visits / 2 + d$confide / 2 + d$social_activities
  d$new_food_score_score <- ifelse(is.na(d$new_food_score), NA_integer_, as.integer(d$new_food_score <=
    1))
  d$social_support_score <- ifelse(is.na(d$social_support), NA_integer_, as.integer(d$social_support >=
    1.5))
  d$sexpartners_score <- ifelse(is.na(d$sexpartners), NA_integer_, as.integer(d$sexpartners >=
    5))
  d$lifestyle <- rowSums(d[, RISK_FACTORS, drop = FALSE])
  d
}

# 4. Read inputs, perform this analysis, and save inspectable results ----
# Check saved inputs ----
# Verify private input/output provenance before reading analysis files.
# Run this part from the repository root. Hashes detect changed files; they are not public data.
stage_name <- "06"
stage_dependencies <- c("05")
stage_code_file <- normalizePath(file.path(getwd(), "scripts", "06_covariates_imputation.R"), mustWork = TRUE)
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
stage_input_files <- c(file.path(results_dir, c("05_All_Outcome_Dates.rds")),
                         file.path(source_dir, c("new-cov_participant.csv")))
optional_input_files <- c(sexual_history_file, file.path(source_dir, "additional_covariates.csv"))
stage_input_files <- c(stage_input_files, optional_input_files[file.exists(optional_input_files)])
stage_absent_input_files <- optional_input_files[!file.exists(optional_input_files)]
stopifnot(all(file.exists(stage_input_files)))
stage_input_files <- normalizePath(stage_input_files, mustWork = TRUE)
stage_input_hashes <- tools::md5sum(stage_input_files)

outcomes <- readRDS(file.path(results_dir, "05_All_Outcome_Dates.rds"))

stopifnot("pre_tie_eligible" %in% names(outcomes))
pre_tie <- outcomes[pre_tie_eligible == TRUE]

covariate_file <- file.path(source_dir, "new-cov_participant.csv")

covariate_header <- names(fread(covariate_file, nrows = 0L))

sex_history_column <- "Age first had sexual intercourse | Instance 0"
method_columns <- intersect(c(unname(RAW_BASELINE_FIELD_MAP), sex_history_column), covariate_header)

raw <- fread(covariate_file, select = method_columns, colClasses = "character", na.strings = c(
  "",
  "NA"
))

stopifnot(!anyDuplicated(raw[["Participant ID"]]))

raw <- raw[match(pre_tie$Participant_ID, raw[["Participant ID"]])]

stopifnot(identical(raw[["Participant ID"]], pre_tie$Participant_ID))

extra_file <- file.path(source_dir, "additional_covariates.csv")

if (file.exists(extra_file)) {
  extra <- fread(extra_file, colClasses = "character", na.strings = c("", "NA"))
  stopifnot(!anyDuplicated(extra[["Participant ID"]]))
  extra <- extra[match(raw[["Participant ID"]], extra[["Participant ID"]])]
  stopifnot(identical(extra[["Participant ID"]], raw[["Participant ID"]]))
  extra_columns <- setdiff(names(extra), names(raw))
  raw <- cbind(raw, extra[, ..extra_columns])
}

intake_columns <- unname(RAW_BASELINE_FIELD_MAP[c(
  "CookedVeg", "RawVeg", "FreshFruit", "DriedFruit",
  "BreadIntake", "CerealIntake"
)])

for (column in intake_columns) {
  rows <- which(raw[[column]] == "Less than one")
  set(raw, i = rows, j = column, value = as.character(less_than_one_value))
}

# Field 2139 is used only to restore the skip-rule zero, not as a model covariate.
if (file.exists(sexual_history_file)) {
  partner_column <- RAW_BASELINE_FIELD_MAP[["SexPartners"]]
  history_columns <- names(fread(sexual_history_file, nrows = 0L))
  stopifnot(all(c("Participant ID", sex_history_column) %in% history_columns))
  history_columns <- intersect(c("Participant ID", sex_history_column, partner_column), history_columns)
  sex_history <- fread(sexual_history_file, select = history_columns,
    colClasses = "character", na.strings = c("", "NA"))
  stopifnot(!anyDuplicated(sex_history[["Participant ID"]]))
  sex_history <- sex_history[match(raw[["Participant ID"]], sex_history[["Participant ID"]])]
  stopifnot(identical(sex_history[["Participant ID"]], raw[["Participant ID"]]))
  set(raw, j = sex_history_column, value = sex_history[[sex_history_column]])
  # Keep both sexual-history responses from the same supplied export when available.
  if (partner_column %in% names(sex_history)) {
    set(raw, j = partner_column, value = sex_history[[partner_column]])
  }
}
if (!sex_history_column %in% names(raw)) {
  stop("Provide baseline Field 2139 through the baseline export or UKB_SEX_HISTORY_FILE")
}
structural_zero_n <- 0L
if (sex_history_column %in% names(raw)) {
  never_sex <- trimws(raw[[sex_history_column]]) %in% c("Never had sex", "-2")
  partner_column <- RAW_BASELINE_FIELD_MAP[["SexPartners"]]
  partner_text <- trimws(raw[[partner_column]])
  partner_text[partner_text %in% c("", "NA", "Do not know", "Prefer not to answer", "-1", "-3")] <- NA_character_
  partner_value <- suppressWarnings(as.numeric(partner_text))
  stopifnot(!any(!is.na(partner_text) & is.na(partner_value)))
  if (any(never_sex & !is.na(partner_value) & partner_value > 0)) {
    stop("Field 2139 conflicts with a positive lifetime partner count; resolve the source conflict")
  }
  structural_zero <- never_sex & is.na(partner_value)
  structural_zero_n <- sum(structural_zero)
  set(raw, i = which(structural_zero), j = partner_column, value = "0")
}
fwrite(data.table(structural_zeros_restored = structural_zero_n,
  field_2139_available = sex_history_column %in% names(raw),
  skip_rule = "Explicit Never had sex plus missing partner count maps to zero"),
  file.path(results_dir, "checks", "06_Partner_Skip_Rule.csv"))

available_map <- RAW_BASELINE_FIELD_MAP[RAW_BASELINE_FIELD_MAP %in% names(raw)]

missing_fields <- setdiff(names(RAW_BASELINE_FIELD_MAP), names(available_map))

stopifnot(all(missing_fields %in% c("BMI", "Race", "Education")))

baseline <- map_raw_baseline(as.data.frame(raw), available_map)

imputation_columns <- intersect(RAW_IMPUTATION_FIELDS, names(baseline))

not_all_missing <- vapply(
  baseline[, imputation_columns, drop = FALSE], function(x) !all(is.na(x)),
  logical(1)
)

imputation_columns <- imputation_columns[not_all_missing]

imputation_input <- baseline[, imputation_columns, drop = FALSE]

stopifnot("Age_base" %in% names(imputation_input), !anyNA(imputation_input$Age_base))
method <- mice::make.method(imputation_input)
method["Age_base"] <- ""

predictor_matrix <- mice::make.predictorMatrix(imputation_input)
predictor_matrix["Age_base", ] <- 0L
fwrite(data.table(variable = names(imputation_input), missing_n = colSums(is.na(imputation_input))),
  file.path(results_dir, "checks", "06_Pre_Imputation_Missingness.csv"))

# Generate five completions but analyse completion 1, as confirmed by the author.
imputed <- mice::mice(imputation_input,
  m = 5L, maxit = 5L, seed = 12345L, method = method,
  predictorMatrix = predictor_matrix, printFlag = FALSE
)

baseline[, imputation_columns] <- mice::complete(imputed, action = 1L)
saveRDS(list(method = imputed$method, predictor_matrix = imputed$predictorMatrix,
  variables = names(imputation_input), m = imputed$m, maxit = 5L, seed = 12345L,
  RNGkind = RNGkind(), analysed_completion = 1L, pooling = FALSE,
  logged_events = imputed$loggedEvents), file.path(results_dir, "06_Imputation_Settings.rds"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "checks", "06_SessionInfo.txt"))

# Derive risk indicators only after baseline imputation.
scored <- encode_lifestyle_components(baseline)

if (!"BMI" %in% names(scored)) scored$BMI <- NA_real_

if (!"Race" %in% names(scored)) scored$Race <- NA_character_

if (!"Education" %in% names(scored)) scored$Education <- NA_character_

if (any(!is.na(scored$Race)) || any(!is.na(scored$Education))) {
  scored <- recode_demographics(scored)
}

selected_columns <- unique(c(
  "Participant_ID", "Townsend", "Age_base", "BMI", "Race", "Education",
  "sexpartners", "new_food_score", "social_support", "lifestyle", RISK_FACTORS, DIET_COMPONENTS,
  SOCIAL_COMPONENTS
))

scored <- scored[, selected_columns, drop = FALSE]

ready <- merge(as.data.frame(pre_tie), scored, by = "Participant_ID", all.x = TRUE, sort = FALSE)

for (disease in c("BPH", "depress", "PCa")) {
  ready[[paste0("self_", disease, "_miss")]] <- as.integer(ready[[paste0(
    "self_", disease,
    "_unknown"
  )]])
}

saveRDS(ready, file.path(results_dir, "06_Pre_Tie_Scored.rds"))

analysis_data <- ready[ready$exclusion_reason == "", , drop = FALSE]

saveRDS(analysis_data, file.path(results_dir, "06_Analysis_Data.rds"))

saveRDS(imputed, file.path(results_dir, "06_Private_MICE_Object.rds"))

fwrite(data.table(field = c("BMI", "Race", "Education"), available = vapply(analysis_data[
  ,
  c("BMI", "Race", "Education")
], function(x) any(!is.na(x)), logical(1))), file.path(
  results_dir,
  "checks", "06_Table1_Covariate_Availability.csv"
))

cat("Step 06 complete. Completion 1 used. Table1 availability is recorded separately.\n")

# Save this stage record ----
# Register only after this part completes; old downstream records are checked on their next use.
stopifnot(identical(stage_input_hashes, tools::md5sum(stage_input_files)),
          identical(stage_code_hash, tools::md5sum(stage_code_file)),
          !any(file.exists(stage_absent_input_files)))
stage_output_files <- c(list.files(results_dir, pattern = "^06_", full.names = TRUE),
  list.files(file.path(results_dir, "checks"), pattern = "^06_", full.names = TRUE))
stage_output_files <- stage_output_files[!dir.exists(stage_output_files)]
stopifnot(length(stage_output_files) > 0L)
stage_output_files <- normalizePath(stage_output_files, mustWork = TRUE)
manifest[[stage_name]] <- list(code = stage_code_hash, inputs = stage_input_hashes,
  outputs = tools::md5sum(stage_output_files), dependencies = stage_dependencies,
  absent_inputs = stage_absent_input_files,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
saveRDS(manifest, manifest_file)
