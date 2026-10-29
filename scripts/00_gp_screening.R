# Complete GP preparation in one standalone script; no additional GP scripts are required.
# All three steps are implemented below: filter codes, assign dates, and retain the earliest date per participant.
# Packages and local paths ----
library(data.table)
source_dir <- Sys.getenv("UKB_SOURCE_DIR", unset = file.path(getwd(), "data", "raw"))
results_dir <- Sys.getenv("UKB_RESULTS_DIR", unset = file.path(getwd(), "results"))
gp_file <- Sys.getenv("UKB_GP_FILE", unset = file.path(source_dir, "GP_clinicalevents.csv"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
data.table::setDTthreads(4L)

# Author code lists ----
# The author-supplied Read v2 and CTV3 code lists are included directly below.
# Match Read v2 and CTV3 codes exactly and preserve their case.
bph_v2 <- c(
  "K20..", "K200.", "K201.", "K202.", "K20z."
)

bph_ctv3 <- c(
  "1AA..", "K20..", "K200.", "K201.", "K202.", "K20z.", "X30Nz",
  "X30O0", "X76ey", "X78YD", "XM1N5", "Xa7Xy", "XE0e6", "XE0e7"
)

depress_v2 <- c(
  "E11z2", "E2B..", "E2B1.", "Eu32.", "Eu320", "Eu321", "Eu322",
  "Eu323", "Eu324", "Eu325", "Eu326", "Eu327", "Eu328", "Eu329",
  "Eu32B", "Eu32y", "Eu32z", "Eu32A", "Eu33.", "Eu330", "Eu331",
  "Eu332", "Eu333", "Eu334", "Eu33y", "Eu33z"
)

depress_ctv3 <- c(
  "E2B..", "X40MH", "2257.", "E0043", "E1121", "E1122", "E1123",
  "E1124", "E1126", "E11y2", "E130.", "E2B0.", "E2B1.", "Eu320",
  "Eu321", "Eu322", "Eu323", "Eu32y", "Eu32z", "X00Qy", "X00SO",
  "X00SQ", "X00SU", "XE1Y0", "XE1YC", "XE1ZY", "XE1ZZ", "XE1Za",
  "XE1Zb", "XSEGJ", "XSGok", "XSGol", "XSGom", "XSGon", "XaB9J",
  "XaCHr", "XaCHs", "XaCIs", "XaCIt", "XaCIu", "XaX53", "XaY2C",
  "E1131", "E1132", "E1133", "E1134", "E1135", "E1136", "E1137",
  "Eu330", "Eu331", "Eu332", "Eu333", "Eu334", "Eu33y", "Eu33z",
  "X00SR", "X00SS", "XE1Y1", "XE1Zc", "XE1Zd", "XE1Ze", "XE1Zf",
  "XM1GC", "Xa0wV", "XaCHo", "XaX54"
)

pca_v2 <- c(
  "B46.."
)

pca_ctv3 <- c(
  "B46..", "X78Y6", "X78id", "Xa0bb", "XaFrM", "XaFrk"
)

# Stage 1: read the four relevant columns and keep target records. ----
# Reading only these columns uses less memory than reading the entire export.
columns <- c("Participant ID(participant - eid)", "Date clinical code was entered",
             "Read v2 code", "CTV3 (Read v3) code")
gp <- fread(gp_file, select = columns, colClasses = "character", na.strings = c("", "NA"))
for (column in columns) {
  value <- gsub("[\n\r\t]", "", gp[[column]])
  value <- gsub("\u3000", " ", value)
  set(gp, j = column, value = trimws(value))
}
baseline <- fread(file.path(source_dir, "base date + censor region.csv"),
                  select = "Participant ID", colClasses = "character")
stopifnot(!anyNA(baseline[["Participant ID"]]), !anyDuplicated(baseline[["Participant ID"]]))
gp <- gp[`Participant ID(participant - eid)` %in% baseline[["Participant ID"]]]
read_v2_keep <- unique(c(bph_v2, depress_v2, pca_v2))
ctv3_keep <- unique(c(bph_ctv3, depress_ctv3, pca_ctv3))
gp <- gp[`Read v2 code` %in% read_v2_keep | `CTV3 (Read v3) code` %in% ctv3_keep]
saveRDS(gp, file.path(results_dir, "00_GP_1_Filtered_Records.rds"))

# Stage 2: tag each target condition and assign its diagnosis date. ----
setnames(gp, "Participant ID(participant - eid)", "Participant ID")
gp[, BPH_match := `Read v2 code` %in% bph_v2 | `CTV3 (Read v3) code` %in% bph_ctv3]
gp[, depress_match := `Read v2 code` %in% depress_v2 | `CTV3 (Read v3) code` %in% depress_ctv3]
gp[, PCa_match := `Read v2 code` %in% pca_v2 | `CTV3 (Read v3) code` %in% pca_ctv3]
date_text <- gp[["Date clinical code was entered"]]
diagnosis_date <- as.IDate(date_text, format = "%Y-%m-%d")

# Keep date-quality flags required by the next preparation step.
# UKB birth-date masks retain pre-baseline evidence without an invented exact day.
baseline_evidence <- date_text %in% c("1902-02-02", "1903-03-03")
unknown_date <- is.na(diagnosis_date) | date_text %in% c("1901-01-01", "2037-07-07")
valid_date <- !baseline_evidence & !unknown_date
for (condition in c("BPH", "depress", "PCa")) {
  matched <- gp[[paste0(condition, "_match")]]
  value <- diagnosis_date
  value[!matched | !valid_date] <- as.IDate(NA_character_)
  set(gp, j = paste0(condition, "-date"), value = value)
  set(gp, j = paste0("GP_", condition, "_baseline"), value = as.integer(matched & baseline_evidence))
  set(gp, j = paste0("GP_", condition, "_unknown"), value = as.integer(matched & unknown_date))
}
date_columns <- c("BPH-date", "depress-date", "PCa-date")
flag_columns <- paste0("GP_", rep(c("BPH", "depress", "PCa"), each = 2), "_", rep(c("baseline", "unknown"), 3))
gp_dates <- gp[, c("Participant ID", date_columns, flag_columns), with = FALSE]
saveRDS(gp_dates, file.path(results_dir, "00_GP_2_Assigned_Dates.rds"))

# Stage 3: one row per participant. All-missing dates remain NA, not Inf. ----
GP_final <- gp_dates[, .(
  `BPH-date` = if (all(is.na(`BPH-date`))) as.IDate(NA_character_) else min(`BPH-date`, na.rm = TRUE),
  `depress-date` = if (all(is.na(`depress-date`))) as.IDate(NA_character_) else min(`depress-date`, na.rm = TRUE),
  `PCa-date` = if (all(is.na(`PCa-date`))) as.IDate(NA_character_) else min(`PCa-date`, na.rm = TRUE),
  GP_BPH_baseline = as.integer(any(GP_BPH_baseline == 1L)),
  GP_BPH_unknown = as.integer(any(GP_BPH_unknown == 1L)),
  GP_depress_baseline = as.integer(any(GP_depress_baseline == 1L)),
  GP_depress_unknown = as.integer(any(GP_depress_unknown == 1L)),
  GP_PCa_baseline = as.integer(any(GP_PCa_baseline == 1L)),
  GP_PCa_unknown = as.integer(any(GP_PCa_unknown == 1L))
), by = "Participant ID"]
stopifnot(!anyDuplicated(GP_final[["Participant ID"]]))
fwrite(GP_final, file.path(results_dir, "00_GP_Dates.csv"), na = "")
cat("GP screening complete. Participants with target records:", nrow(GP_final), "\n")
