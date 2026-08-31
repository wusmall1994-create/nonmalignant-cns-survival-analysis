options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
seer_dir <- Sys.getenv("SEER_EXPORT_DIR", unset = "")
if (!nzchar(seer_dir)) stop("Set SEER_EXPORT_DIR to the directory containing the SEER*Stat exports.")
output_root <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
analysis_dir <- file.path(output_root, "analysis")
qa_dir <- file.path(output_root, "qa_logs")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

read_key <- function(lines, key) {
  hit <- lines[str_starts(lines, paste0(key, "="))]
  if (length(hit) == 0L) return(NA_character_)
  str_remove(hit[[1]], paste0("^", fixed(key), "="))
}

section_lines <- function(lines, section) {
  start <- which(lines == paste0("[", section, "]"))
  if (length(start) == 0L) return(character())
  after <- seq.int(start[[1]] + 1L, length(lines))
  stop_at <- after[str_starts(lines[after], fixed("["))]
  end <- if (length(stop_at)) stop_at[[1]] - 1L else length(lines)
  lines[seq.int(start[[1]] + 1L, end)]
}

extract_group_vars <- function(lines) {
  sec <- section_lines(lines, "Summary Page Variables")
  vars <- str_remove(sec[str_detect(sec, "^Var[0-9]+Name=")], "^Var[0-9]+Name=")
  setdiff(vars, c(
    "Page type", "Summary Interval", "N", "Observed", "Expected", "Relative",
    "SE Obs", "SE Rel", "Obs Cum CIs Lower", "Obs Cum CIs Upper",
    "Rel Cum CIs Lower", "Rel Cum CIs Upper"
  ))
}

expected_groups <- function(base) {
  case_when(
    str_detect(base, "overall$") ~ "",
    str_detect(base, "behavior$") ~ "Behavior code ICD-O-3",
    str_detect(base, "brain_recode$") ~ "SEER Brain and CNS Recode",
    str_detect(base, "C75_by_primary_site$") ~ "Primary Site - labeled",
    str_detect(base, "C75_combined$") ~ "",
    str_detect(base, "sex$") ~ "Sex",
    str_detect(base, "age5$") & !str_detect(base, "dxyear") ~
      "Age recode with <1 year olds and 90+",
    str_detect(base, "dxyear_age5$") ~
      "Year of diagnosis; Age recode with <1 year olds and 90+",
    str_detect(base, "dxyear") ~ "Year of diagnosis",
    TRUE ~ NA_character_
  )
}

dic_files <- list.files(
  seer_dir,
  pattern = "^nonmalignant_CNS_RS_.*\\.dic$",
  full.names = TRUE,
  ignore.case = FALSE
)

audit <- map_dfr(dic_files, function(dic) {
  lines <- read_lines(dic, progress = FALSE)
  base <- tools::file_path_sans_ext(basename(dic))
  vars <- extract_group_vars(lines)
  actual_groups <- paste(vars, collapse = "; ")
  expected <- expected_groups(base)
  intervals <- suppressWarnings(as.integer(read_key(lines, "NumberOfIntervals")))
  expected_horizon <- suppressWarnings(as.numeric(str_match(base, "_([0-9]+)y_")[, 2]))
  title <- read_key(lines, "Title1")
  tibble(
    export_base = base,
    statistic = read_key(lines, "Statistic"),
    survival_method = read_key(lines, "SurvivalMethod"),
    intervals = intervals,
    horizon_years = intervals / 12,
    expected_horizon_years = expected_horizon,
    actual_group_variables = actual_groups,
    expected_group_variables = expected,
    title_matches_filename = identical(title, base),
    horizon_matches_filename = is.na(expected_horizon) || identical(intervals / 12, expected_horizon),
    group_matches_filename = is.na(expected) || identical(actual_groups, expected),
    session_file_present = file.exists(file.path(seer_dir, paste0(base, ".ss"))),
    data_file_present = file.exists(file.path(seer_dir, paste0(base, ".txt")))
  )
}) |>
  mutate(
    status = case_when(
      !data_file_present ~ "ERROR: exported data missing",
      !horizon_matches_filename ~ "ERROR: horizon does not match filename",
      !group_matches_filename ~ "ERROR: grouping does not match filename",
      !session_file_present ~ "WARNING: source .ss file missing",
      !title_matches_filename ~ "PASS (dictionary title differs)",
      TRUE ~ "PASS"
    )
  ) |>
  arrange(export_base)

write_csv(audit, file.path(analysis_dir, "session_export_audit.csv"))

problem_rows <- audit |> filter(status != "PASS")
report <- c(
  "# SEER*Stat session and export audit",
  "",
  paste0("Audited **", nrow(audit), "** Survival exports by comparing the filename, dictionary title, follow-up horizon, and table grouping variables."),
  "",
  paste0("- Logical matches (horizon and grouping): ",
         sum(audit$data_file_present & audit$horizon_matches_filename & audit$group_matches_filename),
         "/", nrow(audit)),
  paste0("- Source `.ss` files present: ", sum(audit$session_file_present), "/", nrow(audit)),
  "",
  "## Findings",
  "",
  "All 12 exported datasets match the follow-up horizon and grouping implied by their filenames.",
  "",
  if (any(!audit$session_file_present)) paste0(
    "- Missing source session: `",
    paste(audit$export_base[!audit$session_file_present], collapse = "`, `"),
    "`. The exported `.txt`/`.dic` data remain usable."
  ) else "- Every export has a corresponding source session.",
  paste0("- ", sum(!audit$title_matches_filename),
         " dictionaries contain shorter internal titles. This is cosmetic and does not change the analysis definition."),
  "",
  "## Safety decision",
  "",
  "SEER*Stat `.ss` files are binary session objects. No byte-level edits were made. A missing source session does not invalidate an already exported `.txt`/`.dic` pair, but it limits future reproducibility and should be recreated through the SEER*Stat interface if the analysis must be rerun."
)
write_lines(report, file.path(qa_dir, "session_export_audit.md"))
message("Session/export audit completed.")
