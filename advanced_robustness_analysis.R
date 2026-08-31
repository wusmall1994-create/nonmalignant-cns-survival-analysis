options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(splines)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
seer_dir <- Sys.getenv("SEER_EXPORT_DIR", unset = "")
if (!nzchar(seer_dir)) stop("Set SEER_EXPORT_DIR to the directory containing the SEER*Stat exports.")
case_file <- file.path(seer_dir, "nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt")
if (!file.exists(case_file)) stop("Case-listing file not found: ", case_file)

output_root <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
analysis_dir <- file.path(output_root, "analysis")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
invisible(lapply(c(analysis_dir, figure_dir, table_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

cols_needed <- c(
  "Patient ID", "Record number recode", "Year of diagnosis",
  "Age recode with single ages and 90+", "Sex",
  "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)",
  "Behavior code ICD-O-3", "Histologic Type ICD-O-3", "ICD-O-3 Hist/behav",
  "SEER Brain and CNS Recode", "Diagnostic Confirmation", "Type of Reporting Source",
  "Survival Days", "Vital status recode (study cutoff used)",
  "RX Summ--Surg Prim Site (2003-2022)", "RX Summ--Surg Prim Site 2023 (2023+)",
  "Chemotherapy recode (yes, no/unk) (2004+)",
  "Median household income inflation adj to 2024"
)

message("Reading individual-level data for advanced robustness analyses...")
raw <- read_tsv(
  case_file,
  col_select = all_of(cols_needed),
  col_types = cols(.default = col_character()),
  na = c("NA", ""),
  progress = FALSE,
  show_col_types = FALSE
)

collapse_brain_group <- function(x) {
  case_when(
    x == "2.1 Meningiomas" ~ "Meningiomas",
    x == "2.2 Tumors of the cranial and paraspinal nerves" ~ "Cranial/paraspinal nerve tumors",
    x == "2.9 Other non-malignant Brain/ONS" ~ "Other non-malignant brain/ONS",
    x == "5. All Other Cancer Types" ~ "Other included sites/types",
    TRUE ~ "Other CNS histologies"
  )
}

collapse_confirmation <- function(x) {
  case_when(
    grepl("histology|cytology|microscopic|immunophenotyping|laboratory test", x,
          ignore.case = TRUE) ~ "Microscopic/laboratory",
    grepl("Radiography|Direct visualization", x, ignore.case = TRUE) ~
      "Radiographic/direct visualization",
    TRUE ~ "Clinical/unknown"
  )
}

classify_surgery <- function(old_code, new_code, year) {
  code <- ifelse(year <= 2022L, old_code, new_code)
  numeric_code <- suppressWarnings(as.integer(sub("^A", "", code)))
  case_when(
    code %in% c("00", "A000") ~ "No surgery",
    !is.na(numeric_code) & numeric_code >= 10L & numeric_code <= 900L &
      !code %in% c("98", "99", "A980", "A990") ~ "Surgery performed",
    TRUE ~ "Unknown/not recorded"
  )
}

classify_income <- function(x) {
  lower_bound <- suppressWarnings(parse_number(x))
  case_when(
    is.na(x) | grepl("Unknown|missing|no match", x, ignore.case = TRUE) ~ "Unknown",
    grepl("^<", x) | (!is.na(lower_bound) & lower_bound < 60000) ~ "<$60,000",
    lower_bound < 75000 ~ "$60,000-$74,999",
    lower_bound < 90000 ~ "$75,000-$89,999",
    lower_bound < 110000 ~ "$90,000-$109,999",
    lower_bound >= 110000 ~ "$110,000+",
    TRUE ~ "Unknown"
  )
}

dat <- raw |>
  transmute(
    patient_id = `Patient ID`,
    record_number = `Record number recode`,
    year = suppressWarnings(as.integer(`Year of diagnosis`)),
    age = parse_number(`Age recode with single ages and 90+`),
    sex = trimws(Sex),
    race = `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`,
    behavior = `Behavior code ICD-O-3`,
    histology_code = `Histologic Type ICD-O-3`,
    histology_label = `ICD-O-3 Hist/behav`,
    brain_group = collapse_brain_group(`SEER Brain and CNS Recode`),
    confirmation = collapse_confirmation(`Diagnostic Confirmation`),
    reporting_source = `Type of Reporting Source`,
    survival_days = suppressWarnings(as.numeric(`Survival Days`)),
    vital_status = `Vital status recode (study cutoff used)`,
    surgery = classify_surgery(
      `RX Summ--Surg Prim Site (2003-2022)`,
      `RX Summ--Surg Prim Site 2023 (2023+)`,
      year
    ),
    chemotherapy = `Chemotherapy recode (yes, no/unk) (2004+)`,
    income_group = classify_income(`Median household income inflation adj to 2024`)
  ) |>
  mutate(
    event = as.integer(vital_status == "Dead"),
    time_days = ifelse(is.na(survival_days), NA_real_, pmax(survival_days, 0.5)),
    dco_autopsy = grepl("Death certificate|Autopsy", reporting_source, ignore.case = TRUE),
    alive_zero_time = event == 0L & (is.na(survival_days) | survival_days <= 0),
    age_group = case_when(
      age < 18 ~ "0-17",
      age < 40 ~ "18-39",
      age < 60 ~ "40-59",
      age < 75 ~ "60-74",
      TRUE ~ "75+"
    ),
    period = ifelse(year <= 2008L, "2004-2008", "2009-2013")
  )

eligible <- dat |>
  filter(
    year >= 2004L,
    year <= 2013L,
    !dco_autopsy,
    !alive_zero_time,
    !is.na(patient_id),
    !is.na(time_days),
    !is.na(age),
    sex %in% c("Male", "Female")
  ) |>
  arrange(patient_id, year, record_number) |>
  group_by(patient_id) |>
  slice_head(n = 1L) |>
  ungroup() |>
  mutate(
    followup_months = time_days / 30.4375,
    time10 = pmin(followup_months, 120),
    event10 = as.integer(event == 1L & followup_months <= 120),
    age_group = factor(age_group, levels = c("0-17", "18-39", "40-59", "60-74", "75+"))
  )

age_weights <- eligible |>
  count(age_group, name = "n") |>
  mutate(weight = n / sum(n)) |>
  select(age_group, weight)
write_csv(age_weights, file.path(analysis_dir, "direct_standardization_age_weights.csv"))

km_cell <- function(df) {
  fit <- survfit(Surv(time10, event10) ~ 1, data = df)
  sm <- summary(fit, times = 120, extend = TRUE)
  tibble(
    n = nrow(df),
    events = sum(df$event10),
    survival_10y = unname(sm$surv[1]),
    se_10y = unname(sm$std.err[1])
  )
}

direct_standardize <- function(df, group_var, domain_label) {
  group_sym <- rlang::sym(group_var)
  cells <- df |>
    filter(!is.na(!!group_sym), !!group_sym != "Unknown") |>
    group_by(!!group_sym, age_group) |>
    group_modify(~ km_cell(.x)) |>
    ungroup() |>
    left_join(age_weights, by = "age_group")
  valid_groups <- cells |>
    count(!!group_sym, name = "age_strata") |>
    filter(age_strata == nrow(age_weights)) |>
    pull(!!group_sym)
  cells |>
    filter((!!group_sym) %in% valid_groups) |>
    group_by(!!group_sym) |>
    summarise(
      n = sum(n),
      events = sum(events),
      standardized_os_10y = sum(weight * survival_10y),
      se = sqrt(sum((weight * se_10y)^2)),
      .groups = "drop"
    ) |>
    mutate(
      ci_lower = pmax(0, standardized_os_10y - 1.959963984540054 * se),
      ci_upper = pmin(1, standardized_os_10y + 1.959963984540054 * se),
      domain = domain_label,
      group = as.character(!!group_sym)
    ) |>
    select(domain, group, n, events, standardized_os_10y, se, ci_lower, ci_upper)
}

race_std <- direct_standardize(
  eligible |> filter(race != "Non-Hispanic Unknown Race"),
  "race",
  "Race/ethnicity"
)
income_std <- direct_standardize(
  eligible |> filter(income_group != "Unknown"),
  "income_group",
  "County income"
)
tumor_std <- direct_standardize(eligible, "brain_group", "Tumor group")
period_std <- direct_standardize(eligible, "period", "Diagnosis period")
standardized_os <- bind_rows(race_std, income_std, tumor_std, period_std)
write_csv(standardized_os, file.path(analysis_dir, "age_standardized_10y_os.csv"))

# Decompose the heterogeneous broad 'other non-malignant brain/ONS' group.
other_histology <- eligible |>
  filter(brain_group == "Other non-malignant brain/ONS") |>
  mutate(
    histology = ifelse(is.na(histology_label), paste0("ICD-O-3 ", histology_code), histology_label),
    histology = sub("^[0-9]{4}/[0-9]:\\s*", "", histology)
  ) |>
  group_by(histology_code, histology) |>
  group_modify(~ km_cell(.x)) |>
  ungroup() |>
  arrange(desc(n)) |>
  mutate(stable_for_reporting = n >= 200 & events >= 20)
write_csv(other_histology, file.path(analysis_dir, "other_group_histology_decomposition.csv"))

model_data <- eligible |>
  filter(
    race != "Non-Hispanic Unknown Race",
    income_group != "Unknown",
    surgery != "Unknown/not recorded"
  ) |>
  mutate(
    race = factor(
      race,
      levels = c(
        "Non-Hispanic White", "Non-Hispanic Black", "Hispanic (All Races)",
        "Non-Hispanic Asian or Pacific Islander",
        "Non-Hispanic American Indian/Alaska Native"
      )
    ),
    income_group = factor(
      income_group,
      levels = c(
        "<$60,000", "$60,000-$74,999", "$75,000-$89,999",
        "$90,000-$109,999", "$110,000+"
      )
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    behavior = factor(behavior, levels = c("Benign", "Borderline malignancy")),
    brain_group = factor(brain_group),
    surgery = factor(surgery),
    confirmation = factor(confirmation)
  ) |>
  drop_na(race, income_group, sex, behavior, brain_group, surgery, confirmation)

fit_period <- function(df, label) {
  fit <- coxph(
    Surv(time10, event10) ~ race + income_group + ns(age, df = 4) + sex + behavior +
      brain_group + surgery + strata(confirmation),
    data = df,
    ties = "efron"
  )
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      period = label,
      n = nrow(df),
      events = sum(df$event10),
      term,
      hazard_ratio = estimate,
      ci_lower = conf.low,
      ci_upper = conf.high,
      p_value = p.value
    )
}

period_models <- bind_rows(
  fit_period(filter(model_data, period == "2004-2008"), "2004-2008"),
  fit_period(filter(model_data, period == "2009-2013"), "2009-2013")
)
write_csv(period_models, file.path(analysis_dir, "period_specific_adjusted_associations.csv"))

# Figure 7 source data from prespecified SEER*Stat survival tables.
behavior_rlel <- read_csv(
  file.path(analysis_dir, "behavior_10y_rlel.csv"), show_col_types = FALSE
) |>
  filter(tau_months == 120) |>
  transmute(
    domain = "Behavior",
    group = `Behavior code ICD-O-3_label`,
    n, rlel_months, rlel_ci_lower, rlel_ci_upper
  )
sex_rlel <- read_csv(file.path(analysis_dir, "sex_10y_rlel.csv"), show_col_types = FALSE) |>
  filter(tau_months == 120, Sex_label %in% c("Female", "Male")) |>
  transmute(domain = "Sex", group = Sex_label, n, rlel_months, rlel_ci_lower, rlel_ci_upper)
brain_rlel <- read_csv(
  file.path(analysis_dir, "brain_recode_10y_rlel.csv"), show_col_types = FALSE
) |>
  filter(
    tau_months == 120,
    n >= 500,
    trimws(`SEER Brain and CNS Recode_label`) != "2. Non-malignant Brain/ONS"
  ) |>
  transmute(
    domain = "Tumor group",
    group = recode(
      sub("^2\\.[0-9]+\\s+", "", trimws(`SEER Brain and CNS Recode_label`)),
      "5. All Other Cancer Types" = "Other included sites/types"
    ),
    n, rlel_months, rlel_ci_lower, rlel_ci_upper
  ) |>
  arrange(desc(rlel_months)) |>
  slice_head(n = 6)
rlel_forest <- bind_rows(behavior_rlel, sex_rlel, brain_rlel)
write_csv(rlel_forest, file.path(analysis_dir, "figure7_rlel_forest_source.csv"))

age_relative <- read_csv(file.path(analysis_dir, "age_10y_rlel.csv"), show_col_types = FALSE) |>
  filter(tau_months == 120) |>
  mutate(
    age_label = `Age recode with <1 year olds and 90+_label`,
    age_numeric = ifelse(grepl("<1", age_label), 0, parse_number(age_label))
  ) |>
  filter(!is.na(age_numeric)) |>
  select(age_label, age_numeric, n, rlel_months, expected_rmst_months, relative_rlel_percent)
write_csv(age_relative, file.path(analysis_dir, "figure7_relative_rlel_age_source.csv"))

theme_set(
  theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(size = 7.2, face = "bold"),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      panel.grid = element_blank()
    )
)

forest_plot <- rlel_forest |>
  mutate(
    label = paste0(group, "  (n=", comma(n), ")"),
    label = factor(label, levels = rev(label))
  ) |>
  ggplot(aes(
    x = rlel_months, y = label, xmin = rlel_ci_lower, xmax = rlel_ci_upper,
    colour = domain
  )) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey70") +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(
    "Behavior" = "#B64342", "Sex" = "#547AA5", "Tumor group" = "#2A7F62"
  )) +
  labs(
    x = "10-year RLEL (months)", y = NULL, colour = NULL,
    title = "10-year RLEL across prespecified groups"
  ) +
  theme(legend.position = "top")

relative_plot <- ggplot(age_relative, aes(x = age_numeric, y = relative_rlel_percent)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_line(colour = "#0F4D92", linewidth = 0.6) +
  geom_point(colour = "#0F4D92", size = 1.1) +
  labs(
    x = "Age group at diagnosis (years)", y = "Relative RLEL (%)",
    title = "Relative RLEL by age"
  )

std_plot_data <- standardized_os |>
  filter(domain %in% c("Race/ethnicity", "County income")) |>
  mutate(
    group = recode(
      group,
      "Non-Hispanic White" = "NH White",
      "Non-Hispanic Black" = "NH Black",
      "Hispanic (All Races)" = "Hispanic",
      "Non-Hispanic Asian or Pacific Islander" = "NH Asian/PI",
      "Non-Hispanic American Indian/Alaska Native" = "NH AI/AN"
    ),
    group = factor(group, levels = rev(unique(group)))
  )
std_plot <- ggplot(
  std_plot_data,
  aes(
    x = standardized_os_10y, y = group, xmin = ci_lower, xmax = ci_upper,
    colour = domain
  )
) +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.45) +
  geom_point(size = 1.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Race/ethnicity" = "#7B61A8", "County income" = "#C17C2F")) +
  labs(
    x = "Age-standardized 10-year OS", y = NULL, colour = NULL,
    title = "Age-standardized 10-year OS"
  ) +
  theme(legend.position = "top")

period_plot_data <- period_models |>
  filter(term %in% c("raceNon-Hispanic Black", "income_group$110,000+")) |>
  mutate(
    contrast = recode(
      term,
      "raceNon-Hispanic Black" = "NH Black vs NH White",
      "income_group$110,000+" = "$110k+ vs <$60k"
    )
  )
period_plot <- ggplot(
  period_plot_data,
  aes(
    x = hazard_ratio, y = contrast, xmin = ci_lower, xmax = ci_upper,
    colour = period, shape = period
  )
) +
  geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey65", linetype = 2) +
  geom_errorbar(
    orientation = "y", width = 0, linewidth = 0.45,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(size = 1.7, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("2004-2008" = "#547AA5", "2009-2013" = "#B64342")) +
  labs(
    x = "Adjusted hazard ratio", y = NULL,
    colour = "Diagnosis period", shape = "Diagnosis period",
    title = "Adjusted associations by period"
  ) +
  theme(legend.position = "top")

figure7 <- forest_plot + (relative_plot / std_plot / period_plot) +
  plot_layout(widths = c(1.28, 1), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 8),
    plot.tag.position = c(0, 1)
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 145, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(file.path(figure_dir, paste0(stem, ".svg")), width = w, height = h)
  print(plot)
  dev.off()
  cairo_pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  agg_tiff(
    file.path(figure_dir, paste0(stem, ".tiff")), width = w, height = h,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  dev.off()
  agg_png(
    file.path(figure_dir, paste0(stem, ".png")), width = w, height = h,
    units = "in", res = 300
  )
  print(plot)
  dev.off()
}
save_pub(figure7, "Figure7_relative_burden_robustness")

top_other <- other_histology |>
  filter(stable_for_reporting) |>
  slice_head(n = 15)
supp_lines <- c(
  "# Supplementary relative-burden and robustness analyses", "",
  "## Relative RLEL definition", "",
  "Relative RLEL was calculated as RLEL divided by the expected restricted mean survival through the same horizon and expressed as a percentage. It is a proportional burden measure, not a conventional relative-survival ratio.", "",
  "## Direct age-standardization", "",
  "Ten-year overall survival was estimated within five age strata and standardized to the age distribution of the complete 2004-2013 person-level cohort. Confidence intervals use Greenwood variances with fixed standardization weights.", "",
  "## Decomposition of other nonmalignant brain/ONS histologies", "",
  "| ICD-O-3 | Histology | N | Deaths by 10 years | 10-year OS | Stable for reporting |",
  "|---|---|---:|---:|---:|---|",
  apply(top_other, 1, function(x) {
    sprintf(
      "| %s | %s | %s | %s | %.1f%% | %s |",
      x[["histology_code"]], x[["histology"]],
      format(as.numeric(x[["n"]]), big.mark = ","),
      format(as.numeric(x[["events"]]), big.mark = ","),
      100 * as.numeric(x[["survival_10y"]]),
      ifelse(as.logical(x[["stable_for_reporting"]]), "Yes", "No")
    )
  }), "",
  "Reporting stability was prespecified as at least 200 individuals and at least 20 deaths by 10 years. These histology-specific estimates are descriptive and do not replace the prespecified SEER recode groups.", "",
  "## Internal robustness", "",
  "Age-standardized estimates were repeated by race/ethnicity, county-income category, tumor group, and diagnosis period. Adjusted Cox models were also fitted separately in 2004-2008 and 2009-2013. Registry-level validation was not possible because registry identifiers were unavailable in the Research Limited-Field data.", "",
  "## Uncertainty boundary", "",
  "A genuine patient-level bootstrap of RLEL could not be performed because individual expected-survival trajectories were not exported. RLEL intervals therefore use the SEER curve-based Greenwood covariance approximation with expected survival treated as fixed; patient-level bootstrap remains limited to overall survival."
)
write_lines(supp_lines, file.path(table_dir, "supplementary_relative_burden_robustness.md"))
write_lines(
  capture.output(sessionInfo()),
  file.path(figure_dir, "R_sessionInfo_advanced_robustness.txt")
)
message("Advanced robustness analyses completed successfully.")
