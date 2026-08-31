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
  library(Matrix)
  library(svglite)
  library(ragg)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
seer_dir <- Sys.getenv("SEER_EXPORT_DIR", unset = "")
if (!nzchar(seer_dir)) stop("Set SEER_EXPORT_DIR to the directory containing the SEER*Stat exports.")
case_file <- file.path(
  seer_dir,
  "nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt"
)
if (!file.exists(case_file)) stop("Case-listing file not found: ", case_file)

output_root <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
analysis_dir <- file.path(output_root, "analysis")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
export_dir <- file.path(output_root, "reports")
qa_dir <- file.path(output_root, "qa_logs")
invisible(lapply(
  c(analysis_dir, figure_dir, table_dir, export_dir, qa_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

set.seed(20260829)
bootstrap_reps <- as.integer(Sys.getenv("BOOT_REPS", unset = "1000"))
if (is.na(bootstrap_reps) || bootstrap_reps < 200L) {
  stop("BOOT_REPS must be an integer of at least 200")
}

cols_needed <- c(
  "Patient ID",
  "Record number recode",
  "Year of diagnosis",
  "Age recode with single ages and 90+",
  "Sex",
  "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)",
  "Behavior code ICD-O-3",
  "SEER Brain and CNS Recode",
  "Diagnostic Confirmation",
  "Type of Reporting Source",
  "Survival months",
  "Survival Days",
  "Vital status recode (study cutoff used)",
  "RX Summ--Surg Prim Site (2003-2022)",
  "RX Summ--Surg Prim Site 2023 (2023+)",
  "Chemotherapy recode (yes, no/unk) (2004+)",
  "Median household income inflation adj to 2024"
)

message("Reading individual-level case listing...")
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
    brain_group = collapse_brain_group(`SEER Brain and CNS Recode`),
    confirmation = collapse_confirmation(`Diagnostic Confirmation`),
    reporting_source = `Type of Reporting Source`,
    survival_months = suppressWarnings(as.numeric(`Survival months`)),
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
    time_days = pmax(survival_days, 0.5, na.rm = TRUE),
    dco_autopsy = grepl("Death certificate|Autopsy", reporting_source, ignore.case = TRUE),
    alive_zero_time = event == 0L & (is.na(survival_days) | survival_days <= 0),
    age_group = case_when(
      age < 18 ~ "0-17",
      age < 40 ~ "18-39",
      age < 60 ~ "40-59",
      age < 75 ~ "60-74",
      TRUE ~ "75+"
    ),
    dx_year_per5 = (year - 2004) / 5
  )

# Patient ID links multiple tumor records. For person-level analyses, retain the
# earliest eligible nonmalignant CNS record per patient, breaking ties by record number.
eligible_record <- dat |>
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
  ungroup()

cohort_flow <- tibble(
  step = c(
    "All exported tumor records",
    "Diagnosis years 2004-2013",
    "After DCO/autopsy and alive-zero-time exclusions",
    "Known age and binary sex",
    "Earliest eligible record per Patient ID",
    "Known race/origin (demographic model)",
    "Known surgery (treatment-associated model)",
    "Known surgery and county income (extended model)"
  ),
  n = c(
    nrow(dat),
    sum(dat$year >= 2004L & dat$year <= 2013L, na.rm = TRUE),
    sum(dat$year >= 2004L & dat$year <= 2013L & !dat$dco_autopsy & !dat$alive_zero_time,
        na.rm = TRUE),
    sum(dat$year >= 2004L & dat$year <= 2013L & !dat$dco_autopsy & !dat$alive_zero_time &
          !is.na(dat$age) & dat$sex %in% c("Male", "Female"), na.rm = TRUE),
    nrow(eligible_record),
    sum(eligible_record$race != "Non-Hispanic Unknown Race", na.rm = TRUE),
    sum(eligible_record$race != "Non-Hispanic Unknown Race" &
          eligible_record$surgery != "Unknown/not recorded", na.rm = TRUE),
    sum(eligible_record$race != "Non-Hispanic Unknown Race" &
          eligible_record$surgery != "Unknown/not recorded" &
          eligible_record$income_group != "Unknown", na.rm = TRUE)
  )
)
write_csv(cohort_flow, file.path(analysis_dir, "individual_analysis_cohort_flow.csv"))

bootstrap_cohort <- eligible_record |>
  mutate(
    followup_months = time_days / 30.4375,
    time10 = pmin(followup_months, 120),
    event10 = as.integer(event == 1L & followup_months <= 120)
  )

message("Running patient-level nonparametric bootstrap (B=", bootstrap_reps, ")...")
max_day <- as.integer(ceiling(120 * 30.4375))
day_index <- pmin(pmax(as.integer(ceiling(bootstrap_cohort$time_days)), 1L), max_day)
n_patients <- nrow(bootstrap_cohort)
at_risk_matrix <- sparseMatrix(
  i = day_index,
  j = seq_len(n_patients),
  x = 1,
  dims = c(max_day, n_patients)
)
death_rows <- which(bootstrap_cohort$event == 1L & bootstrap_cohort$time_days <= max_day)
death_matrix <- sparseMatrix(
  i = day_index[death_rows],
  j = death_rows,
  x = 1,
  dims = c(max_day, n_patients)
)

km_from_multiplicity <- function(multiplicity) {
  totals_by_day <- as.numeric(at_risk_matrix %*% multiplicity)
  deaths_by_day <- as.numeric(death_matrix %*% multiplicity)
  risk <- rev(cumsum(rev(totals_by_day)))
  decrement <- ifelse(risk > 0, pmax(0, 1 - deaths_by_day / risk), 1)
  surv <- cumprod(decrement)
  target_days <- pmin(as.integer(round(c(60, 120) * 30.4375)), max_day)
  surv[target_days]
}

observed_os <- km_from_multiplicity(rep.int(1L, n_patients))
bootstrap_values <- matrix(NA_real_, nrow = bootstrap_reps, ncol = 2L)
colnames(bootstrap_values) <- c("os_5y", "os_10y")
for (b in seq_len(bootstrap_reps)) {
  sampled <- sample.int(n_patients, size = n_patients, replace = TRUE)
  multiplicity <- tabulate(sampled, nbins = n_patients)
  bootstrap_values[b, ] <- km_from_multiplicity(multiplicity)
  if (b %% 100L == 0L) message("  bootstrap replicate ", b, "/", bootstrap_reps)
}

bootstrap_replicates <- as_tibble(bootstrap_values) |>
  mutate(replicate = row_number(), .before = 1L)
write_csv(bootstrap_replicates, file.path(analysis_dir, "bootstrap_replicates_os.csv"))

bootstrap_summary <- tibble(
  horizon_years = c(5L, 10L),
  n_patients = n_patients,
  estimate = observed_os,
  bootstrap_se = apply(bootstrap_values, 2L, sd, na.rm = TRUE),
  ci_lower = apply(bootstrap_values, 2L, quantile, probs = 0.025, na.rm = TRUE),
  ci_upper = apply(bootstrap_values, 2L, quantile, probs = 0.975, na.rm = TRUE),
  bootstrap_reps = bootstrap_reps,
  method = "Patient-level nonparametric percentile bootstrap; earliest eligible record per Patient ID"
)
write_csv(bootstrap_summary, file.path(analysis_dir, "bootstrap_os_5y_10y.csv"))

model_base <- eligible_record |>
  filter(race != "Non-Hispanic Unknown Race") |>
  mutate(
    followup_months = time_days / 30.4375,
    time10 = pmin(followup_months, 120),
    event10 = as.integer(event == 1L & followup_months <= 120),
    race = factor(
      race,
      levels = c(
        "Non-Hispanic White",
        "Non-Hispanic Black",
        "Hispanic (All Races)",
        "Non-Hispanic Asian or Pacific Islander",
        "Non-Hispanic American Indian/Alaska Native"
      )
    ),
    age_group = factor(age_group, levels = c("40-59", "0-17", "18-39", "60-74", "75+")),
    sex = factor(sex, levels = c("Female", "Male")),
    behavior = factor(behavior, levels = c("Benign", "Borderline malignancy")),
    brain_group = factor(
      brain_group,
      levels = c(
        "Meningiomas",
        "Cranial/paraspinal nerve tumors",
        "Other non-malignant brain/ONS",
        "Other included sites/types",
        "Other CNS histologies"
      )
    ),
    confirmation = factor(
      confirmation,
      levels = c("Microscopic/laboratory", "Radiographic/direct visualization", "Clinical/unknown")
    ),
    surgery = factor(surgery, levels = c("No surgery", "Surgery performed", "Unknown/not recorded")),
    chemotherapy = factor(chemotherapy, levels = c("No/Unknown", "Yes")),
    income_group = factor(
      income_group,
      levels = c(
        "<$60,000", "$60,000-$74,999", "$75,000-$89,999",
        "$90,000-$109,999", "$110,000+", "Unknown"
      )
    )
  )

model_demographic_data <- model_base |>
  drop_na(race, age_group, sex, behavior, brain_group, confirmation, dx_year_per5)

model_treatment_data <- model_demographic_data |>
  filter(surgery != "Unknown/not recorded") |>
  droplevels()

model_enhanced_data <- model_demographic_data |>
  filter(
    surgery != "Unknown/not recorded",
    income_group != "Unknown",
    !is.na(chemotherapy)
  ) |>
  droplevels()

message("Fitting adjusted Cox models...")
model_demographic <- coxph(
  Surv(time10, event10) ~ race + age_group + sex + dx_year_per5 + behavior +
    brain_group + strata(confirmation),
  data = model_demographic_data,
  ties = "efron",
  x = TRUE,
  model = TRUE
)
model_treatment <- coxph(
  Surv(time10, event10) ~ race + sex + dx_year_per5 + surgery +
    strata(age_group, behavior, brain_group, confirmation),
  data = model_treatment_data,
  ties = "efron",
  x = TRUE,
  model = TRUE
)
model_enhanced <- coxph(
  Surv(time10, event10) ~ race + sex + dx_year_per5 + chemotherapy + income_group +
    strata(age_group, behavior, brain_group, confirmation, surgery),
  data = model_enhanced_data,
  ties = "efron",
  x = TRUE,
  model = TRUE
)

tidy_cox <- function(model, model_name) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      model = model_name,
      term,
      hazard_ratio = estimate,
      ci_lower = conf.low,
      ci_upper = conf.high,
      p_value = p.value
    )
}

cox_results <- bind_rows(
  tidy_cox(model_demographic, "Demographic/tumor-adjusted"),
  tidy_cox(model_treatment, "Demographic/tumor + surgery"),
  tidy_cox(model_enhanced, "Extended SES/chemotherapy-adjusted")
)
write_csv(cox_results, file.path(analysis_dir, "adjusted_cox_models.csv"))
write_csv(
  tidy_cox(model_enhanced, "Extended SES/chemotherapy-adjusted"),
  file.path(analysis_dir, "enhanced_ses_chemotherapy_cox.csv")
)

model_performance <- bind_rows(
  glance(model_demographic) |>
    mutate(model = "Demographic/tumor-adjusted", n = nrow(model_demographic_data),
           events = sum(model_demographic_data$event10)),
  glance(model_treatment) |>
    mutate(model = "Demographic/tumor + surgery", n = nrow(model_treatment_data),
           events = sum(model_treatment_data$event10)),
  glance(model_enhanced) |>
    mutate(model = "Extended SES/chemotherapy-adjusted", n = nrow(model_enhanced_data),
           events = sum(model_enhanced_data$event10))
) |>
  select(model, n, events, concordance, std.error.concordance, statistic.log, p.value.log,
         statistic.sc, p.value.sc)
write_csv(model_performance, file.path(analysis_dir, "adjusted_cox_model_performance.csv"))

ph_to_tibble <- function(model, model_name) {
  z <- cox.zph(model, transform = "km")
  as.data.frame(z$table) |>
    tibble::rownames_to_column("term") |>
    as_tibble() |>
    rename(chisq = chisq, df = df, p_value = p) |>
    mutate(model = model_name, .before = 1L)
}

ph_results <- bind_rows(
  ph_to_tibble(model_demographic, "Demographic/tumor-adjusted"),
  ph_to_tibble(model_treatment, "Demographic/tumor + surgery"),
  ph_to_tibble(model_enhanced, "Extended SES/chemotherapy-adjusted")
)
write_csv(ph_results, file.path(analysis_dir, "cox_proportional_hazards_tests.csv"))

# The global proportional-hazards test is expected to be sensitive in a cohort
# this large. Estimate interval-specific associations while stratifying the
# baseline hazard by the covariates with the strongest non-proportionality.
fit_interval_model <- function(data, start_month, end_month, interval_label) {
  interval_data <- data |>
    filter(followup_months > start_month) |>
    mutate(
      interval_time = pmin(followup_months, end_month) - start_month,
      interval_event = as.integer(event == 1L & followup_months <= end_month)
    )
  fit <- coxph(
    Surv(interval_time, interval_event) ~ race + sex + dx_year_per5 + surgery +
      strata(age_group, behavior, brain_group, confirmation),
    data = interval_data,
    ties = "efron"
  )
  list(
    results = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      transmute(
        interval = interval_label,
        start_month,
        end_month,
        n = nrow(interval_data),
        events = sum(interval_data$interval_event),
        term,
        hazard_ratio = estimate,
        ci_lower = conf.low,
        ci_upper = conf.high,
        p_value = p.value
      ),
    fit = fit
  )
}

interval_fits <- list(
  fit_interval_model(model_treatment_data, 0, 12, "0-1 year"),
  fit_interval_model(model_treatment_data, 12, 60, "1-5 years"),
  fit_interval_model(model_treatment_data, 60, 120, "5-10 years")
)
time_varying_results <- bind_rows(lapply(interval_fits, `[[`, "results"))
write_csv(time_varying_results, file.path(analysis_dir, "time_varying_cox_results.csv"))

fit_enhanced_interval_model <- function(data, start_month, end_month, interval_label) {
  interval_data <- data |>
    filter(followup_months > start_month) |>
    mutate(
      interval_time = pmin(followup_months, end_month) - start_month,
      interval_event = as.integer(event == 1L & followup_months <= end_month)
    )
  fit <- coxph(
    Surv(interval_time, interval_event) ~ race + sex + dx_year_per5 +
      chemotherapy + income_group +
      strata(age_group, behavior, brain_group, confirmation, surgery),
    data = interval_data,
    ties = "efron"
  )
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    transmute(
      interval = interval_label,
      start_month,
      end_month,
      n = nrow(interval_data),
      events = sum(interval_data$interval_event),
      term,
      hazard_ratio = estimate,
      ci_lower = conf.low,
      ci_upper = conf.high,
      p_value = p.value
    )
}

enhanced_time_varying_results <- bind_rows(
  fit_enhanced_interval_model(model_enhanced_data, 0, 12, "0-1 year"),
  fit_enhanced_interval_model(model_enhanced_data, 12, 60, "1-5 years"),
  fit_enhanced_interval_model(model_enhanced_data, 60, 120, "5-10 years")
)
write_csv(
  enhanced_time_varying_results,
  file.path(analysis_dir, "enhanced_time_varying_cox_results.csv")
)

race_counts <- model_demographic_data |>
  count(race, name = "n") |>
  mutate(percent = 100 * n / sum(n))
surgery_counts <- model_treatment_data |>
  count(surgery, name = "n") |>
  mutate(percent = 100 * n / sum(n))
income_counts <- model_enhanced_data |>
  count(income_group, name = "n") |>
  mutate(percent = 100 * n / sum(n))
chemotherapy_counts <- model_enhanced_data |>
  count(chemotherapy, name = "n") |>
  mutate(percent = 100 * n / sum(n))
write_csv(race_counts, file.path(analysis_dir, "adjusted_model_race_counts.csv"))
write_csv(surgery_counts, file.path(analysis_dir, "adjusted_model_surgery_counts.csv"))
write_csv(income_counts, file.path(analysis_dir, "enhanced_model_income_counts.csv"))
write_csv(chemotherapy_counts, file.path(analysis_dir, "enhanced_model_chemotherapy_counts.csv"))

term_labels <- c(
  "raceNon-Hispanic Black" = "Race: non-Hispanic Black",
  "raceHispanic (All Races)" = "Race: Hispanic",
  "raceNon-Hispanic Asian or Pacific Islander" = "Race: non-Hispanic Asian/Pacific Islander",
  "raceNon-Hispanic American Indian/Alaska Native" = "Race: non-Hispanic American Indian/Alaska Native",
  "sexMale" = "Sex: male"
)

forest_main <- cox_results |>
  filter(model == "Demographic/tumor + surgery", term %in% names(term_labels)) |>
  mutate(
    label = unname(term_labels[term]),
    association_group = "Adjusted average association"
  ) |>
  select(label, hazard_ratio, ci_lower, ci_upper, p_value, association_group)

forest_surgery <- time_varying_results |>
  filter(term == "surgerySurgery performed") |>
  transmute(
    label = paste0("Surgery: ", interval),
    hazard_ratio,
    ci_lower,
    ci_upper,
    p_value,
    association_group = "Time-specific surgery association"
  )

forest_data <- bind_rows(forest_main, forest_surgery) |>
  mutate(
    label = factor(
      label,
      levels = rev(c(
        unname(term_labels),
        "Surgery: 0-1 year", "Surgery: 1-5 years", "Surgery: 5-10 years"
      ))
    )
  )
write_csv(forest_data, file.path(analysis_dir, "figure5_adjusted_associations_source.csv"))

blue <- "#0F4D92"
rose <- "#B64342"
neutral <- "#707070"
light <- "#D9E2F2"

theme_set(
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.8),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.8, colour = neutral),
      plot.tag = element_text(size = 9, face = "bold"),
      panel.grid = element_blank()
    )
)

p5a <- ggplot(forest_data, aes(x = hazard_ratio, y = label, colour = association_group)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = 2, colour = "#A0A0A0") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0, linewidth = 0.55) +
  geom_point(size = 2.2) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", hazard_ratio, ci_lower, ci_upper)),
            hjust = -0.03, nudge_x = 0.03, size = 2.15, colour = neutral,
            show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Adjusted average association" = blue,
    "Time-specific surgery association" = rose
  )) +
  scale_x_log10(expand = expansion(mult = c(0.04, 0.34))) +
  labs(
    title = "Adjusted mortality associations",
    x = "Adjusted hazard ratio for all-cause death (log scale)",
    y = NULL,
    colour = NULL,
    subtitle = "10-year follow-up; surgery HRs are interval-specific"
  ) +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "top")

bootstrap_plot <- bootstrap_summary |>
  mutate(horizon = factor(paste0(horizon_years, " years"), levels = c("5 years", "10 years")))

p5b <- ggplot(bootstrap_plot, aes(x = horizon, y = estimate)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.14,
                linewidth = 0.6, colour = rose) +
  geom_point(size = 2.5, colour = rose) +
  geom_text(aes(label = percent(estimate, accuracy = 0.1)), vjust = -1,
            size = 2.4, colour = neutral) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.08, 0.16))) +
  labs(
    title = "Bootstrap overall survival",
    x = "Follow-up horizon",
    y = "Observed overall survival",
    subtitle = paste0("Patient level, B=", comma(bootstrap_reps))
  )

figure5 <- p5a + p5b + plot_annotation(tag_levels = "a") +
  plot_layout(widths = c(1.8, 1))

save_pub <- function(plot, stem, width_mm = 183, height_mm = 105, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite(file.path(figure_dir, paste0(stem, ".svg")), width = width_in, height = height_in)
  print(plot)
  dev.off()
  cairo_pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = width_in, height = height_in,
            family = "sans")
  print(plot)
  dev.off()
  agg_tiff(file.path(figure_dir, paste0(stem, ".tiff")), width = width_in, height = height_in,
           units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  agg_png(file.path(figure_dir, paste0(stem, ".png")), width = width_in, height = height_in,
          units = "in", res = 300)
  print(plot)
  dev.off()
}

save_pub(figure5, "Figure5_adjusted_associations_bootstrap")

key_terms <- cox_results |>
  filter(model == "Demographic/tumor + surgery", term %in% names(term_labels)) |>
  mutate(label = unname(term_labels[term])) |>
  select(label, hazard_ratio, ci_lower, ci_upper, p_value)

fmt_hr <- function(label_text) {
  row <- key_terms |> filter(label == label_text)
  if (nrow(row) != 1L) return("not estimable")
  sprintf("HR %.2f (95%% CI %.2f-%.2f; p=%s)", row$hazard_ratio, row$ci_lower,
          row$ci_upper, format.pval(row$p_value, digits = 3, eps = 0.001))
}

global_ph <- ph_results |>
  filter(model == "Demographic/tumor + surgery", term == "GLOBAL") |>
  pull(p_value)

surgery_interval_terms <- time_varying_results |>
  filter(term == "surgerySurgery performed")

enhanced_term_labels <- c(
  "raceNon-Hispanic Black" = "Race: non-Hispanic Black",
  "chemotherapyYes" = "Recorded chemotherapy: yes",
  "income_group$60,000-$74,999" = "County income: $60,000-$74,999",
  "income_group$75,000-$89,999" = "County income: $75,000-$89,999",
  "income_group$90,000-$109,999" = "County income: $90,000-$109,999",
  "income_group$110,000+" = "County income: $110,000+"
)

enhanced_key_terms <- cox_results |>
  filter(model == "Extended SES/chemotherapy-adjusted", term %in% names(enhanced_term_labels)) |>
  mutate(label = unname(enhanced_term_labels[term])) |>
  select(label, hazard_ratio, ci_lower, ci_upper, p_value)

enhanced_interval_terms <- enhanced_time_varying_results |>
  filter(term %in% c("chemotherapyYes", "income_group$110,000+")) |>
  mutate(
    label = recode(
      term,
      chemotherapyYes = "Recorded chemotherapy: yes",
      `income_group$110,000+` = "County income: $110,000+"
    )
  )

fmt_surgery_interval <- function(interval_text) {
  row <- surgery_interval_terms |> filter(interval == interval_text)
  if (nrow(row) != 1L) return("not estimable")
  sprintf("%s: HR %.2f (95%% CI %.2f-%.2f)", interval_text, row$hazard_ratio,
          row$ci_lower, row$ci_upper)
}

report_lines <- c(
  "# Additional individual-level analyses",
  "",
  "## Completed analyses",
  "",
  paste0("- Person-level analytic cohort: **", format(n_patients, big.mark = ","),
         "** patients diagnosed in 2004-2013 after standard exclusions and retention of the earliest eligible record per Patient ID."),
  paste0("- Individual-level nonparametric bootstrap: **", bootstrap_reps,
         "** replicates; 5-year OS ", percent(bootstrap_summary$estimate[1], accuracy = 0.1),
         " (95% percentile CI ", percent(bootstrap_summary$ci_lower[1], accuracy = 0.1), "-",
         percent(bootstrap_summary$ci_upper[1], accuracy = 0.1), "); 10-year OS ",
         percent(bootstrap_summary$estimate[2], accuracy = 0.1), " (",
         percent(bootstrap_summary$ci_lower[2], accuracy = 0.1), "-",
         percent(bootstrap_summary$ci_upper[2], accuracy = 0.1), ")."),
  paste0("- Adjusted race association (Black vs White): ", fmt_hr("Race: non-Hispanic Black"), "."),
  paste0("- Adjusted race association (Hispanic vs White): ", fmt_hr("Race: Hispanic"), "."),
  paste0("- Time-specific surgery associations: ",
         paste(vapply(c("0-1 year", "1-5 years", "5-10 years"),
                      fmt_surgery_interval, character(1)), collapse = "; "),
         ". These are associations, not causal treatment effects."),
  paste0("- Global proportional-hazards test for the surgery-adjusted model: p=",
         format.pval(global_ph, digits = 3, eps = 0.001),
         "; because of the very large sample, term-specific diagnostics and effect sizes should be emphasized over significance alone."),
  paste0("- Extended complete-case cohort with county income and chemotherapy: **",
         format(nrow(model_enhanced_data), big.mark = ","),
         "** patients. Relative to county income <$60,000, the adjusted HR was ",
         sprintf("%.2f (95%% CI %.2f-%.2f)",
                 enhanced_key_terms$hazard_ratio[enhanced_key_terms$label == "County income: $110,000+"],
                 enhanced_key_terms$ci_lower[enhanced_key_terms$label == "County income: $110,000+"],
                 enhanced_key_terms$ci_upper[enhanced_key_terms$label == "County income: $110,000+"]),
         " for county income $110,000+. Recorded chemotherapy yes versus no/unknown had adjusted HR ",
         sprintf("%.2f (95%% CI %.2f-%.2f)",
                 enhanced_key_terms$hazard_ratio[enhanced_key_terms$label == "Recorded chemotherapy: yes"],
                 enhanced_key_terms$ci_lower[enhanced_key_terms$label == "Recorded chemotherapy: yes"],
                 enhanced_key_terms$ci_upper[enhanced_key_terms$label == "Recorded chemotherapy: yes"]),
         "."),
  "",
  "## Interpretation boundary",
  "",
  "- These models estimate all-cause mortality associations, not relative survival or RLEL regression.",
  "- Surgery is vulnerable to confounding by indication, tumor severity, resectability, and incomplete clinical detail.",
  "- Patient IDs with multiple tumor records were reduced to the earliest eligible record for person-level inference.",
  "- County median household income and recorded chemotherapy were included in an extended model; county income is an ecological measure and the chemotherapy comparator combines no treatment with unknown treatment.",
  "- Rurality and radiotherapy were not present in the revised case-listing export and therefore were not modeled.",
  "- A genuine patient-level RLEL bootstrap remains unavailable because individual expected-survival trajectories were not exported.",
  "",
  "## SEER*Stat sensitivity export",
  "",
  "Pohar-Perme exports are available and are processed by `pohar_perme_analysis.R`.",
  "",
  "Available comparisons include overall 5-, 10-, and 15-year estimates; behavior and sex at 10 years; and diagnosis-year and age-by-diagnosis-year estimates at 5 years.",
  "",
  "## Data still unavailable",
  "",
  "- Rural-Urban Continuum Codes.",
  "- Radiation recode."
)
write_lines(report_lines, file.path(export_dir, "additional_individual_analysis_report.md"))

table_lines <- c(
  "# Supplementary adjusted and bootstrap results",
  "",
  "## Patient-level bootstrap overall survival",
  "",
  paste0("| Horizon | Patients | Estimate | Bootstrap SE | 95% percentile CI | Replicates |"),
  "|---|---:|---:|---:|---:|---:|",
  apply(bootstrap_summary, 1L, function(x) {
    sprintf("| %s years | %s | %.3f | %.4f | %.3f-%.3f | %s |",
            x[["horizon_years"]], format(as.numeric(x[["n_patients"]]), big.mark = ","),
            as.numeric(x[["estimate"]]), as.numeric(x[["bootstrap_se"]]),
            as.numeric(x[["ci_lower"]]), as.numeric(x[["ci_upper"]]),
            x[["bootstrap_reps"]])
  }),
  "",
  "## Selected surgery-adjusted Cox estimates",
  "",
  "| Contrast | Adjusted HR | 95% CI | p value |",
  "|---|---:|---:|---:|",
  apply(key_terms, 1L, function(x) {
    sprintf("| %s | %.2f | %.2f-%.2f | %s |", x[["label"]],
            as.numeric(x[["hazard_ratio"]]), as.numeric(x[["ci_lower"]]),
            as.numeric(x[["ci_upper"]]),
            format.pval(as.numeric(x[["p_value"]]), digits = 3, eps = 0.001))
  }),
  "",
  "## Time-specific surgery associations",
  "",
  "| Interval | Patients at interval entry | Deaths | Adjusted HR | 95% CI | p value |",
  "|---|---:|---:|---:|---:|---:|",
  apply(surgery_interval_terms, 1L, function(x) {
    sprintf("| %s | %s | %s | %.2f | %.2f-%.2f | %s |", x[["interval"]],
            format(as.numeric(x[["n"]]), big.mark = ","),
            format(as.numeric(x[["events"]]), big.mark = ","),
            as.numeric(x[["hazard_ratio"]]), as.numeric(x[["ci_lower"]]),
            as.numeric(x[["ci_upper"]]),
            format.pval(as.numeric(x[["p_value"]]), digits = 3, eps = 0.001))
  }),
  "",
  "## Extended county-income and chemotherapy-adjusted Cox estimates",
  "",
  "The model additionally included county median household income and recorded chemotherapy, while retaining race/ethnicity, sex, diagnosis year, and stratification by age group, behavior, tumor group, diagnostic confirmation, and surgery. Income <$60,000 and chemotherapy no/unknown were reference categories.",
  "",
  "| Contrast | Adjusted HR | 95% CI | p value |",
  "|---|---:|---:|---:|",
  apply(enhanced_key_terms, 1L, function(x) {
    sprintf("| %s | %.2f | %.2f-%.2f | %s |", x[["label"]],
            as.numeric(x[["hazard_ratio"]]), as.numeric(x[["ci_lower"]]),
            as.numeric(x[["ci_upper"]]),
            format.pval(as.numeric(x[["p_value"]]), digits = 3, eps = 0.001))
  }),
  "",
  "## Time-specific extended-model sensitivity estimates",
  "",
  "| Contrast | Interval | Patients at interval entry | Deaths | Adjusted HR | 95% CI | p value |",
  "|---|---|---:|---:|---:|---:|---:|",
  apply(enhanced_interval_terms, 1L, function(x) {
    sprintf("| %s | %s | %s | %s | %.2f | %.2f-%.2f | %s |", x[["label"]],
            x[["interval"]], format(as.numeric(x[["n"]]), big.mark = ","),
            format(as.numeric(x[["events"]]), big.mark = ","),
            as.numeric(x[["hazard_ratio"]]), as.numeric(x[["ci_lower"]]),
            as.numeric(x[["ci_upper"]]),
            format.pval(as.numeric(x[["p_value"]]), digits = 3, eps = 0.001))
  }),
  "",
  "County income is ecological. Recorded chemotherapy compares yes with no/unknown and is sparse; all treatment estimates are prognostic associations, not causal effects."
)
write_lines(table_lines, file.path(table_dir, "supplementary_adjusted_bootstrap_tables.md"))

write_lines(capture.output(sessionInfo()), file.path(figure_dir, "R_sessionInfo_individual_analysis.txt"))
message("Individual-level analyses completed successfully.")
