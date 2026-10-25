# Authorized local inputs

Input files are not distributed. Each survival export needs its same-stem `.dic` so numeric category codes can be decoded. Filenames are case-sensitive on non-Windows systems.

## SEER_EXPORT_DIR

- `nonmalignant_CNS_RS_10y_overall.txt`
- `nonmalignant_CNS_RS_5y_overall.txt`
- `nonmalignant_CNS_RS_15y_overall.txt`
- `nonmalignant_CNS_RS_10y_brain_recode.txt`
- `nonmalignant_CNS_RS_10y_C75_by_primary_site.txt`
- `nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt`

## SEER_SURVIVAL_EXPORT_DIR

- `Nonmalignant_CNS_10y_cohort_audit.txt`
- `nonmalignant_CNS_RS_10y_brain_recode_age5_EdererI.txt`
- `nonmalignant_CNS_RS_10y_brain_recode_age5_EdererII.txt`
- `Nonmalignant_CNS_RS_10y_brain_recode_allage_EdererI.txt`
- `nonmalignant_CNS_RS_10y_brain_recode_confirmation_age5_EdererI.txt`
- `nonmalignant_CNS_RS_10y_C75_site_allage_EdererI.txt`
- `Nonmalignant_CNS_RS_EdererII_10y_C75_site_age5.txt`
- `Nonmalignant_CNS_RS_EdererII_5y_dxyear_brain_recode_age5.txt`
- `Nonmalignant_CNS_RS_EdererII_5y_dxyear_confirmation_age5.txt`
- `Nonmalignant_CNS_RS_EdererI_10y_C75_site_age5.txt`
- `Nonmalignant_CNS_RS_EdererI_10y_dx2004_2013_overall.txt`
- `Nonmalignant_CNS_RS_EdererI_15y_dx2004_2008_overall.txt`
- `Nonmalignant_CNS_RS_EdererI_5y_dx2004_2018_overall_Life.txt`
- `Nonmalignant_CNS_RS_EdererI_5y_dx2004_2018_overall_Summary.txt`
- `Nonmalignant_CNS_RS_EdererI_5y_dxyear_brain_recode_age5.txt`
- `Nonmalignant_CNS_RS_EdererI_5y_dxyear_confirmation_age5.txt`

## Export and matching conventions

SEER*Stat 9.0.43.0; SEER Research Limited-Field Data, 21 Registries (excluding IL), November 2025 submission, 2000–2023, linked county attributes. Primary diagnosis years 2004–2013; administrative cutoff December 2023; actuarial monthly survival. Overall five- and fifteen-year contexts use diagnosis years 2004–2018 and 2004–2008. Primary selected behavior is /0 or /1; registry sequence selection retains only one state registry-defined neoplasm or the first of two or more. Death-certificate/autopsy, unknown or life-table-unsupported age, and alive zero-survival observations are excluded by the survival session.

Expected survival references: U.S. by SES/geography/race (NHW,NHB,NHAIAN,NHAPI,HISP), 1992–2023, ages 0–99, state-county (modeled by varied state-county-SES). Survival exports are tab-delimited UTF-8 with headers, numeric category codes and dictionary formats, proportions, no thousands separators, and NA for missing values. Cohort listings use labels. Export the full result matrix, including time zero and all required intervals.

## Required fields

Survival matrix: `Page type`, grouping columns, `Interval`, `Observed Survival (Cum)`, `Expected Survival (Cum)`, `Alive at Start`, `Died`, `Lost to Followup`, `Adjusted Alive`. Grouping columns include `Age recode with <1 year olds and 90+`, `SEER Brain and CNS Recode`, `Primary Site - labeled`, `Year of diagnosis`, and `Diagnostic Confirmation`, as indicated by each filename. Primary curves require months 0–120; five- and fifteen-year overall exports require 0–60 and 0–180.

Exact cohort listing: `Patient ID`, `Record number recode`, `Age recode with single ages and 90+`, `Number of Intervals (Calculated)`, `End Calc Vital Status (Adjusted)`, `Primary Site`, `Behavior code ICD-O-3`, `SEER Brain and CNS Recode`, `Diagnostic Confirmation`. The local case listing supplies `Histologic Type ICD-O-3` and `ICD-O-3 Hist/behav` for record-level classification. Keys are normalized locally and the merge is checked one-to-one; no key values are published.

The all-age Ederer I files must have one grouping variable, with filenames matching `*brain_recode*allage*EdererI.txt` or `*C75*allage*EdererI.txt`. Keep exactly one intended matching file per pattern. Confirmation-stratified survival exports retain tumor group, five-year age and diagnostic confirmation. Parent microscopy categories overlap their child categories and are excluded from leaf-category sums.
