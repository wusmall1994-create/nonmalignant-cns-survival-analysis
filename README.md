# Statistical analysis code for nonmalignant brain and CNS tumor survival

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22207634.svg)](https://doi.org/10.5281/zenodo.22207634)

This repository contains the statistical analysis code used to estimate restricted loss in life expectancy (RLEL), compare Ederer II and Pohar–Perme relative-survival estimators, evaluate subgroup heterogeneity, perform individual-level bootstrap and adjusted survival analyses, and generate publication figures.

## Privacy and data availability

No manuscript, raw data, case-level records, SEER\*Stat session files, exported results, patient identifiers, credentials, or personal information are included. The code expects locally generated SEER\*Stat exports. SEER data remain subject to the SEER Research Data Agreement and must be obtained independently by an eligible registered user.

Do not commit files with the extensions `.txt`, `.dic`, or `.ss` to this repository. The repository uses a whitelist `.gitignore` so that only source code and package metadata are tracked.

## Software

- SEER\*Stat 9.0.43.0 for cohort construction and survival exports
- Python 3.12 with the packages listed in `requirements.txt`
- R 4.6.1 with the packages listed in `R-packages.txt`
- PowerShell 7 or Windows PowerShell 5.1 for the orchestration script

## Configuration

Set the following environment variables before running the workflow:

- `SEER_EXPORT_DIR` (required): directory containing the locally generated SEER\*Stat exports.
- `ANALYSIS_OUTPUT_DIR` (optional): destination for generated analyses, figures, tables, reports, and QA logs. Defaults to `outputs/` inside the repository, which is ignored by Git.
- `PYTHON_EXE` (optional): full path to the Python executable. If omitted, `python.exe` on `PATH` is used.
- `RSCRIPT_EXE` (optional): full path to `Rscript.exe`. If omitted, `Rscript.exe` on `PATH` is used.
- `BOOT_REPS` (optional): number of patient-level bootstrap replicates; default 1000 and minimum 200.

Example in PowerShell:

```powershell
$env:SEER_EXPORT_DIR = 'X:\path\to\authorized\seer_exports'
$env:ANALYSIS_OUTPUT_DIR = 'X:\path\to\analysis_outputs'
$env:PYTHON_EXE = 'X:\path\to\python.exe'
$env:RSCRIPT_EXE = 'X:\path\to\Rscript.exe'
$env:BOOT_REPS = '1000'
.\run_all.ps1
```

## Script order

1. `analyze.py` parses SEER\*Stat relative-survival exports, integrates observed and expected survival curves, estimates RLEL, computes subgroup contrasts and age-standardized trends, and writes source-file checksums.
2. `individual_analysis.R` performs patient-level cohort construction, bootstrap overall-survival estimates, adjusted Cox models, proportional-hazards diagnostics, and time-varying sensitivity models.
3. `audit_sessions.R` verifies exported dictionaries, grouping variables, time horizons, estimator labels, and the presence of matching local export/session files.
4. `pohar_perme_analysis.R` compares Ederer II and Pohar–Perme relative-survival estimates and prepares estimator-sensitivity tables and figures.
5. `advanced_robustness_analysis.R` performs direct age standardization, morphology decomposition, diagnosis-period models, and relative RLEL analyses.
6. `make_figures.R` generates the main RLEL and subgroup figures from the analysis outputs.

`run_all.ps1` executes these scripts in the above dependency order.

## Required local exports

The code references the following SEER\*Stat export families. These files are inputs only and are intentionally absent from this repository:

- Ederer II relative-survival exports for 5-, 10-, and 15-year overall analyses
- Ederer II subgroup exports for behavior, brain/CNS recode, C75 primary site, sex, age, and diagnosis year
- Matching Pohar–Perme exports for overall, behavior, sex, and diagnosis-year analyses
- `nonmalignant_CNS_case_listing_2004_2023.txt`
- `nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt`
- Corresponding `.dic` files and, for the local audit only, `.ss` session files

Exact input stems and required columns are defined in the scripts, making missing or incompatible exports fail with an explicit error.

## Outputs

Generated files are organized under `ANALYSIS_OUTPUT_DIR`:

- `analysis/`: analysis tables and checksums
- `figures/`: SVG, PDF, TIFF, and PNG figures plus R session information
- `tables/`: supplementary Markdown tables
- `reports/`: analysis summaries
- `qa_logs/`: session/export audit report

Generated outputs are excluded from version control because they may derive from restricted SEER data.

## Reproducibility and limitations

The repository reproduces the statistical workflow when supplied with authorized exports created using the documented variable names and SEER\*Stat settings. SEER data and SEER\*Stat itself cannot be redistributed here. Software versions are pinned to the environment used for the archived release; later package versions may require minor compatibility changes.

## License

The source code is released under the MIT License. The license applies only to the code and does not grant rights to SEER data or other third-party materials.

## Citation

For an exact software version, cite release `v1.0.0` using DOI [`10.5281/zenodo.22207635`](https://doi.org/10.5281/zenodo.22207635). To cite the software independently of version, use the concept DOI [`10.5281/zenodo.22207634`](https://doi.org/10.5281/zenodo.22207634).
