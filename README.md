# Nonmalignant CNS tumor survival analysis

Final statistical analysis code for ten-year expected minus observed restricted mean survival time (RMST) and relative survival in nonmalignant central nervous system tumor groups.

## Contents

- `analyze.py`: export parsing; actuarial RMST and Greenwood covariance; common internal age standardization; all-age and actual age profiles; age-range sensitivity; histology/behavior composition; landmark relative survival; conditional survival; exploratory diagnosis-period contrasts; fixed-reference resampling; and diagnostic-confirmation checks.
- `make_figures.py`: main and supplementary figures from locally generated aggregate tables.
- `run_analysis.py`: runs both scripts in order.
- `INPUTS.md`: authorized input filenames, required fields and export conventions.
- `requirements.txt`: tested Python dependencies (Python 3.12.14).

## Run

Use a Python environment with the listed dependencies. Obtain the required SEER exports independently. In PowerShell:

```powershell
$env:SEER_EXPORT_DIR = 'X:\authorized_exports'
$env:SEER_SURVIVAL_EXPORT_DIR = 'X:\authorized_exports\survival_exports'
$env:ANALYSIS_OUTPUT_DIR = 'X:\local_results'
$env:PRIVATE_WORK_DIR = 'X:\private_intermediates'
python run_analysis.py
```

`SEER_EXPORT_DIR` is required. The other defaults are its `survival_exports` subdirectory, `outputs`, and `work/private`, respectively. Relative paths resolve against the working directory. Private intermediate pickles contain case-level records and must remain local. Figure outputs are TIFF/SVG, with local PNG previews.

## Statistical definitions

The primary horizon is 120 months. Monthly trapezoidal integration estimates observed and expected RMST separately. The difference is expected minus observed RMST. Identical empirical five-year age weights from the primary cohort at ages 20–79 standardize both components. Proportional difference is 100 times the standardized difference divided by standardized expected RMST. Landmark relative survival is the ratio of separately standardized survival probabilities, rather than an average of age-specific ratios.

Confidence intervals use Greenwood covariance for observed survival, conditional on fixed expected curves and age weights. Fixed-reference resampling uses 2,000 multinomial draws with seed 20261024. It does not repeat individual life-table matching. Incomplete age-specific curves are not extrapolated. Confirmation-specific comparisons use each subgroup's own Ederer I reference. Conditional survivor results are unstandardized all-cause outcomes; period comparisons are exploratory.

## Data access and privacy

Only source code and its operating documentation are distributed. No manuscript, original or aggregate dataset, individual record, identifier value, session file, dictionary, result, credential or local user path is included. Field names such as `Patient ID` describe the required input schema and contain no actual identifiers. SEER data and life-table access remain subject to SEER requirements: https://seer.cancer.gov/data/access.html . Execution requires authorized local inputs; this repository cannot reproduce numerical results without them.

Generated records and results are excluded by the allowlist `.gitignore`. Assertions target the documented study cohort, rather than arbitrary replacement cohorts. The public license applies to code only and grants no rights to third-party data.

## License and citation

MIT License. Cite the exact archived release using `CITATION.cff` and the DOI displayed in the release record.
