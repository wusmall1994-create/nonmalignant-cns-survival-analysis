from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


if "SEER_EXPORT_DIR" not in os.environ:
    raise RuntimeError("Set SEER_EXPORT_DIR to the directory containing the SEER*Stat exports.")
SEER_DIR = Path(os.environ["SEER_EXPORT_DIR"])
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = Path(os.environ.get("ANALYSIS_OUTPUT_DIR", SCRIPT_DIR / "outputs"))
OUT = OUTPUT_ROOT / "analysis"
OUT.mkdir(parents=True, exist_ok=True)


EXPORTS = {
    "overall_5y": "nonmalignant_CNS_RS_5y_overall",
    "overall_10y": "nonmalignant_CNS_RS_10y_overall",
    "overall_15y": "nonmalignant_CNS_RS_15y_overall",
    "behavior_10y": "nonmalignant_CNS_RS_10y_behavior",
    "brain_recode_10y": "nonmalignant_CNS_RS_10y_brain_recode",
    "c75_site_10y": "nonmalignant_CNS_RS_10y_C75_by_primary_site",
    "sex_10y": "nonmalignant_CNS_RS_10y_sex",
    "age_10y": "nonmalignant_CNS_RS_10y_age5",
    "dxyear_5y": "nonmalignant_CNS_RS_5y_dxyear_v2",
    "dxyear_age_5y": "nonmalignant_CNS_RS_5y_dxyear_age5",
}

POHAR_PERME_EXPORTS = [
    "nonmalignant_CNS_NS_PP_5y_overall",
    "nonmalignant_CNS_NS_PP_10y_overall",
    "nonmalignant_CNS_NS_PP_15y_overall",
    "nonmalignant_CNS_NS_PP_10y_behavior",
    "nonmalignant_CNS_NS_PP_10y_sex",
    "nonmalignant_CNS_NS_PP_5y_dxyear",
    "nonmalignant_CNS_NS_PP_5y_dxyear_age5",
]


@dataclass
class ExportData:
    summary: pd.DataFrame
    life: pd.DataFrame
    formats: dict[str, dict[str, str]]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_formats(dic_path: Path) -> dict[str, dict[str, str]]:
    formats: dict[str, dict[str, str]] = {}
    current: str | None = None
    with dic_path.open("r", encoding="utf-8-sig", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            match = re.match(r"\[Format=(.+)\]", line)
            if match:
                current = match.group(1)
                formats[current] = {}
                continue
            if line.startswith("["):
                current = None
                continue
            if current and "=" in line:
                code, label = line.split("=", 1)
                formats[current][code.strip()] = label.strip().strip('"')
    return formats


def parse_export(base_name: str) -> ExportData:
    txt_path = SEER_DIR / f"{base_name}.txt"
    dic_path = SEER_DIR / f"{base_name}.dic"
    rows: dict[str, list[dict[str, str]]] = {"summary": [], "life": []}
    header: list[str] | None = None
    kind: str | None = None
    with txt_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        for row in reader:
            if not row:
                continue
            if row[0] == "Page type":
                header = [x.replace("~", " ") for x in row]
                if "Summary Interval" in header:
                    kind = "summary"
                elif "Interval" in header:
                    kind = "life"
                else:
                    kind = None
                continue
            if header is None or kind is None or len(row) != len(header):
                continue
            if row[0] not in {"0", "1"}:
                continue
            rows[kind].append(dict(zip(header, row)))

    summary = pd.DataFrame(rows["summary"])
    life = pd.DataFrame(rows["life"])
    formats = parse_formats(dic_path)
    for frame in (summary, life):
        if frame.empty:
            continue
        for col in frame.columns:
            if col in formats:
                frame[f"{col} label"] = frame[col].map(formats[col]).fillna(frame[col])
    return ExportData(summary=summary, life=life, formats=formats)


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = frame.copy()
    for col in columns:
        if col in result:
            result[col] = pd.to_numeric(result[col], errors="coerce")
    return result


def trap_weights(n: int) -> np.ndarray:
    w = np.ones(n, dtype=float)
    if n > 1:
        w[0] = 0.5
        w[-1] = 0.5
    return w


def rlel_from_curve(curve: pd.DataFrame, tau: int) -> dict[str, float]:
    cols = [
        "Interval",
        "Observed Survival (Cum)",
        "Expected Survival (Cum)",
        "Relative Survival (Cum)",
        "Observed SE (Cum)",
        "Relative SE (Cum)",
        "Alive at Start",
    ]
    curve = numeric(curve, cols).sort_values("Interval")
    curve = curve[(curve["Interval"] >= 0) & (curve["Interval"] <= tau)].copy()
    if curve.empty or int(curve["Interval"].max()) != tau:
        raise ValueError(f"Curve does not contain complete 0-{tau} month follow-up")
    if not np.array_equal(curve["Interval"].to_numpy(), np.arange(tau + 1)):
        raise ValueError(f"Curve has missing monthly intervals within 0-{tau}")

    obs = curve["Observed Survival (Cum)"].to_numpy(float)
    exp = curve["Expected Survival (Cum)"].to_numpy(float)
    se_obs = curve["Observed SE (Cum)"].fillna(0).to_numpy(float)
    w = trap_weights(len(curve))
    observed_rmst = float(np.dot(w, obs))
    expected_rmst = float(np.dot(w, exp))
    rlel = float(np.dot(w, exp - obs))

    # Greenwood covariance approximation for nested cumulative survival estimates.
    # Expected survival comes from population life tables and is treated as fixed.
    greenwood = np.divide(se_obs, obs, out=np.zeros_like(se_obs), where=obs > 0) ** 2
    min_idx = np.minimum.outer(np.arange(len(obs)), np.arange(len(obs)))
    cov = np.outer(obs, obs) * greenwood[min_idx]
    variance = float(w @ cov @ w)
    se = math.sqrt(max(variance, 0.0))
    z = 1.959963984540054
    endpoint = curve.iloc[-1]
    return {
        "tau_months": tau,
        "n": int(curve.iloc[0]["Alive at Start"]),
        "observed": float(endpoint["Observed Survival (Cum)"]),
        "expected": float(endpoint["Expected Survival (Cum)"]),
        "relative": float(endpoint["Relative Survival (Cum)"]),
        "observed_rmst_months": observed_rmst,
        "expected_rmst_months": expected_rmst,
        "rlel_months": rlel,
        "relative_rlel_percent": 100 * rlel / expected_rmst if expected_rmst > 0 else np.nan,
        "rlel_se": se,
        "rlel_ci_lower": rlel - z * se,
        "rlel_ci_upper": rlel + z * se,
    }


def group_columns(life: pd.DataFrame) -> list[str]:
    interval_pos = life.columns.get_loc("Interval")
    return [c for c in life.columns[1:interval_pos] if not c.endswith(" label")]


def rlel_table(data: ExportData, taus: list[int]) -> pd.DataFrame:
    life = data.life[data.life["Page type"] == "0"].copy()
    groups = group_columns(life)
    records: list[dict[str, object]] = []
    iterator = [((), life)] if not groups else life.groupby(groups, dropna=False, sort=False)
    for key, curve in iterator:
        keys = key if isinstance(key, tuple) else (key,)
        base: dict[str, object] = {}
        for col, value in zip(groups, keys):
            base[col] = value
            base[f"{col}_label"] = data.formats.get(col, {}).get(str(value), str(value))
        n0 = pd.to_numeric(curve["Alive at Start"], errors="coerce").max()
        if not np.isfinite(n0) or n0 <= 0:
            continue
        max_interval = int(pd.to_numeric(curve["Interval"], errors="coerce").max())
        for tau in taus:
            if max_interval < tau:
                continue
            try:
                records.append(base | rlel_from_curve(curve, tau))
            except ValueError:
                continue
    return pd.DataFrame(records)


def add_pairwise_differences(table: pd.DataFrame, group_label: str, tau: int) -> pd.DataFrame:
    sub = table[table["tau_months"] == tau].reset_index(drop=True)
    records = []
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            a, b = sub.iloc[i], sub.iloc[j]
            diff = float(a["rlel_months"] - b["rlel_months"])
            se = math.sqrt(float(a["rlel_se"] ** 2 + b["rlel_se"] ** 2))
            zstat = diff / se if se > 0 else np.nan
            p = float(math.erfc(abs(zstat) / math.sqrt(2))) if np.isfinite(zstat) else np.nan
            records.append(
                {
                    "group_a": a[group_label],
                    "group_b": b[group_label],
                    "tau_months": tau,
                    "difference_months": diff,
                    "se": se,
                    "ci_lower": diff - 1.959963984540054 * se,
                    "ci_upper": diff + 1.959963984540054 * se,
                    "z": zstat,
                    "p_value": p,
                }
            )
    return pd.DataFrame(records)


def age_standardized_trend(cell_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    year_col = "Year of diagnosis_label"
    age_col = "Age recode with <1 year olds and 90+_label"
    cells = cell_table[cell_table["tau_months"] == 60].copy()
    cells["year"] = pd.to_numeric(cells[year_col], errors="coerce")
    cells = cells[(cells["year"] >= 2004) & (cells["year"] <= 2018)]
    cells = cells[~cells[age_col].str.contains("Unknown", na=False)]
    age_totals = cells.groupby(age_col)["n"].sum()
    weights = (age_totals / age_totals.sum()).rename("standard_weight")
    cells = cells.join(weights, on=age_col)
    records = []
    for year, part in cells.groupby("year"):
        if len(part) != len(weights):
            continue
        w = part["standard_weight"].to_numpy(float)
        obs = float(np.dot(w, part["observed"]))
        exp = float(np.dot(w, part["expected"]))
        rlel = float(np.dot(w, part["rlel_months"]))
        variance = float(np.dot(w**2, part["rlel_se"].to_numpy(float) ** 2))
        se = math.sqrt(variance)
        records.append(
            {
                "year": int(year),
                "standardized_observed_5y": obs,
                "standardized_expected_5y": exp,
                "standardized_relative_5y": obs / exp,
                "standardized_rlel_5y_months": rlel,
                "se": se,
                "ci_lower": rlel - 1.959963984540054 * se,
                "ci_upper": rlel + 1.959963984540054 * se,
            }
        )
    trend = pd.DataFrame(records).sort_values("year")
    x = sm.add_constant(trend["year"].astype(float))
    fit_wls = sm.WLS(
        trend["standardized_rlel_5y_months"], x, weights=1 / trend["se"].clip(lower=1e-8) ** 2
    ).fit()
    fit_hc3 = sm.OLS(trend["standardized_rlel_5y_months"], x).fit(cov_type="HC3")
    first, last = trend.iloc[0], trend.iloc[-1]
    first_last_se = math.sqrt(float(first["se"] ** 2 + last["se"] ** 2))
    stats = pd.DataFrame(
        [
            {
                "model": "inverse-variance WLS",
                "slope_months_per_year": fit_wls.params["year"],
                "se": fit_wls.bse["year"],
                "ci_lower": fit_wls.conf_int().loc["year", 0],
                "ci_upper": fit_wls.conf_int().loc["year", 1],
                "p_value": fit_wls.pvalues["year"],
            },
            {
                "model": "OLS with HC3 robust SE",
                "slope_months_per_year": fit_hc3.params["year"],
                "se": fit_hc3.bse["year"],
                "ci_lower": fit_hc3.conf_int().loc["year", 0],
                "ci_upper": fit_hc3.conf_int().loc["year", 1],
                "p_value": fit_hc3.pvalues["year"],
            },
            {
                "model": "2018 minus 2004",
                "slope_months_per_year": last["standardized_rlel_5y_months"]
                - first["standardized_rlel_5y_months"],
                "se": first_last_se,
                "ci_lower": last["standardized_rlel_5y_months"]
                - first["standardized_rlel_5y_months"]
                - 1.959963984540054 * first_last_se,
                "ci_upper": last["standardized_rlel_5y_months"]
                - first["standardized_rlel_5y_months"]
                + 1.959963984540054 * first_last_se,
                "p_value": math.erfc(
                    abs(
                        (last["standardized_rlel_5y_months"] - first["standardized_rlel_5y_months"])
                        / first_last_se
                    )
                    / math.sqrt(2)
                ),
            },
        ]
    )
    weights.reset_index().to_csv(OUT / "standard_age_weights.csv", index=False)
    return trend, stats


def case_listing_audit() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = SEER_DIR / "nonmalignant_CNS_case_listing_2004_2023.txt"
    usecols = [
        "Year of diagnosis",
        "Age recode with single ages and 90+",
        "Sex",
        "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)",
        "Primary Site - labeled",
        "Behavior code ICD-O-3",
        "SEER Brain and CNS Recode",
        "Diagnostic Confirmation",
        "Type of Reporting Source",
        "Survival months",
        "Vital status recode (study cutoff used)",
    ]
    cohort_counts = {"all_2004_2023": 0, "5y_2004_2018": 0, "10y_2004_2013": 0, "15y_2004_2008": 0}
    cohort_counts_excl = {k: 0 for k in cohort_counts}
    table1_parts = []
    missing_counts = {c: 0 for c in usecols}
    total = 0
    for chunk in pd.read_csv(path, sep="\t", quotechar='"', usecols=usecols, dtype=str, chunksize=75_000):
        total += len(chunk)
        year = pd.to_numeric(chunk["Year of diagnosis"], errors="coerce")
        surv = pd.to_numeric(chunk["Survival months"], errors="coerce")
        reporting = chunk["Type of Reporting Source"].fillna("")
        vital = chunk["Vital status recode (study cutoff used)"].fillna("")
        excluded = reporting.str.contains("Death certificate|Autopsy", case=False, regex=True) | (
            (vital == "Alive") & (surv.fillna(-1) <= 0)
        )
        masks = {
            "all_2004_2023": year.between(2004, 2023),
            "5y_2004_2018": year.between(2004, 2018),
            "10y_2004_2013": year.between(2004, 2013),
            "15y_2004_2008": year.between(2004, 2008),
        }
        for key, mask in masks.items():
            cohort_counts[key] += int(mask.sum())
            cohort_counts_excl[key] += int((mask & ~excluded).sum())
        main = chunk[masks["10y_2004_2013"] & ~excluded].copy()
        table1_parts.append(main)
        for col in usecols:
            missing_counts[col] += int(chunk[col].isna().sum() + chunk[col].isin(["", "Blank(s)", "Unknown"]).sum())

    main = pd.concat(table1_parts, ignore_index=True)
    age = pd.to_numeric(main["Age recode with single ages and 90+"].str.extract(r"(\d+)")[0], errors="coerce")
    table1_records = [
        {"characteristic": "N", "level": "10-year analytic cohort", "n": len(main), "percent": 100.0},
        {"characteristic": "Age", "level": "Median (IQR)", "n": np.nan, "percent": np.nan,
         "value": f"{age.median():.0f} ({age.quantile(.25):.0f}-{age.quantile(.75):.0f})"},
    ]
    for variable in [
        "Sex",
        "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)",
        "Behavior code ICD-O-3",
        "Diagnostic Confirmation",
        "SEER Brain and CNS Recode",
    ]:
        counts = main[variable].fillna("Missing").value_counts(dropna=False)
        for level, count in counts.items():
            table1_records.append(
                {
                    "characteristic": variable,
                    "level": level,
                    "n": int(count),
                    "percent": 100 * count / len(main),
                }
            )
    cohort = pd.DataFrame(
        [
            {
                "cohort": k,
                "case_listing_count_before_standard_exclusions": cohort_counts[k],
                "case_listing_count_after_reconstructed_exclusions": cohort_counts_excl[k],
            }
            for k in cohort_counts
        ]
    )
    missing = pd.DataFrame(
        [
            {"variable": col, "missing_or_unknown_n": count, "percent_of_all_listing": 100 * count / total}
            for col, count in missing_counts.items()
        ]
    )
    return cohort, pd.DataFrame(table1_records), missing


def main() -> None:
    parsed = {key: parse_export(base) for key, base in EXPORTS.items()}
    tables: dict[str, pd.DataFrame] = {}
    tables["overall_5y"] = rlel_table(parsed["overall_5y"], [12, 60])
    tables["overall_10y"] = rlel_table(parsed["overall_10y"], [12, 60, 120])
    tables["overall_15y"] = rlel_table(parsed["overall_15y"], [12, 60, 120, 180])
    for key in ["behavior_10y", "brain_recode_10y", "c75_site_10y", "sex_10y", "age_10y"]:
        tables[key] = rlel_table(parsed[key], [60, 120])
    tables["dxyear_5y"] = rlel_table(parsed["dxyear_5y"], [60])
    tables["dxyear_age_5y"] = rlel_table(parsed["dxyear_age_5y"], [60])

    for name, table in tables.items():
        table.to_csv(OUT / f"{name}_rlel.csv", index=False)

    behavior_diffs = add_pairwise_differences(tables["behavior_10y"], "Behavior code ICD-O-3_label", 120)
    sex_diffs = add_pairwise_differences(
        tables["sex_10y"][tables["sex_10y"]["Sex"] != "0"], "Sex_label", 120
    )
    behavior_diffs.to_csv(OUT / "behavior_pairwise_10y.csv", index=False)
    sex_diffs.to_csv(OUT / "sex_pairwise_10y.csv", index=False)

    trend, trend_stats = age_standardized_trend(tables["dxyear_age_5y"])
    trend.to_csv(OUT / "age_standardized_dxyear_5y.csv", index=False)
    trend_stats.to_csv(OUT / "age_standardized_trend_tests.csv", index=False)

    cohort, table1, missing = case_listing_audit()
    cohort.to_csv(OUT / "cohort_reconciliation.csv", index=False)
    table1.to_csv(OUT / "table1_10y_cohort.csv", index=False)
    missing.to_csv(OUT / "case_listing_missingness.csv", index=False)

    # Exact descriptive composition of the 10-year survival cohort from
    # mutually exclusive SEER*Stat survival pages. This avoids presenting the
    # larger case-listing export as if it were the final survival denominator.
    analytic_table1 = []
    total_n = int(tables["overall_10y"].iloc[0]["n"])
    analytic_table1.append(
        {"characteristic": "Cohort", "level": "Overall", "n": total_n, "percent": 100.0}
    )
    for key, code_col, label_col, excluded_codes in [
        ("sex_10y", "Sex", "Sex_label", {"0"}),
        ("behavior_10y", "Behavior code ICD-O-3", "Behavior code ICD-O-3_label", set()),
        (
            "age_10y",
            "Age recode with <1 year olds and 90+",
            "Age recode with <1 year olds and 90+_label",
            set(),
        ),
    ]:
        part = tables[key]
        part = part[part["tau_months"] == 120].drop_duplicates(code_col)
        part = part[~part[code_col].astype(str).isin(excluded_codes)]
        for _, row in part.iterrows():
            analytic_table1.append(
                {
                    "characteristic": code_col,
                    "level": row[label_col],
                    "n": int(row["n"]),
                    "percent": 100 * int(row["n"]) / total_n,
                }
            )
    pd.DataFrame(analytic_table1).to_csv(OUT / "table1_exact_10y_survival_cohort.csv", index=False)

    # Figure source data. Drawing is performed exclusively in R.
    overall_curve = numeric(
        parsed["overall_15y"].life[parsed["overall_15y"].life["Page type"] == "0"],
        [
            "Interval",
            "Observed Survival (Cum)",
            "Expected Survival (Cum)",
            "Observed Lower CI (Cum)",
            "Observed Upper CI (Cum)",
        ],
    )
    overall_curve = overall_curve[overall_curve["Interval"].between(0, 180)][
        [
            "Interval",
            "Observed Survival (Cum)",
            "Expected Survival (Cum)",
            "Observed Lower CI (Cum)",
            "Observed Upper CI (Cum)",
        ]
    ]
    overall_curve.to_csv(OUT / "figure1_overall_survival_curve_source.csv", index=False)

    heterogeneity_parts = []
    behavior = tables["behavior_10y"]
    behavior = behavior[behavior["tau_months"] == 120].copy()
    behavior["domain"] = "Behavior"
    behavior["group"] = behavior["Behavior code ICD-O-3_label"]
    heterogeneity_parts.append(behavior)
    brain = tables["brain_recode_10y"]
    brain["_label_clean"] = brain["SEER Brain and CNS Recode_label"].str.strip()
    keep_brain = [
        "2.1 Meningiomas",
        "2.2 Tumors of the cranial and paraspinal nerves",
        "2.3 Ependymal tumors",
        "2.4 Choroid plexus tumors",
        "2.5 Neuronal and mixed neuronal-glial tumors",
        "2.6 Mesenchymal, nonmeningothelial tumors",
        "2.7 Other astrocytic tumors",
        "2.9 Other non-malignant Brain/ONS",
    ]
    brain = brain[(brain["tau_months"] == 120) & brain["_label_clean"].isin(keep_brain)].copy()
    brain["domain"] = "Brain/CNS recode"
    brain["group"] = brain["_label_clean"].str.replace(r"^2\.\d+\s+", "", regex=True)
    heterogeneity_parts.append(brain)
    c75 = tables["c75_site_10y"]
    c75 = c75[c75["tau_months"] == 120].copy()
    c75["domain"] = "C75 site"
    c75["group"] = c75["Primary Site - labeled_label"].str.replace(r"^C75\.\d-", "", regex=True)
    heterogeneity_parts.append(c75)
    heterogeneity = pd.concat(heterogeneity_parts, ignore_index=True)
    heterogeneity[
        [
            "domain", "group", "n", "relative", "observed_rmst_months",
            "expected_rmst_months", "rlel_months", "relative_rlel_percent",
            "rlel_ci_lower", "rlel_ci_upper"
        ]
    ].to_csv(OUT / "figure2_heterogeneity_source.csv", index=False)

    manifest = []
    for base in list(EXPORTS.values()) + POHAR_PERME_EXPORTS:
        for suffix in (".txt", ".dic", ".ss"):
            path = SEER_DIR / f"{base}{suffix}"
            if path.exists():
                manifest.append(
                    {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                )
    for listing_name in (
        "nonmalignant_CNS_case_listing_2004_2023.txt",
        "nonmalignant_CNS_case_listing_2004_2023_SES_treatment.txt",
    ):
        listing = SEER_DIR / listing_name
        if listing.exists():
            manifest.append(
                {"file": listing.name, "bytes": listing.stat().st_size, "sha256": sha256(listing)}
            )
    pd.DataFrame(manifest).to_csv(OUT / "source_file_manifest.csv", index=False)

    key_results = {
        "method": "Trapezoidal integration of expected minus observed monthly cumulative survival; relative RLEL is RLEL divided by expected restricted mean survival through the same horizon; Greenwood covariance approximation; expected survival treated as fixed.",
        "overall": {
            name: table.to_dict(orient="records") for name, table in tables.items() if name.startswith("overall")
        },
        "age_standardized_trend": trend.to_dict(orient="records"),
        "trend_tests": trend_stats.to_dict(orient="records"),
        "cohort_reconciliation": cohort.to_dict(orient="records"),
    }
    (OUT / "key_results.json").write_text(json.dumps(key_results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output": str(OUT), "tables": len(tables), "status": "ok"}, indent=2))


if __name__ == "__main__":
    main()
