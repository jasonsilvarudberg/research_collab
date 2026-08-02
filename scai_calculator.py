# Assign SCAI cardiogenic shock stage using data from the first 24 hours
# after ICU admission.
#
# Jason Silva-Rudberg (Limited English Proficiency and End-of-Life Decision-Making Among Patients with Cardiogenic Shock in the Cardiac Intensive Care Unit)
#
# Usage:
#   python scai_calculator.py -i input.xlsx -o scai_output.xlsx
#   python scai_calculator.py -i input.csv  -o scai_output.csv
#
# The subject ID column must be named: subj
#
# Input variable definitions (must be in input.xlsx or input.csv)
# sbp<60: SBP <60 mmHg (0/1) in first 24 hr
# map<50: MAP <50 mmHg (0/1) in first 24 hr
# sbp<90: SBP <90 mmHg (0/1) in first 24 hr
# map<65: MAP <65 mmHg (0/1) in first 24 hr
# max_lactate_24hr: highest lactate in first 24 hr
# max_alt_24hr: highest ALT in first 24hr
# lowest_pH: lowest arterial pH in first 24 hr
# impella: Impella present (0/1) in first 24 hr
# ecmo: ECMO present (0/1) in first 24 hr
# iabp: IABP present (0/1) in first 24 hr
# 24hr_n_vasoactives: maximum number of concurrent vasopressors and inotropes in  first 24 hr

# ── SCAI staging (hierarchical; first True condition wins) ────────────────
#
# Limitations
#   1. Out-of-hospital cardiac arrest (a Stage E OR-criterion) is not available
#      in our dataset and is therefore not captured.
#   2. Stage D path 1 requires BOTH confirmed hypotension AND D-level hypoperfusion
#      (lactate >5–10 or ALT >500). Patients with severe isolated hypoperfusion
#      but missing/absent BP data will not meet this path and may be
#      underclassified. Review Stage A cases with high lactate manually.
#   3. 24-hour worst-value approach may upstage patients with transient
#      abnormalities.

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


NUMERIC_COLS = [
    "sbp<60",
    "map<50",
    "sbp<90",
    "map<65",
    "max_lactate_24hr",
    "max_alt_24hr",
    "lowest_pH",
    "impella",
    "ecmo",
    "iabp",
    "24hr_n_vasoactives",
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Assign SCAI cardiogenic shock stage."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input Excel or CSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the output Excel or CSV file.",
    )
    return parser.parse_args()


def read_input_file(file_path):
    extension = file_path.suffix.lower()

    if extension == ".csv":
        return pd.read_csv(file_path)

    if extension in {".xlsx", ".xls"}:
        return pd.read_excel(file_path)

    raise ValueError(
        "Unsupported input format. Input must be a .csv, .xlsx, or .xls file."
    )


def save_output_file(df, file_path):
    extension = file_path.suffix.lower()

    file_path.parent.mkdir(parents=True, exist_ok=True)

    if extension == ".csv":
        df.to_csv(file_path, index=False)
        return

    if extension in {".xlsx", ".xls"}:
        df.to_excel(file_path, index=False)
        return

    raise ValueError(
        "Unsupported output format. Output must be a .csv, .xlsx, or .xls file."
    )


def assign_scai_stage(df):
    if "subj" not in df.columns:
        raise ValueError("Subject ID column not found. It must be named 'subj'.")

    # Convert relevant columns to numeric. Nonnumeric values become missing.
    for col in NUMERIC_COLS:
        if col not in df.columns:
            print(f"Warning: '{col}' not found; treating it as missing.")
            df[col] = np.nan
            continue

        missing_before = df[col].isna().sum()
        df[col] = pd.to_numeric(df[col], errors="coerce")
        missing_after = df[col].isna().sum()

        newly_missing = missing_after - missing_before
        if newly_missing > 0:
            print(
                f"Warning: '{col}' had {newly_missing} nonnumeric "
                "value(s) converted to missing."
            )

    ################# Hemodynamic criteria###################3
    sbp_60_90 = (df["sbp<90"] == 1) & (df["sbp<60"] == 0)
    map_50_65 = (df["map<65"] == 1) & (df["map<50"] == 0)

    hypo_bcd = sbp_60_90 | map_50_65
    hypo_e = (df["sbp<60"] == 1) | (df["map<50"] == 1)

    # Perfusion criteria
    lactate = df["max_lactate_24hr"]
    alt = df["max_alt_24hr"]
    ph = df["lowest_pH"]

    perf_bc = (
        ((lactate >= 2) & (lactate <= 5))
        | ((alt >= 200) & (alt <= 500))
    )

    perf_d = (
        ((lactate > 5) & (lactate <= 10))
        | (alt > 500)
    )

    perf_e = (lactate > 10) | (ph < 7.2)

    # #########Treatment intensity ################
    # ############Missing device indicators and vasoactive counts are assumed to be zero. ##############
    n_devices = (
        df[["impella", "ecmo", "iabp"]]
        .fillna(0)
        .sum(axis=1)
    )
    n_drugs = df["24hr_n_vasoactives"].fillna(0)
    n_total = n_devices + n_drugs

    tx_0 = n_total == 0
    tx_1 = n_total == 1
    tx_d = (n_total >= 2) & (n_total <= 5)
    tx_e = (n_drugs >= 3) | (n_devices >= 3)

    # SCAI staging is hierarchical; the first true condition is assigned.

    # Stage E: extremis or refractory shock
    cond_e = hypo_e | perf_e | tx_e

    # Stage D: failure to stabilize with initial therapy
    cond_d = (
        (hypo_bcd & perf_d)
        | tx_d
        | (tx_1 & (hypo_bcd | perf_bc | perf_d))
    )

    # Stage C: hypoperfusion consistent with shock
    cond_c = (
        (hypo_bcd & perf_bc & tx_0)
        | (tx_1 & ~hypo_bcd & ~perf_bc & ~perf_d)
    )

    # Stage B: hemodynamic or perfusion abnormality without treatment
    cond_b = (hypo_bcd | perf_bc) & tx_0

    df["scai"] = np.select(
        [cond_e, cond_d, cond_c, cond_b],
        ["E", "D", "C", "B"],
        default="A",
    )

    return df[["subj", "scai"]].copy()


def main():
    args = parse_arguments()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = read_input_file(input_path)
    print(f"Loaded {len(df):,} rows × {df.shape[1]} columns")

    output = assign_scai_stage(df)
    save_output_file(output, output_path)

    print(f"\nSaved to: {output_path}")
    print(f"N = {len(output):,}")

    print("\nSCAI stage distribution:")
    distribution = output["scai"].value_counts().sort_index()

    for stage in ["A", "B", "C", "D", "E"]:
        n = int(distribution.get(stage, 0))
        percent = 100 * n / len(output) if len(output) else 0
        print(f"  Stage {stage}: {n:>5,} ({percent:.1f}%)")


if __name__ == "__main__":
    main()
