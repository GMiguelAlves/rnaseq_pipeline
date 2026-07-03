#!/usr/bin/env python3
"""Validate the final RNA-seq metadata contract and common analysis fields."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {"dataset", "sample_id", "run_accession"}
DEFAULT_RECOMMENDED_COLUMNS = [
    "condition",
    "stage",
    "tissue",
    "sex",
    "batch",
    "replicate",
]
UNKNOWN_VALUES = {"", "unknown", "na", "n/a", "none", "null", "nan", "NA"}


def read_table(path: Path) -> pd.DataFrame:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    sep = "\t" if ".tsv" in suffixes else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def is_unknown(value: object) -> bool:
    return str(value).strip() in UNKNOWN_VALUES


def warn(message: str) -> None:
    print(f"[WARN] {message}")


def fail_or_warn(message: str, strict: bool, failures: list[str]) -> None:
    if strict:
        failures.append(message)
    else:
        warn(message)


def usable_levels(series: pd.Series) -> list[str]:
    values = sorted({str(value).strip() for value in series if not is_unknown(value)})
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--allow-missing", action="store_true")
    parser.add_argument(
        "--recommended-columns",
        default=",".join(DEFAULT_RECOMMENDED_COLUMNS),
        help="Comma-separated columns checked when present in downstream analyses.",
    )
    parser.add_argument(
        "--min-samples-per-project",
        type=int,
        default=2,
        help="Warn when a project has fewer samples than this.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat recommended-field warnings as errors.",
    )
    args = parser.parse_args()

    if not args.metadata.exists():
        if args.allow_missing:
            print(f"[WARN] Metadata not found yet: {args.metadata}")
            return
        raise FileNotFoundError(f"Metadata not found: {args.metadata}")

    metadata = read_table(args.metadata)
    failures: list[str] = []
    missing = sorted(REQUIRED_COLUMNS - set(metadata.columns))
    if missing:
        raise SystemExit("[ERROR] Missing required columns: " + ", ".join(missing))

    empty = [
        column
        for column in sorted(REQUIRED_COLUMNS)
        if metadata[column].astype(str).str.strip().eq("").any()
    ]
    if empty:
        raise SystemExit("[ERROR] Empty values in required columns: " + ", ".join(empty))

    duplicated = metadata["run_accession"].duplicated()
    if duplicated.any():
        examples = ", ".join(metadata.loc[duplicated, "run_accession"].head(10))
        raise SystemExit(f"[ERROR] Duplicated run_accession values, e.g. {examples}")

    recommended_columns = split_csv(args.recommended_columns)
    missing_recommended = [column for column in recommended_columns if column not in metadata.columns]
    if missing_recommended:
        fail_or_warn(
            "Recommended columns missing for downstream analyses: "
            + ", ".join(missing_recommended)
            + ". Add them in the project YAML or remove them from DEG/WGCNA/Mfuzz settings.",
            args.strict,
            failures,
        )

    for column in recommended_columns:
        if column not in metadata.columns:
            continue
        unknown_count = metadata[column].map(is_unknown).sum()
        if unknown_count == len(metadata):
            fail_or_warn(
                f"Column '{column}' is entirely unknown/empty; comparisons using it will be skipped.",
                args.strict,
                failures,
            )
        elif unknown_count:
            warn(f"Column '{column}' has {unknown_count}/{len(metadata)} unknown or empty values.")

    project_counts = metadata.groupby("dataset")["sample_id"].nunique().sort_index()
    small_projects = project_counts[project_counts < args.min_samples_per_project]
    if not small_projects.empty:
        fail_or_warn(
            "Projects with fewer than "
            f"{args.min_samples_per_project} samples: "
            + ", ".join(f"{project}={count}" for project, count in small_projects.items()),
            args.strict,
            failures,
        )

    if "replicate" in metadata.columns:
        missing_replicate = metadata["replicate"].map(is_unknown).sum()
        if missing_replicate:
            warn(
                f"Column 'replicate' has {missing_replicate}/{len(metadata)} unknown or empty values; "
                "DEG/DTU/splicing plans rely on sample counts per group, but explicit replicate labels help auditing."
            )

    usable_variables: list[str] = []
    for column in recommended_columns:
        if column not in metadata.columns:
            continue
        levels = usable_levels(metadata[column])
        if len(levels) >= 2:
            usable_variables.append(f"{column}({len(levels)} levels)")

    if not usable_variables:
        fail_or_warn(
            "No recommended biological column has at least two non-unknown levels. "
            "DEG/DTU/splicing comparisons will probably produce no plan.",
            args.strict,
            failures,
        )

    if failures:
        raise SystemExit("[ERROR] Metadata validation failed:\n- " + "\n- ".join(failures))

    print(f"[OK] Metadata validated: {args.metadata}")
    print(f"[OK] Rows: {len(metadata)}")
    print(f"[OK] Projects: {metadata['dataset'].nunique()}")
    print(f"[OK] Samples: {metadata['sample_id'].nunique()}")
    print("[OK] Samples per project:")
    for project, count in project_counts.items():
        print(f"  - {project}: {count}")
    if usable_variables:
        print("[OK] Usable metadata variables: " + ", ".join(usable_variables))
    print("[OK] Tip: run with --strict before production submissions.")


if __name__ == "__main__":
    main()
