#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


EVENT_FILES = {
    "SE": "SE.MATS.JC.txt",
    "A5SS": "A5SS.MATS.JC.txt",
    "A3SS": "A3SS.MATS.JC.txt",
    "MXE": "MXE.MATS.JC.txt",
    "RI": "RI.MATS.JC.txt",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize rMATS event tables.")
    parser.add_argument("--rmats-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--significant-output", default="")
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def numeric_series(df, column):
    if column not in df.columns:
        return pd.Series(dtype="float64")
    return pd.to_numeric(df[column], errors="coerce")


def main():
    args = parse_args()
    rmats_dir = Path(args.rmats_dir)
    rows = []
    significant_tables = []

    for event_type, filename in EVENT_FILES.items():
        path = rmats_dir / filename
        if not path.exists():
            rows.append(
                {
                    "event_type": event_type,
                    "file": str(path),
                    "file_found": 0,
                    "total_events": 0,
                    "significant_events": 0,
                    "min_fdr": "",
                    "median_abs_inc_level_difference": "",
                    "max_abs_inc_level_difference": "",
                }
            )
            continue

        df = pd.read_csv(path, sep="\t")
        fdr = numeric_series(df, "FDR")
        sig_mask = fdr <= args.fdr
        inc_diff = numeric_series(df, "IncLevelDifference").abs()
        sig_df = df.loc[sig_mask].copy()
        if not sig_df.empty:
            sig_df.insert(0, "event_type", event_type)
            significant_tables.append(sig_df)

        rows.append(
            {
                "event_type": event_type,
                "file": str(path),
                "file_found": 1,
                "total_events": int(len(df)),
                "significant_events": int(sig_mask.sum()),
                "min_fdr": "" if fdr.dropna().empty else float(fdr.min()),
                "median_abs_inc_level_difference": "" if inc_diff.dropna().empty else float(inc_diff.median()),
                "max_abs_inc_level_difference": "" if inc_diff.dropna().empty else float(inc_diff.max()),
            }
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False, compression="infer")

    if args.significant_output:
        sig_output = Path(args.significant_output)
        sig_output.parent.mkdir(parents=True, exist_ok=True)
        if significant_tables:
            pd.concat(significant_tables, ignore_index=True).to_csv(sig_output, sep="\t", index=False, compression="infer")
        else:
            pd.DataFrame(columns=["event_type"]).to_csv(sig_output, sep="\t", index=False, compression="infer")

    print(f"[OK] Resumo rMATS: {output}")


if __name__ == "__main__":
    main()
