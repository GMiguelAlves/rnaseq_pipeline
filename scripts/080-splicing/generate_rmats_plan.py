#!/usr/bin/env python3

import argparse
import csv
import gzip
import itertools
import re
from collections import Counter
from pathlib import Path


UNKNOWN_VALUES = {"", "na", "nan", "none", "unknown", "not_available"}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate an rMATS splicing analysis plan.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--star-root", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--projects", default="auto")
    parser.add_argument("--include-all", action="store_true")
    parser.add_argument("--test-variables", default="condition,stage,sex,tissue,infection_mode")
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args()


def sanitize(value):
    value = str(value or "unknown")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.strip("_")
    return value or "unknown"


def split_csv(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, newline="", encoding="utf-8-sig")


def detect_delimiter(path):
    with open_text(path) as handle:
        first = ""
        for line in handle:
            if line.strip():
                first = line
                break
    return "\t" if first.count("\t") > first.count(",") else ","


def read_metadata(path):
    delimiter = detect_delimiter(path)
    with open_text(path) as handle:
        lines = [line for line in handle if line.strip()]
    return list(csv.DictReader(lines, delimiter=delimiter))


def normalized_level(value):
    value = str(value or "").strip()
    if value.lower() in UNKNOWN_VALUES:
        return ""
    return value


def bam_path(star_root, row):
    return Path(star_root) / row["dataset"] / row["sample_id"] / "Aligned.sortedByCoord.out.bam"


def write_bam_list(path, bams):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(",".join(str(bam) for bam in bams) + "\n", encoding="utf-8")


def add_comparison(rows, args, sample_rows, scope, project, variable, level_a, level_b):
    selected_a = [row for row in sample_rows if normalized_level(row.get(variable)) == level_a]
    selected_b = [row for row in sample_rows if normalized_level(row.get(variable)) == level_b]

    bams_a = [(row, bam_path(args.star_root, row)) for row in selected_a]
    bams_b = [(row, bam_path(args.star_root, row)) for row in selected_b]
    missing = [path for _, path in bams_a + bams_b if not path.exists()]

    if missing and not args.allow_missing:
        shown = ", ".join(str(path) for path in missing[:20])
        raise SystemExit(
            f"[ERRO] BAMs ausentes para {scope}/{variable}/{level_a}_vs_{level_b}: {shown}. "
            "Use --allow-missing apenas se quiser usar o subconjunto existente."
        )

    if args.allow_missing:
        bams_a = [(row, path) for row, path in bams_a if path.exists()]
        bams_b = [(row, path) for row, path in bams_b if path.exists()]

    if len(bams_a) < args.min_replicates or len(bams_b) < args.min_replicates:
        return

    scope_dir = project if scope == "project" else "all_projects"
    contrast = f"{sanitize(level_a)}_vs_{sanitize(level_b)}"
    analysis_id = f"{sanitize(scope_dir)}_{sanitize(variable)}_{contrast}"
    out_dir = Path(args.output_root) / scope_dir / sanitize(variable) / contrast
    tmp_dir = out_dir / "tmp"
    input_dir = Path(args.output_root) / scope_dir / "inputs"
    b1_file = input_dir / f"{analysis_id}_b1.txt"
    b2_file = input_dir / f"{analysis_id}_b2.txt"

    bams_a_paths = [path for _, path in bams_a]
    bams_b_paths = [path for _, path in bams_b]
    write_bam_list(b1_file, bams_a_paths)
    write_bam_list(b2_file, bams_b_paths)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    rows.append(
        {
            "analysis_id": analysis_id,
            "scope": scope,
            "project": project,
            "variable": variable,
            "level_a": level_a,
            "level_b": level_b,
            "n_level_a": len(bams_a_paths),
            "n_level_b": len(bams_b_paths),
            "b1_file": str(b1_file),
            "b2_file": str(b2_file),
            "gtf": str(args.gtf),
            "output_dir": str(out_dir),
            "tmp_dir": str(tmp_dir),
            "bams_level_a": ",".join(str(path) for path in bams_a_paths),
            "bams_level_b": ",".join(str(path) for path in bams_b_paths),
        }
    )


def add_scope(rows, args, metadata, scope, project, variables):
    if scope == "project":
        sample_rows = [row for row in metadata if row.get("dataset") == project]
    else:
        sample_rows = list(metadata)

    if not sample_rows:
        return

    for variable in variables:
        if variable not in sample_rows[0]:
            continue
        levels = [normalized_level(row.get(variable)) for row in sample_rows]
        counts = Counter(level for level in levels if level)
        valid_levels = sorted(level for level, count in counts.items() if count >= args.min_replicates)
        for level_a, level_b in itertools.combinations(valid_levels, 2):
            add_comparison(rows, args, sample_rows, scope, project, variable, level_a, level_b)


def main():
    args = parse_args()
    metadata = read_metadata(args.metadata)
    required = {"dataset", "sample_id"}
    if not metadata:
        raise SystemExit("[ERRO] Metadata vazio.")
    missing_cols = sorted(required - set(metadata[0]))
    if missing_cols:
        raise SystemExit(f"[ERRO] Metadata sem colunas obrigatorias: {', '.join(missing_cols)}")

    metadata = [
        row for row in metadata
        if row.get("dataset") and row.get("sample_id")
    ]

    if args.projects == "auto":
        projects = sorted({row["dataset"] for row in metadata})
    else:
        projects = split_csv(args.projects)

    variables = split_csv(args.test_variables)
    rows = []
    for project in projects:
        add_scope(rows, args, metadata, "project", project, variables)
    if args.include_all:
        add_scope(rows, args, metadata, "all_projects", "all_projects", variables)

    if not rows:
        raise SystemExit("[ERRO] Nenhuma comparacao entrou no plano. Verifique metadata, BAMs, variaveis e replicatas.")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Plano rMATS: {output}")
    print(f"[OK] Analises: {len(rows)}")


if __name__ == "__main__":
    main()
