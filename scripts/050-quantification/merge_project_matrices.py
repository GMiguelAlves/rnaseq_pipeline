#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge project-level quantification matrices into all-project outputs."
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--projects", default="auto", help="Comma-separated projects or auto for existing project files.")
    parser.add_argument("--method", choices=["salmon", "star"], default="salmon")
    parser.add_argument("--counts-name", required=True)
    parser.add_argument("--expression-name", required=True)
    parser.add_argument("--sample-table-name", required=True)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args()


def strip_table_suffix(name: str) -> str:
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".tsv"):
        name = name[:-4]
    return name


def table_candidates(path: Path) -> list[Path]:
    candidates = [path]
    text = str(path)
    if text.endswith(".tsv.gz"):
        candidates.append(Path(text[:-3]))
    elif text.endswith(".tsv"):
        candidates.append(Path(text + ".gz"))
    return candidates


def first_existing(path: Path) -> Path | None:
    for candidate in table_candidates(path):
        if candidate.is_file():
            return candidate
    return None


def discover_projects(output_dir: Path, method: str) -> list[str]:
    suffix = "_star_cpm_matrix" if method == "star" else "_tpm_matrix"
    projects: set[str] = set()
    for path in output_dir.glob("*_counts_matrix.tsv*"):
        name = strip_table_suffix(path.name)
        if name.endswith("_counts_matrix"):
            projects.add(name[: -len("_counts_matrix")])
    for path in output_dir.glob(f"*{suffix}.tsv*"):
        name = strip_table_suffix(path.name)
        if name.endswith(suffix):
            projects.add(name[: -len(suffix)])
    return sorted(projects)


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.shape[1] < 2:
        raise ValueError(f"Matrix has fewer than two columns: {path}")
    first = df.columns[0]
    if first != "gene_id":
        df = df.rename(columns={first: "gene_id"})
    return df


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def project_paths(output_dir: Path, project: str, method: str) -> tuple[Path | None, Path | None, Path | None]:
    counts = first_existing(output_dir / f"{project}_counts_matrix.tsv.gz")
    expression_suffix = "star_cpm_matrix" if method == "star" else "tpm_matrix"
    expression = first_existing(output_dir / f"{project}_{expression_suffix}.tsv.gz")
    samples = first_existing(output_dir / f"{project}_quant_samples.tsv.gz")
    return counts, expression, samples


def combined_import_id(project: str, sample_table: pd.DataFrame, column: str) -> str:
    if column.startswith(project + "__"):
        return column
    lookup_cols = []
    if "import_id" in sample_table.columns:
        lookup_cols.append("import_id")
    if "sample_id" in sample_table.columns:
        lookup_cols.append("sample_id")
    for col in lookup_cols:
        matched = sample_table.loc[sample_table[col] == column]
        if not matched.empty:
            row = matched.iloc[0]
            dataset = row.get("dataset", project) or project
            sample_id = row.get("sample_id", column) or column
            return f"{dataset}__{sample_id}"
    return f"{project}__{column}"


def prepare_matrix(df: pd.DataFrame, project: str, samples: pd.DataFrame) -> pd.DataFrame:
    rename = {
        col: combined_import_id(project, samples, col)
        for col in df.columns
        if col != "gene_id"
    }
    return df.rename(columns=rename)


def prepare_samples(project: str, samples: pd.DataFrame) -> pd.DataFrame:
    samples = samples.copy()
    if "dataset" not in samples.columns:
        samples["dataset"] = project
    samples["dataset"] = samples["dataset"].replace("", project)
    if "sample_id" not in samples.columns:
        if "import_id" in samples.columns:
            samples["sample_id"] = samples["import_id"]
        else:
            raise ValueError(f"Sample table for {project} needs sample_id or import_id.")
    samples["import_id"] = samples.apply(
        lambda row: row["import_id"]
        if str(row.get("import_id", "")).startswith(str(row["dataset"]) + "__")
        else f"{row['dataset']}__{row['sample_id']}",
        axis=1,
    )
    return samples


def merge_matrices(matrices: list[pd.DataFrame]) -> pd.DataFrame:
    merged = matrices[0]
    for matrix in matrices[1:]:
        merged = merged.merge(matrix, on="gene_id", how="outer")
    value_cols = [col for col in merged.columns if col != "gene_id"]
    merged[value_cols] = merged[value_cols].fillna("0")
    return merged.sort_values("gene_id")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    if args.projects == "auto":
        projects = discover_projects(output_dir, args.method)
    else:
        projects = [p.strip() for p in args.projects.split(",") if p.strip()]

    count_matrices: list[pd.DataFrame] = []
    expression_matrices: list[pd.DataFrame] = []
    sample_tables: list[pd.DataFrame] = []
    missing: list[str] = []

    for project in projects:
        counts_path, expression_path, samples_path = project_paths(output_dir, project, args.method)
        if counts_path is None or expression_path is None or samples_path is None:
            missing.append(project)
            continue
        samples = prepare_samples(project, read_table(samples_path))
        count_matrices.append(prepare_matrix(read_matrix(counts_path), project, samples))
        expression_matrices.append(prepare_matrix(read_matrix(expression_path), project, samples))
        sample_tables.append(samples)

    if missing and not args.allow_missing:
        raise SystemExit("[ERRO] Project quantification outputs missing for: " + ", ".join(missing))
    if not count_matrices:
        raise SystemExit("[ERRO] No project quantification outputs found to merge.")

    counts = merge_matrices(count_matrices)
    expression = merge_matrices(expression_matrices)
    samples = pd.concat(sample_tables, ignore_index=True, sort=False)
    samples = samples.drop_duplicates("import_id").sort_values(["dataset", "sample_id", "import_id"])

    write_table(counts, output_dir / args.counts_name)
    write_table(expression, output_dir / args.expression_name)
    write_table(samples, output_dir / args.sample_table_name)

    print(f"[OK] Counts globais: {output_dir / args.counts_name}")
    print(f"[OK] Expressao global: {output_dir / args.expression_name}")
    print(f"[OK] Amostras globais: {output_dir / args.sample_table_name}")
    if missing:
        print("[WARN] Projetos ignorados sem saidas completas: " + ", ".join(missing))


if __name__ == "__main__":
    main()
