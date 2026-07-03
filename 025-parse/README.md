# 025-parse

Download, validate, parse, and merge metadata with `metaQC`.

## Before You Run

For each project listed in `PIPELINE_PROJECTS`, create a parser YAML:

```text
025-parse/020-metadata_parsers/<PROJECT>/configs/<PROJECT>.yaml
```

Recommended starter command from the repository root:

```bash
python scripts/bootstrap_project.py \
  --project PRJNA000000 \
  --organism "Example organism" \
  --scratch-root /scratch/my_user/rnaseq_project \
  --conda-base /path/to/miniconda3
```

This creates a generic YAML with editable `regex_map` sections for
`condition`, `stage`, `tissue`, `sex`, and `replicate`.

Optional files:

```text
025-parse/020-metadata_parsers/<PROJECT>/configs/<PROJECT>_enrich.yaml
025-parse/020-metadata_parsers/<PROJECT>/author_metadata.tsv
```

Start from:

```text
025-parse/020-metadata_parsers/TEMPLATE_project.yaml
```

The final metadata must contain:

- `dataset`
- `sample_id`
- `run_accession`

Add biological columns such as `condition`, `stage`, `tissue`, `sex`, `batch`,
and covariates in the YAML.

Only projects listed in `PIPELINE_PROJECTS` are used by `rnaseq_pipeline.sh`. Any bundled
project parser folders are examples or legacy inputs, not pipeline defaults.

## Run

Normally `rnaseq_pipeline.sh` submits this step. To run one project manually
from the repository root:

```bash
sbatch --chdir=025-parse scripts/025-parse/run_metaqc.sh PRJXXXX
sbatch --chdir=025-parse scripts/025-parse/run_parse.sh PRJXXXX
```

Without Slurm:

```bash
bash scripts/025-parse/run_metaqc.sh PRJXXXX
bash scripts/025-parse/run_parse.sh PRJXXXX
```

After all projects have been parsed:

```bash
sbatch --chdir=025-parse scripts/025-parse/run_merge.sh
```

Without Slurm:

```bash
bash scripts/025-parse/run_merge.sh
```

Then validate the merged metadata:

```bash
python scripts/validate_metadata.py \
  --metadata 025-parse/030-metadata_final/AllProjects_metadata_new.csv \
  --strict
```

## Outputs

- `025-parse/010-raw_metadata/<PROJECT>.tsv`
- `025-parse/015-intermediate_folder/<PROJECT>_base.csv`
- `025-parse/015-intermediate_folder/<PROJECT>_enriched.csv` when enrich is used
- `025-parse/020-metadata_parsers/Allprojects/<PROJECT>_parsed.csv`
- `025-parse/030-metadata_final/AllProjects_metadata.csv`
