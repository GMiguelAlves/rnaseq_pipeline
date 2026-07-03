# Generic Organism Starter Example

This folder shows the minimal files a new organism/project needs before the
pipeline can run. It intentionally avoids organism-specific gene IDs, stages,
and project accessions.

Use the bootstrap helper from the repository root:

```bash
python scripts/bootstrap_project.py \
  --project PRJNA000000 \
  --organism "Example organism" \
  --scratch-root /scratch/my_user/example_rnaseq \
  --conda-base /path/to/miniconda3 \
  --quant-method salmon
```

Add `--dry-run` to preview the generated file paths without writing anything.

Then edit the generated files:

```text
config/user_settings.sh
020-data-download/datasets/PRJNA000000/config.yaml
020-data-download/datasets/PRJNA000000/ena-file-download-read_run-PRJNA000000-fastq_ftp.sh
025-parse/020-metadata_parsers/PRJNA000000/configs/PRJNA000000.yaml
```

The parser YAML is the only project-specific biological mapping. Change the
`regex_map` values to match the sample names or ENA metadata fields for your
organism.

Before submitting jobs:

```bash
bash scripts/validate_config.sh config/pipeline_config.sh
bash rnaseq_pipeline.sh --all --dry-run
```
