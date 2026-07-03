# Minimal Example Skeleton

This folder mirrors the layout expected by the RNA-seq pipeline without shipping
large FASTQ or reference files.

For a real test, provide:

```text
examples/minimal/fastq/
examples/minimal/reference/
examples/minimal/metadata.tsv
```

Then copy and edit:

```bash
cp config/user_settings_template.sh config/user_settings.sh
```

Set paths in `user_settings.sh` to point to the files in this example or to
real files on your Slurm server. For a local smoke test, set:

```bash
export PIPELINE_EXECUTOR="local"
```

Then run from the repository root:

```bash
bash rnaseq_pipeline.sh --all --local --dry-run
```
