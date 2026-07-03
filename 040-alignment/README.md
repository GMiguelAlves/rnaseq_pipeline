# 040-alignment

Create a sample-level plan from the step 030 QC plan and submit one job per
biological sample.

This step has two modes:

- `QUANT_METHOD=salmon`: runs Salmon against the transcriptome and writes
  `quant.sf`.
- `QUANT_METHOD=star`: runs STAR against the annotated genome index, writes
  `ReadsPerGene.out.tab` for gene counts. By default it does not write BAM
  files; set `STAR_WRITE_BAM=1` only when sorted BAMs are needed.

## Before You Run

Confirm these exist:

```text
030-qc-fastq/work/<PROJECT>_qc_plan.csv
010-reference/salmon_index/       # Salmon mode
010-reference/star_index_gtf/      # STAR mode
```

Relevant defaults in `config/pipeline_config.sh`:

- `SALMON_INDEX_DIR`
- `QUANT_DIR`
- `SALMON_CONCURRENCY`
- `STAR_QUANT_INDEX_DIR`
- `STAR_QUANT_DIR`
- `STAR_QUANT_CONCURRENCY`
- `STAR_GENECOUNT_COLUMN`
- `STAR_READ_FILES_COMMAND`
- `STAR_WRITE_BAM`
- `STAR_LIMIT_BAM_SORT_RAM`
- `STAR_EXTRA_ARGS`
- `PIPELINE_STORAGE_MODE`

## Run

Normally `rnaseq_pipeline.sh` submits this step and chooses the script from
`QUANT_METHOD`.

Run Salmon manually:

```bash
bash scripts/040-alignment/run_alignment_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv"
```

Run STAR manually:

```bash
bash scripts/040-alignment/run_star_quant_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv"
```

Without Slurm, add `--local`:

```bash
bash scripts/040-alignment/run_alignment_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv" \
  --local

bash scripts/040-alignment/run_star_quant_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv" \
  --local
```

Preview only:

```bash
bash scripts/040-alignment/run_alignment_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv" \
  --dry-run

bash scripts/040-alignment/run_star_quant_project.sh \
  PRJXXXX \
  "$PWD/030-qc-fastq/work/PRJXXXX_qc_plan.csv" \
  --dry-run
```

## Outputs

- `040-alignment/work/<PROJECT>_salmon_plan.csv`
- `${QUANT_DIR}/<PROJECT>/<sample_id>/quant.sf`
- `040-alignment/work/<PROJECT>_star_plan.csv`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/ReadsPerGene.out.tab`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/Aligned.sortedByCoord.out.bam` only when `STAR_WRITE_BAM=1`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/Log.final.out`

For ordinary counts/CPM/DEG/report workflows, keep:

```bash
export STAR_WRITE_BAM=0
```

If you need sorted BAMs for rMATS/splicing or manual inspection:

```bash
export STAR_WRITE_BAM=1
export STAR_LIMIT_BAM_SORT_RAM=24000000000
```

When STAR reports `not enough memory for BAM sorting`, set
`STAR_LIMIT_BAM_SORT_RAM` to at least the number printed in the error message.

## Disk Cleanup

This step can clean large scratch intermediates automatically after the Salmon
or STAR job array finishes successfully.

Set in `config/user_settings.sh`:

```bash
export PIPELINE_STORAGE_MODE="full"      # keep everything
export PIPELINE_STORAGE_MODE="balanced"  # remove fastqc_raw/, fastqc_trimmed_runs/, fastqc_merged/, trimmed_runs/
export PIPELINE_STORAGE_MODE="minimal"   # also remove fastq_ftp/, trimmed_merged/, STAR BAMs, and STAR _STARtmp dirs
```

`multiqc_030/`, quantification outputs, imported matrices, DEG outputs, and
gene reports are kept. Use `minimal` only when you do not need to rerun
QC/alignment from existing FASTQs.

Preview cleanup for one project:

```bash
bash scripts/040-alignment/cleanup_project_storage.sh PRJXXXX --storage-mode balanced --dry-run
```

Run cleanup manually for an already finished project:

```bash
bash scripts/040-alignment/cleanup_project_storage.sh PRJXXXX --storage-mode balanced
```

The old `salmon_quant.sh` lives under `040-alignment/legacy/` and is ignored.
Use scripts under `scripts/040-alignment/` for new manual runs.
