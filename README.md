# Generic RNA-seq Pipeline for Slurm or Local Runs

This repository contains a modular RNA-seq pipeline designed for reproducible
analysis on Slurm/HPC systems or on a local workstation. It is organism
agnostic: genome FASTA, transcript FASTA, GFF3/GTF annotation, project
accessions, metadata parser YAML files, scratch paths, Conda environments, and
analysis defaults are configured by the user.

The pipeline starts from ENA/SRA FASTQ links and produces parsed metadata,
FastQC/MultiQC reports, trimmed and merged FASTQs, Salmon or STAR gene-count
quantification, gene-level count and TPM/CPM matrices, optional batch-corrected
matrices, DESeq2 outputs, optional transcript-usage/splicing/co-expression
analyses, and an optional candidate gene report.

## Directory Layout

- `config/`: central configuration, example config, and metadata template
- `scripts/`: active executables for every pipeline step, validation utilities, and shared functions
- `slurm/`: Slurm notes
- `envs/`: Conda environment templates
- `examples/`: minimal example skeleton
- `000-logs/`: high-level Slurm logs created by the orchestrator
- `010-reference/`: reference FASTA/GFF/GTF files and Salmon/STAR indexes
- `020-data-download/`: FASTQ download configs and link files
- `025-parse/`: ENA metadata download, project YAML parsers, and merged metadata
- `030-qc-fastq/`: FastQC, trimming, merged FASTQs, and MultiQC
- `040-alignment/`: Salmon or STAR plans and quantification outputs
- `050-quantification/`: imported count and expression matrices
- `055-batch-correction/`: optional batch assessment/correction
- `060-deg-analysis/`: DESeq2 plans, contrasts, plots, and summaries
- `070-dtu-analysis/`: optional differential transcript usage screening
- `080-splicing/`: optional alternative splicing analysis with rMATS
- `085-wgcna/`: optional weighted gene co-expression network analysis
- `086-mfuzz/`: optional soft clustering of expression profiles
- `090-search-gene/`: optional candidate gene/group report

## Requirements

Create the required Conda environments from `envs/` or provide equivalent tools
on `PATH`.

Core tools for both modes:

- Python 3 with pandas
- R with tximport, DESeq2, rtracklayer, ggplot2, and related packages
- FastQC and MultiQC
- Trim Galore
- Salmon for the default transcript quantification mode
- STAR and gffread for STAR genome-alignment quantification or STAR indexes
- rMATS for optional alternative splicing
- WGCNA and Mfuzz for optional expression-pattern modules
- wget

Additional tools for Slurm mode:

- `sbatch` and `squeue`

## Configure

For most users, start with the setup helper:

```bash
python scripts/bootstrap_project.py \
  --project PRJNA000000 \
  --organism "Example organism" \
  --scratch-root /scratch/my_user/rnaseq_project \
  --conda-base /path/to/miniconda3 \
  --quant-method salmon
```

This creates `config/user_settings.sh`, the per-project download config, the
FASTQ link placeholder, and the metadata parser YAML. Then edit the generated
reference paths, FASTQ URLs, and parser regular expressions.
Add `--dry-run` to preview the files without writing them.

Manual setup is still possible:

```bash
cp config/user_settings_template.sh config/user_settings.sh
nano config/user_settings.sh
```

Required user-specific inputs:

- `ORGANISM_NAME`
- `PIPELINE_PROJECTS`
- `SCRATCH_ROOT`
- `CONDA_BASE`
- reference URLs or local reference file paths

Choose where jobs run:

```bash
export PIPELINE_EXECUTOR="slurm"   # default
export PIPELINE_EXECUTOR="local"   # no sbatch/squeue
```

Choose the quantification method:

```bash
export QUANT_METHOD="salmon"  # default; outputs TPM
export QUANT_METHOD="star"    # genome alignment; outputs STAR counts and CPM
```

Use Salmon when you have a transcript FASTA. Use STAR when you want
genome-aligned gene counts through `ReadsPerGene.out.tab`. By default STAR
does not write BAM files, which keeps disk use much smaller. Set
`STAR_WRITE_BAM=1` only when you need sorted BAMs, for example for rMATS,
splicing analysis, or manual genome-browser inspection.

Optionally customize quantification outputs in `config/user_settings.sh`:

```bash
export QUANT_DIR="${PROJECT_DIR}/040-alignment/quants"
export STAR_QUANT_DIR="${PROJECT_DIR}/040-alignment/star_quant"
export QUANTIFICATION_DIR="${PROJECT_DIR}/050-quantification"
export QUANT_COUNTS_MATRIX_NAME="counts_matrix.tsv"
export SALMON_TPM_MATRIX_NAME="tpm_matrix.tsv"
export STAR_CPM_MATRIX_NAME="star_cpm_matrix.tsv"
export QUANT_SAMPLES_NAME="quant_samples.tsv"
```

Optional STAR BAM output:

```bash
export STAR_WRITE_BAM=0               # default: GeneCounts only, no BAM
export STAR_WRITE_BAM=1               # sorted BAMs for splicing/rMATS/inspection
export STAR_LIMIT_BAM_SORT_RAM=24000000000
```

Choose how much generated data to keep:

```bash
export PIPELINE_STORAGE_MODE="full"      # keep everything; best for debugging/reruns
export PIPELINE_STORAGE_MODE="balanced"  # remove per-run trimmed FASTQs and individual FastQC dirs after Salmon/STAR
export PIPELINE_STORAGE_MODE="minimal"   # also remove raw FASTQs, merged trimmed FASTQs, STAR BAMs, and STAR temp dirs
```

Use `balanced` for most production runs when disk is limited. It keeps the
MultiQC report, raw FASTQs, merged trimmed FASTQs, and quantification outputs.
Use `minimal` only when you accept that rerunning QC/alignment will require
downloading and processing FASTQs again.

Enable optional downstream analyses only when the required inputs exist:

```bash
export RUN_DTU_ANALYSIS=1       # Salmon quant.sf + tx2gene.tsv(.gz)
export RUN_SPLICING_ANALYSIS=1  # STAR BAMs + GTF + rMATS
export RUN_WGCNA_ANALYSIS=1     # expression matrix + metadata traits
export RUN_MFUZZ_ANALYSIS=1     # expression matrix + ordered/profile metadata
```

DTU works from Salmon transcript-level quantification. Splicing works from STAR
`Aligned.sortedByCoord.out.bam` files and needs `QUANT_METHOD=star` in the same
run or pre-existing STAR BAMs under `STAR_QUANT_DIR`. WGCNA and Mfuzz work from
the imported expression matrix in step 050; Mfuzz additionally needs a metadata
column such as `stage` to order or group profiles.

`config/pipeline_config.sh` contains advanced defaults, directory variables,
and helper functions. Most users should leave it alone after creating
`config/user_settings.sh`.

Run commands from the repository root. The numbered directories are for inputs,
logs, work files, and results. Active scripts are centralized under `scripts/`.

For each project, create:

```text
020-data-download/datasets/<PROJECT>/config.yaml
025-parse/020-metadata_parsers/<PROJECT>/configs/<PROJECT>.yaml
```

## Metadata

The final parsed metadata must contain:

```text
dataset sample_id run_accession
```

Recommended columns for downstream analyses:

```text
condition stage tissue sex batch replicate
```

Use `config/metadata_template.tsv` as a minimal contract reference. Project
YAML files in `025-parse/020-metadata_parsers/<PROJECT>/configs/` are
responsible for creating these columns.

Validate an existing metadata table:

```bash
python scripts/validate_metadata.py \
  --metadata 025-parse/030-metadata_final/AllProjects_metadata_new.csv
```

Use `--strict` before production submissions to turn recommended-field
warnings into errors.

## Run

Before submitting jobs, audit the setup:

```bash
bash scripts/validate_config.sh config/pipeline_config.sh
```

Complete run with Slurm dependencies:

```bash
bash rnaseq_pipeline.sh --all
```

Complete run locally, without Slurm:

```bash
bash rnaseq_pipeline.sh --all --local
```

Inspect the job graph without submitting:

```bash
bash rnaseq_pipeline.sh --all --dry-run
bash rnaseq_pipeline.sh --all --local --dry-run
```

Run one or more coarse steps:

```bash
bash rnaseq_pipeline.sh --step reference
bash rnaseq_pipeline.sh --step metadata --step qc
bash rnaseq_pipeline.sh --step salmon --step tximport --step deg
bash rnaseq_pipeline.sh --step star --step tximport --step splicing --step report
bash rnaseq_pipeline.sh --step tximport --step wgcna --step mfuzz
```

Supported steps:

```text
reference download metadata qc salmon star tximport batch deg dtu splicing wgcna mfuzz report
```

`main.sh` is kept as a compatibility wrapper around `rnaseq_pipeline.sh`; new
runs should use `rnaseq_pipeline.sh` directly.

## Workflow

Full execution order:

```text
reference
download + metadata -> qc
qc + reference -> salmon
qc + reference -> star (when QUANT_METHOD=star)
salmon/star -> quantification import
quantification import -> batch (optional)
quantification import or batch -> deg
quantification import -> dtu (optional; Salmon)
star alignment with STAR_WRITE_BAM=1 -> splicing (optional; STAR BAMs)
quantification import -> wgcna (optional)
quantification import -> mfuzz (optional)
deg + dtu + splicing + wgcna + mfuzz -> report (optional)
```

In Slurm mode, sample/project-level steps are submitted independently where
possible and downstream steps use `--dependency=afterok`. In local mode, the
same steps run in dependency order on the current machine; Slurm arrays are
simulated sequentially with `SLURM_ARRAY_TASK_ID`.

## Outputs

Key files and directories:

- `010-reference/salmon_index/`
- `${SCRATCH_ROOT}/<PROJECT>/fastq_ftp/`
- `025-parse/030-metadata_final/AllProjects_metadata_new.csv`
- `${SCRATCH_ROOT}/<PROJECT>/multiqc_030/`
- `${QUANT_DIR}/<PROJECT>/<sample_id>/quant.sf`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/ReadsPerGene.out.tab`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/Aligned.sortedByCoord.out.bam` only when `STAR_WRITE_BAM=1`
- `${QUANT_COUNTS_MATRIX_FILE}`
- `${SALMON_TPM_MATRIX_FILE}`
- `${STAR_CPM_MATRIX_FILE}`
- `${QUANT_SAMPLES_FILE}`
- `055-batch-correction/all_projects/counts_batch_corrected.tsv`
- `060-deg-analysis/all_projects/raw/DEGs_all_results.tsv`
- `070-dtu-analysis/all_projects/dtu_results.tsv.gz`
- `080-splicing/all_projects/<variable>/<level_a>_vs_<level_b>/splicing_summary.tsv.gz`
- `085-wgcna/all_projects/module_assignments.tsv.gz`
- `085-wgcna/all_projects/hub_genes.tsv.gz`
- `086-mfuzz/all_projects/mfuzz_membership.tsv.gz`
- `090-search-gene/results/gene_set_report.html`
- `090-search-gene/results/tables/stage_specificity_summary.tsv.gz`

When `PIPELINE_STORAGE_MODE` is `balanced` or `minimal`, some scratch
intermediates are removed automatically after step 040 succeeds for each
project. The final matrices, DEG outputs, MultiQC summary, and report outputs
are kept.

## Recovering Failed Jobs

1. Inspect `000-logs/` and the step-specific `logs/` directory.
2. Fix the cause, such as a missing FASTQ, bad YAML parser, missing reference,
   missing Conda environment, or insufficient Slurm memory/time.
3. Rerun the failed coarse step:

```bash
bash rnaseq_pipeline.sh --step <step>
```

## Notes

- Use Slurm mode on the remote HPC server for production-sized datasets.
- Use local mode for small datasets, smoke tests, teaching, and machines
  without `sbatch`.
- Large generated data should not be committed.
- Bundled project-specific folders are examples or legacy inputs. Only
  `PIPELINE_PROJECTS` controls what is submitted.
