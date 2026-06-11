# Configuration

Most users should edit only:

```text
config/user_settings.sh
```

Create it from the template:

```bash
cp config/user_settings_template.sh config/user_settings.sh
```

Then edit:

- `PIPELINE_NAME`
- `ORGANISM_NAME`
- `PIPELINE_PROJECTS`
- `SCRATCH_ROOT`
- `CONDA_BASE`
- reference URLs or local reference file paths
- `QUANT_METHOD`: use `salmon` or `star`
- quantification output paths, if you want them outside the default numbered
  directories:
  - `QUANT_DIR`
  - `STAR_QUANT_DIR`
  - `QUANTIFICATION_DIR`
- `PIPELINE_EXECUTOR`: use `slurm` on an HPC server or `local` without Slurm
- `LOCAL_CPUS_PER_TASK`: number of local CPU threads used by tools that honor it
- `PIPELINE_STORAGE_MODE`: use `full`, `balanced`, or `minimal` to control
  automatic cleanup of large generated intermediates

`CONDA_BASE` should point to the directory that contains
`etc/profile.d/conda.sh`. Absolute paths are safest. Relative paths are
resolved from the directory containing `config/user_settings.sh`; for example,
`export CONDA_BASE="../miniconda3"` points to `<project>/miniconda3`.

Reference requirements depend on the quantification method:

- `QUANT_METHOD=salmon`: provide transcript FASTA through `TRANSCRIPTS_URL`
  or `REF_TRANSCRIPTS_FA`; the report uses `GFF3_URL`/`GTF_URL` for gene names
  and genomic coordinates when available.
- `QUANT_METHOD=star`: provide genome FASTA plus GTF/GFF3 through
  `GENOME_URL`/`REF_GENOME_FA` and `GTF_URL`/`REF_GTF` or
  `GFF3_URL`/`REF_GFF3`; the pipeline builds/uses `STAR_INDEX_GTF_DIR` and
  imports STAR `ReadsPerGene.out.tab`.

Useful STAR options:

```bash
export QUANT_METHOD="star"
export STAR_GENECOUNT_COLUMN="unstranded"  # or stranded_forward / stranded_reverse
export STAR_QUANT_CONCURRENCY=2
export STAR_READ_FILES_COMMAND="zcat"      # set to "" only for uncompressed FASTQ
```

Quantification output options:

```bash
# Per-sample quantification outputs
export QUANT_DIR="${PROJECT_DIR}/040-alignment/quants"
export STAR_QUANT_DIR="${PROJECT_DIR}/040-alignment/star_quant"

# Imported matrices/tables
export QUANTIFICATION_DIR="${PROJECT_DIR}/050-quantification"
export QUANT_COUNTS_MATRIX_NAME="counts_matrix.tsv"
export SALMON_TPM_MATRIX_NAME="tpm_matrix.tsv"
export STAR_CPM_MATRIX_NAME="star_cpm_matrix.tsv"
export QUANT_SAMPLES_NAME="quant_samples.tsv"
export TX2GENE_NAME="tx2gene.tsv"
```

Relative output paths are interpreted from the repository root. Downstream
steps use the derived files `QUANT_COUNTS_MATRIX_FILE`,
`SALMON_TPM_MATRIX_FILE`, `STAR_CPM_MATRIX_FILE`, `QUANT_SAMPLES_FILE`, and
`EXPRESSION_MATRIX_FILE`.

Storage options:

```bash
export PIPELINE_STORAGE_MODE="full"      # keep everything; safest for reruns
export PIPELINE_STORAGE_MODE="balanced"  # remove per-run trimmed FASTQs and individual FastQC dirs after Salmon/STAR
export PIPELINE_STORAGE_MODE="minimal"   # also remove raw FASTQs, merged trimmed FASTQs, and STAR BAMs
```

`balanced` keeps `multiqc_030/`, raw downloads, merged trimmed FASTQs, and
quantification outputs. `minimal` saves more disk, but rerunning QC/alignment
will require downloading and processing FASTQs again.

For a project that already finished quantification, preview cleanup manually:

```bash
bash scripts/040-alignment/cleanup_project_storage.sh PRJXXXX --storage-mode balanced --dry-run
```

Do not replace `config/pipeline_config.sh`. It contains advanced defaults and
helper functions used by the step scripts. It also defines the active script
directories under `scripts/`; most users do not need to edit those paths.

Run this before submitting jobs:

```bash
bash scripts/validate_config.sh config/pipeline_config.sh
```

Run the pipeline from the repository root:

```bash
bash rnaseq_pipeline.sh --all --dry-run
bash rnaseq_pipeline.sh --all
bash rnaseq_pipeline.sh --all --local
```
