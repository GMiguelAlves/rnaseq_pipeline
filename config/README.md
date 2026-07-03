# Configuration

Most users should edit only:

```text
config/user_settings.sh
```

The easiest way to create it and the per-project starter files is:

```bash
python scripts/bootstrap_project.py \
  --project PRJNA000000 \
  --organism "Example organism" \
  --scratch-root /scratch/my_user/rnaseq_project \
  --conda-base /path/to/miniconda3
```

This writes `config/user_settings.sh`, the dataset download config, the FASTQ
link placeholder, and the metadata parser YAML for each project. You can pass
`--project` more than once. Add `--dry-run` to preview the files without
writing them.

Manual setup is also supported:

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
- `PIPELINE_COMPRESS_RESULTS`: keep as `1` to write large tabular outputs as
  `.tsv.gz`; downstream steps read `.tsv` and `.tsv.gz`
- `RUN_DTU_ANALYSIS`, `RUN_SPLICING_ANALYSIS`, `RUN_WGCNA_ANALYSIS`, and
  `RUN_MFUZZ_ANALYSIS`: optional downstream modules, disabled by default

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
export STAR_WRITE_BAM=0                    # default: GeneCounts only, no BAM
export STAR_LIMIT_BAM_SORT_RAM=24000000000 # used only when STAR_WRITE_BAM=1
```

For ordinary STAR quantification, keep `STAR_WRITE_BAM=0`. This writes
`ReadsPerGene.out.tab` and logs without creating large sorted BAM files. Use
`STAR_WRITE_BAM=1` only for rMATS/splicing or manual inspection; if STAR then
reports `not enough memory for BAM sorting`, increase `STAR_LIMIT_BAM_SORT_RAM`
to the value requested in the STAR error message.

Quantification output options:

```bash
# Per-sample quantification outputs
export QUANT_DIR="${PROJECT_DIR}/040-alignment/quants"
export STAR_QUANT_DIR="${PROJECT_DIR}/040-alignment/star_quant"

# Imported matrices/tables
export QUANTIFICATION_DIR="${PROJECT_DIR}/050-quantification"
export QUANT_COUNTS_MATRIX_NAME="counts_matrix.tsv.gz"
export SALMON_TPM_MATRIX_NAME="tpm_matrix.tsv.gz"
export STAR_CPM_MATRIX_NAME="star_cpm_matrix.tsv.gz"
export QUANT_SAMPLES_NAME="quant_samples.tsv.gz"
export TX2GENE_NAME="tx2gene.tsv.gz"
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
export PIPELINE_COMPRESS_RESULTS=1       # gzip TSV-like result tables
```

`balanced` keeps `multiqc_030/`, raw downloads, merged trimmed FASTQs, and
quantification outputs. `minimal` saves more disk, but rerunning QC/alignment
will require downloading and processing FASTQs again.

For a project that already finished quantification, preview cleanup manually:

```bash
bash scripts/040-alignment/cleanup_project_storage.sh PRJXXXX --storage-mode balanced --dry-run
```

Optional DTU/splicing modules:

```bash
# Transcript usage / DTU screening from Salmon quant.sf files
export RUN_DTU_ANALYSIS=1
export DTU_TEST_VARIABLES="condition,stage,sex,tissue,infection_mode"

# Alternative splicing with rMATS from STAR sorted BAM files
export RUN_SPLICING_ANALYSIS=1
export QUANT_METHOD="star"
export STAR_WRITE_BAM=1
export SPLICING_READ_LENGTH=100
export SPLICING_LIB_TYPE="fr-unstranded"  # fr-firststrand / fr-secondstrand if needed
```

DTU requires Salmon transcript-level quantifications plus `tx2gene.tsv(.gz)`
from step 050. Splicing requires STAR `Aligned.sortedByCoord.out.bam` files and
a GTF annotation. When splicing is enabled and `STAR_WRITE_BAM` is not set, the
pipeline enables BAM output by default; `minimal` storage mode keeps STAR BAMs
because rMATS needs them.

Optional WGCNA/Mfuzz modules:

```bash
# Co-expression network modules from the imported expression matrix
export RUN_WGCNA_ANALYSIS=1
export WGCNA_TRAIT_COLUMNS="condition,stage,sex,tissue,batch,dataset"
export WGCNA_MIN_SAMPLES=12

# Soft clustering of temporal/profile expression patterns
export RUN_MFUZZ_ANALYSIS=1
export MFUZZ_TIME_VARIABLE="stage"
export MFUZZ_TIME_LEVELS=""
export MFUZZ_CLUSTERS=6
```

WGCNA and Mfuzz use expression matrices from step 050. WGCNA correlates module
eigengenes against configured metadata traits. Mfuzz needs the selected
`MFUZZ_TIME_VARIABLE` in metadata; set `MFUZZ_TIME_LEVELS` when you want an
explicit biological order rather than the order found in the sample table.

Do not replace `config/pipeline_config.sh`. It contains advanced defaults and
helper functions used by the step scripts. It also defines the active script
directories under `scripts/`; most users do not need to edit those paths.

Run this before submitting jobs:

```bash
bash scripts/validate_config.sh config/pipeline_config.sh
```

After metadata parsing finishes, audit the biological columns:

```bash
python scripts/validate_metadata.py \
  --metadata 025-parse/030-metadata_final/AllProjects_metadata_new.csv \
  --strict
```

Run the pipeline from the repository root:

```bash
bash rnaseq_pipeline.sh --all --dry-run
bash rnaseq_pipeline.sh --all
bash rnaseq_pipeline.sh --all --local
```
