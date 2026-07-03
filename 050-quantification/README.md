# 050-quantification

Import quantification files and generate gene-level matrices.

This step follows `QUANT_METHOD`:

- `salmon`: imports Salmon `quant.sf` with `tximport` and writes counts plus
  TPM.
- `star`: imports STAR `ReadsPerGene.out.tab` and writes counts plus CPM.

When `--all` is used, the step first writes/updates project-specific matrices
for the projects in `PIPELINE_PROJECTS`, then rebuilds the global matrices from
all complete project-specific outputs already present in `QUANTIFICATION_DIR`.
This preserves earlier projects when new projects are processed later.

## Before You Run

Confirm these exist:

```text
${QUANT_DIR}/<PROJECT>/<sample_id>/quant.sf                  # Salmon mode
${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/ReadsPerGene.out.tab # STAR mode
${REF_GTF}                                                   # Salmon mode
025-parse/030-metadata_final/AllProjects_metadata_new.csv
```

Relevant defaults in `config/pipeline_config.sh`:

- `METADATA_FINAL_NEW`
- `QUANT_METHOD`
- `QUANT_DIR`
- `STAR_QUANT_DIR`
- `STAR_GENECOUNT_COLUMN`
- `REF_GTF`
- `QUANTIFICATION_DIR`
- `QUANT_COUNTS_MATRIX_NAME`
- `SALMON_TPM_MATRIX_NAME`
- `STAR_CPM_MATRIX_NAME`
- `QUANT_SAMPLES_NAME`
- `TX2GENE_NAME`
- `PIPELINE_COMPRESS_RESULTS`

Set these in `config/user_settings.sh` to move or rename quantification
outputs:

```bash
export QUANT_DIR="${PROJECT_DIR}/040-alignment/quants"
export STAR_QUANT_DIR="${PROJECT_DIR}/040-alignment/star_quant"
export QUANTIFICATION_DIR="${PROJECT_DIR}/050-quantification"
export PIPELINE_COMPRESS_RESULTS=1
export QUANT_COUNTS_MATRIX_NAME="counts_matrix.tsv.gz"
export SALMON_TPM_MATRIX_NAME="tpm_matrix.tsv.gz"
export STAR_CPM_MATRIX_NAME="star_cpm_matrix.tsv.gz"
export QUANT_SAMPLES_NAME="quant_samples.tsv.gz"
export TX2GENE_NAME="tx2gene.tsv.gz"
```

Relative paths are resolved from the repository root.

## Run

Normally `rnaseq_pipeline.sh` submits this step. To run manually through Slurm:

```bash
bash scripts/050-quantification/run_quantification_slurm.sh --all
```

Force a method manually, if needed:

```bash
bash scripts/050-quantification/run_quantification_slurm.sh --all --method salmon
bash scripts/050-quantification/run_quantification_slurm.sh --all --method star
```

Without Slurm:

```bash
bash scripts/050-quantification/run_quantification_slurm.sh --all --local
```

Run one project:

```bash
bash scripts/050-quantification/run_quantification_slurm.sh PRJXXXX
```

Preview only:

```bash
bash scripts/050-quantification/run_quantification_slurm.sh --all --sbatch-dry-run
bash scripts/050-quantification/run_quantification_slurm.sh --all --local --sbatch-dry-run
```

## Outputs

- `${QUANT_COUNTS_MATRIX_FILE}`
- `${SALMON_TPM_MATRIX_FILE}` when `QUANT_METHOD=salmon`
- `${STAR_CPM_MATRIX_FILE}` when `QUANT_METHOD=star`
- `${QUANT_SAMPLES_FILE}`
- `${TX2GENE_FILE}` when `QUANT_METHOD=salmon`
- project-specific files such as `<PROJECT>_counts_matrix.tsv.gz`,
  `<PROJECT>_tpm_matrix.tsv.gz`, and `<PROJECT>_quant_samples.tsv.gz`

The sample table includes `quant_method`, `expression_unit`, and, for STAR,
the selected `star_count_column`.
