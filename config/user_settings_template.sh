#!/usr/bin/env bash

# Copy this file to config/user_settings.sh and edit only this small block.
# Leave config/pipeline_config.sh alone unless you need advanced behavior.

# 1) Name your analysis.
export PIPELINE_NAME="rnaseq_project"
export ORGANISM_NAME="My organism"

# 2) List the ENA/SRA projects to process.
# Use spaces or commas.
export PIPELINE_PROJECTS="PRJXXXX PRJYYYY"

# 3) Set the scratch/work directory on the Slurm server.
export SCRATCH_ROOT="/scratch/my_user/rnaseq_project"

# 4) Point to Conda on the Slurm server.
# Absolute paths are safest. Relative paths are resolved from this file's
# directory, so "../miniconda3" means "<project>/miniconda3".
export CONDA_BASE="/path/to/miniconda3"

# 5) Provide reference inputs.
# Option A: URLs. The pipeline downloads and decompresses them.
export GENOME_URL="https://example.org/genome.fa.gz"
export TRANSCRIPTS_URL="https://example.org/transcripts.fa.gz"
export GFF3_URL="https://example.org/annotation.gff3.gz"
# If your organism provides GTF instead of GFF3, use this line instead.
# export GTF_URL="https://example.org/annotation.gtf.gz"

# Option B: local files already present on the server.
# If using local files, uncomment and edit these instead of the URLs above.
# export REF_GENOME_FA="/path/to/genome.fa"
# export REF_TRANSCRIPTS_FA="/path/to/transcripts.fa"
# export REF_GFF3="/path/to/annotation.gff3"
# export REF_GTF="/path/to/annotation.gtf"

# 6) Choose the quantification method.
# salmon: transcript-level quantification, imported with tximport.
# star: genome alignment plus STAR GeneCounts. Needs genome + GTF/GFF3.
export QUANT_METHOD="salmon"       # Use "salmon" or "star".
export STAR_GENECOUNT_COLUMN="unstranded"

# Optional: customize where quantification outputs are written.
# Relative paths are resolved from the repository root.
#
# Salmon per-sample quant.sf files:
# export QUANT_DIR="${PROJECT_DIR}/040-alignment/quants"
#
# STAR per-sample ReadsPerGene/BAM/log files:
# export STAR_QUANT_DIR="${PROJECT_DIR}/040-alignment/star_quant"
#
# Imported matrices and sample tables:
# export QUANTIFICATION_DIR="${PROJECT_DIR}/050-quantification"
#
# Optional: customize all-project matrix/table file names written by step 050.
# Project-specific manual imports still use <PROJECT>_counts_matrix.tsv and
# <PROJECT>_quant_samples.tsv unless you pass names on the command line.
# export QUANT_COUNTS_MATRIX_NAME="counts_matrix.tsv"
# export SALMON_TPM_MATRIX_NAME="tpm_matrix.tsv"
# export STAR_CPM_MATRIX_NAME="star_cpm_matrix.tsv"
# export QUANT_SAMPLES_NAME="quant_samples.tsv"
# export TX2GENE_NAME="tx2gene.tsv"

# 7) Choose where jobs run.
export PIPELINE_EXECUTOR="slurm"   # Use "local" to run without sbatch/squeue.
export LOCAL_CPUS_PER_TASK=8       # Used only when PIPELINE_EXECUTOR="local".

# 8) Choose how much intermediate data to keep.
# full: keep everything. Best for debugging and reruns.
# balanced: after Salmon/STAR succeeds, remove individual FastQC folders and
#           run-level trimmed FASTQs. Keeps raw FASTQs, merged trimmed FASTQs,
#           MultiQC and quantification outputs.
# minimal: after Salmon/STAR succeeds, remove raw FASTQs, all trimmed FASTQs,
#          individual FastQC folders and STAR BAMs. Saves the most disk, but
#          rerunning QC/alignment requires downloading/processing again.
export PIPELINE_STORAGE_MODE="full"  # Use "balanced" or "minimal" to save disk.

# 9) Usually keep these defaults.
if [[ "$QUANT_METHOD" == "star" ]]; then
  export RUN_SALMON_INDEX=0
  export RUN_STAR_GTF_INDEX=1
else
  export RUN_SALMON_INDEX=1
  export RUN_STAR_GTF_INDEX=0
fi
export RUN_STAR_INDEX=0
export RUN_BATCH_CORRECTION=0
export RUN_GENE_REPORT=0

# 10) Conda environment names. Change only if your server uses other names.
export RNA_TOOLS_ENV="rna-tools"
export PYTHON_ENV="python-list"
export R_ANALYSIS_ENV="r-analysis"
export BATCH_CORRECTION_ENV="batch-correction"
