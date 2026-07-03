#!/usr/bin/env bash

export PIPELINE_NAME="example_rnaseq"
export ORGANISM_NAME="Example organism"
export PIPELINE_PROJECTS="PRJNA000000"

export SCRATCH_ROOT="/scratch/my_user/example_rnaseq"
export CONDA_BASE="/path/to/miniconda3"

export TRANSCRIPTS_URL="https://example.org/transcripts.fa.gz"
export GTF_URL="https://example.org/annotation.gtf.gz"

export QUANT_METHOD="salmon"
export FASTQ_LAYOUT="paired"
export PIPELINE_EXECUTOR="slurm"
export PIPELINE_STORAGE_MODE="full"
export PIPELINE_COMPRESS_RESULTS=1

export RUN_SALMON_INDEX=1
export RUN_STAR_GTF_INDEX=0
export RUN_STAR_INDEX=0

export RUN_BATCH_CORRECTION=0
export RUN_DTU_ANALYSIS=0
export RUN_SPLICING_ANALYSIS=0
export RUN_WGCNA_ANALYSIS=0
export RUN_MFUZZ_ANALYSIS=0
export RUN_GENE_REPORT=0
