# 080 - Alternative Splicing

Optional alternative splicing analysis with rMATS from STAR sorted BAM files.

Enable it in `config/user_settings.sh`:

```bash
export RUN_SPLICING_ANALYSIS=1
export QUANT_METHOD="star"
export SPLICING_READ_LENGTH=100
export SPLICING_LIB_TYPE="fr-unstranded"
```

Inputs:

- `025-parse/030-metadata_final/AllProjects_metadata_new.csv`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/Aligned.sortedByCoord.out.bam`
- `REF_GTF` or `GTF_URL`

Outputs are written per comparison:

```text
080-splicing/<PROJECT>/<variable>/<level_a>_vs_<level_b>/
080-splicing/<PROJECT>/<variable>/<level_a>_vs_<level_b>/splicing_summary.tsv.gz
080-splicing/<PROJECT>/<variable>/<level_a>_vs_<level_b>/significant_events.tsv.gz
080-splicing/work/rmats_plan.csv
```

Manual run:

```bash
bash scripts/080-splicing/run_splicing_analysis_slurm.sh --include-all
```

The plan generator creates comparisons only when each group has at least
`SPLICING_MIN_REPLICATES` samples with available BAM files. If storage mode is
`minimal`, the central config keeps STAR BAMs when splicing is enabled because
rMATS cannot run without them.
