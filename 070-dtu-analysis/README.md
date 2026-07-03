# 070 - DTU Analysis

Optional differential transcript usage screening from Salmon `quant.sf` files.

Enable it in `config/user_settings.sh`:

```bash
export RUN_DTU_ANALYSIS=1
```

Inputs:

- `025-parse/030-metadata_final/AllProjects_metadata_new.csv`
- `${QUANT_DIR}/<PROJECT>/<sample_id>/quant.sf`
- `050-quantification/tx2gene.tsv` or `tx2gene.tsv.gz`

Outputs are written per project and, with `--include-all`, for all projects:

```text
070-dtu-analysis/<PROJECT>/dtu_results.tsv.gz
070-dtu-analysis/<PROJECT>/dtu_significant.tsv.gz
070-dtu-analysis/<PROJECT>/transcript_usage_long.tsv.gz
070-dtu-analysis/<PROJECT>/dtu_summary.tsv.gz
070-dtu-analysis/all_projects/dtu_results.tsv.gz
```

Manual run:

```bash
bash scripts/070-dtu-analysis/run_dtu_analysis_slurm.sh --include-all
```

The current implementation is a robust transcript-usage screen using a
Kruskal-Wallis test across metadata groups. It is intentionally lightweight and
keeps the pipeline stable; the output tables can later feed a stricter DRIMSeq
or IsoformSwitchAnalyzeR workflow if needed.
