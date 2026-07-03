# 086 - Mfuzz

Optional soft clustering of expression profiles with Mfuzz.

Enable it in `config/user_settings.sh`:

```bash
export RUN_MFUZZ_ANALYSIS=1
export MFUZZ_TIME_VARIABLE="stage"
```

Inputs:

- `${EXPRESSION_MATRIX_FILE}` and project-specific expression matrices from step 050
- `${QUANT_SAMPLES_FILE}` and project-specific sample tables
- parsed metadata containing the selected `MFUZZ_TIME_VARIABLE`

Outputs are written per project and, with `--include-all`, for all projects:

```text
086-mfuzz/<PROJECT>/mfuzz_summary.tsv.gz
086-mfuzz/<PROJECT>/mfuzz_membership.tsv.gz
086-mfuzz/<PROJECT>/cluster_centers.tsv.gz
086-mfuzz/<PROJECT>/cluster_summary.tsv.gz
086-mfuzz/<PROJECT>/expression_context_matrix.tsv.gz
086-mfuzz/<PROJECT>/plots/cluster_centers.png
```

Manual run:

```bash
bash scripts/086-mfuzz/run_mfuzz_analysis_slurm.sh --include-all
```

Useful settings:

```bash
export MFUZZ_TIME_VARIABLE="stage"
export MFUZZ_TIME_LEVELS="egg,miracidium,sporocyst,cercaria,schistosomulum,adult"
export MFUZZ_GROUP_COLUMNS=""   # e.g. condition,tissue if you want context profiles
export MFUZZ_CLUSTERS=6
export MFUZZ_M=0                # 0 estimates fuzzifier automatically
```

Scopes with too few samples, genes, or time/context levels are marked as
`skipped` in `mfuzz_summary.tsv.gz`.
