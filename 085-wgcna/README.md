# 085 - WGCNA

Optional weighted gene co-expression network analysis from the imported
expression matrices.

Enable it in `config/user_settings.sh`:

```bash
export RUN_WGCNA_ANALYSIS=1
```

Inputs:

- `${EXPRESSION_MATRIX_FILE}` and project-specific expression matrices from step 050
- `${QUANT_SAMPLES_FILE}` and project-specific sample tables
- parsed metadata with biological traits such as `condition`, `stage`, `sex`,
  `tissue`, `batch`, or `dataset`

Outputs are written per project and, with `--include-all`, for all projects:

```text
085-wgcna/<PROJECT>/wgcna_summary.tsv.gz
085-wgcna/<PROJECT>/soft_threshold.tsv.gz
085-wgcna/<PROJECT>/module_assignments.tsv.gz
085-wgcna/<PROJECT>/module_eigengenes.tsv.gz
085-wgcna/<PROJECT>/module_trait_correlations.tsv.gz
085-wgcna/<PROJECT>/hub_genes.tsv.gz
085-wgcna/<PROJECT>/plots/
```

Manual run:

```bash
bash scripts/085-wgcna/run_wgcna_analysis_slurm.sh --include-all
```

Useful settings:

```bash
export WGCNA_TRAIT_COLUMNS="condition,stage,sex,tissue,batch,dataset"
export WGCNA_MIN_SAMPLES=12
export WGCNA_MIN_GENES=500
export WGCNA_POWER=0          # 0 selects power automatically
```

Scopes with too few samples or genes are marked as `skipped` in
`wgcna_summary.tsv.gz` instead of stopping the whole pipeline.
