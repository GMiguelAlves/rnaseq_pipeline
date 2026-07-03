# 090-search-gene

Optional exploratory report for genes or gene groups of interest. It combines
the configured expression matrix, sample metadata, GFF/GTF annotation, and DEG
results. When available, it also includes candidate-gene hits from DTU,
alternative splicing, WGCNA, and Mfuzz modules.

DTU hits are read from every `070-dtu-analysis/**/dtu_significant.tsv(.gz)`
file. The overview reports the global DTU evidence found in step 070, while
the individual gene sections show only DTU hits for genes listed in
`genes.txt`.

WGCNA hits are read from `085-wgcna/**/module_assignments.tsv(.gz)` and
`085-wgcna/**/hub_genes.tsv(.gz)`. Mfuzz hits are read from
`086-mfuzz/**/mfuzz_membership.tsv(.gz)`.

The expression matrix is selected by `EXPRESSION_MATRIX_FILE`:

- Salmon mode uses `SALMON_TPM_MATRIX_FILE` and reports TPM.
- STAR mode uses `STAR_CPM_MATRIX_FILE` and reports CPM.

When a GFF3 or GTF annotation is available, the report tables include gene
name, biotype, description, chromosome, start, end, strand, and genomic
location.

The final HTML includes search plus filters for project, section type, groups,
genes, figures, and tables. It reads both `.tsv` and `.tsv.gz` outputs from
upstream steps.

Group and single-gene sections include global `all_projects` figures and, when
multiple datasets are present, project-filtered figures generated from only
that dataset's samples. For example, the `PRJEB14695` group heatmap is written
with a project-specific suffix and does not include samples from `PRJEB32839`.

## Before You Run

Confirm these exist:

```text
${EXPRESSION_MATRIX_FILE}
${QUANT_SAMPLES_FILE}
060-deg-analysis/
070-dtu-analysis/      # optional
080-splicing/          # optional
085-wgcna/             # optional
086-mfuzz/             # optional
010-reference/data/<annotation.gff3-or.gtf>
```

Create a gene list:

```text
090-search-gene/genes.txt
```

Format:

```text
Group A: GENE0001, GENE0002
Group B: annotated_gene_name, another_gene
```

Gene entries can be IDs present in the expression matrix or names available in
the GFF/GTF annotation.

## Optional Organism-Specific Ordering

The report is generic by default. To impose an organism-specific stage order,
set these in `config/pipeline_config.sh`:

```bash
export LIFE_STAGE_LEVELS="stage1,stage2,stage3,unknown"
export STAGE_SYNONYM_MAP="regex1=stage1,regex2=stage2"
```

When `ORGANISM_NAME` contains `Schistosoma mansoni`, the pipeline uses this
life-cycle order by default:

```text
eggs, miracidium, sporocyst_1d, sporocyst_5d, sporocyst_32d, cercariae,
schistosomula_2d, adult_26d, adult, unknown
```

The gene list remains externally supplied through `genes.txt`. The report then
adds a stage-specificity summary for those genes using three evidence types:

- expression concentrated in one/few stages, measured by tau
- significant DESeq2 contrasts where the tested variable is `stage`
- strong Mfuzz cluster membership when temporal clustering was run by `stage`

Relevant thresholds:

```bash
export GENE_REPORT_STAGE_TAU_THRESHOLD=0.60
export GENE_REPORT_STAGE_MIN_EXPRESSION=1
export GENE_REPORT_MFUZZ_MEMBERSHIP_THRESHOLD=0.70
```

Set `ORGANISM_SPECIFIC_REPORTS=1` only when the optional hard-coded panel is
meaningful for the organism and metadata.

## Run

Normally `rnaseq_pipeline.sh` submits this step when `RUN_GENE_REPORT=1`. To run manually:

```bash
bash scripts/090-search-gene/run_gene_report_slurm.sh \
  --genes "$PWD/090-search-gene/genes.txt" \
  --title "Candidate genes"
```

The wrapper uses `EXPRESSION_MATRIX_FILE`, `EXPRESSION_UNIT`, and
`GENE_REPORT_ANNOTATION_FILE` from `config/pipeline_config.sh`. Override them
only when running a custom report:

```bash
bash scripts/090-search-gene/run_gene_report_slurm.sh \
  --genes "$PWD/090-search-gene/genes.txt" \
  --tpm "$PWD/050-quantification/star_cpm_matrix.tsv" \
  --expression-unit CPM \
  --gff "$PWD/010-reference/data/annotation.gtf"
```

Without Slurm:

```bash
bash scripts/090-search-gene/run_gene_report_slurm.sh \
  --genes "$PWD/090-search-gene/genes.txt" \
  --title "Candidate genes" \
  --local
```

Preview only:

```bash
bash scripts/090-search-gene/run_gene_report_slurm.sh \
  --genes "$PWD/090-search-gene/genes.txt" \
  --sbatch-dry-run
bash scripts/090-search-gene/run_gene_report_slurm.sh \
  --genes "$PWD/090-search-gene/genes.txt" \
  --local \
  --sbatch-dry-run
```

## Outputs

- `090-search-gene/results/gene_set_report.html`
- `090-search-gene/results/tables/gene_catalog.tsv.gz`
- `090-search-gene/results/tables/gene_expression_summary.tsv.gz`
- `090-search-gene/results/tables/stage_specificity_summary.tsv.gz`
- `090-search-gene/results/tables/expression_long.tsv.gz`
- `090-search-gene/results/tables/expression_summary_by_context.tsv.gz`
- `090-search-gene/results/tables/deg_hits.tsv.gz`
- `090-search-gene/results/tables/dtu_hits.tsv.gz` with all significant DTU hits found upstream
- `090-search-gene/results/tables/splicing_hits.tsv.gz`
- `090-search-gene/results/tables/wgcna_hits.tsv.gz`
- `090-search-gene/results/tables/mfuzz_hits.tsv.gz`
- `090-search-gene/results/plots/`
- `090-search-gene/results/groups/<group>/`
- `090-search-gene/results/genes/<group>/<gene>/`
