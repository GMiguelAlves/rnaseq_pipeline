# Etapas, entradas e saidas

Este documento resume cada etapa numerada do pipeline, com suas entradas e
saidas principais.

## 010-reference: referencia

Objetivo:

Preparar os arquivos de referencia usados na quantificacao ou alinhamento.

Entradas:

- genoma FASTA;
- transcriptoma FASTA;
- anotacao GFF3/GTF;
- URLs ou arquivos locais definidos na configuracao.

Ferramentas:

- Salmon, para indice de transcriptoma;
- STAR, para indice de genoma;
- gffread, quando necessario para preparar anotacoes.

Saidas:

- `010-reference/salmon_index/`
- `010-reference/star_index/`
- `010-reference/star_index_gtf/`

## 020-data-download: download de dados

Objetivo:

Baixar FASTQs de projetos ENA/SRA definidos em `PIPELINE_PROJECTS`.

Entradas:

- configuracao do projeto;
- arquivos de link ou URLs ENA/SRA.

Saidas:

- `${SCRATCH_ROOT}/<PROJECT>/fastq_ftp/`
- manifestos de renomeacao, quando usados.

## 025-parse: metadados

Objetivo:

Baixar, validar, enriquecer, padronizar e unir metadados das amostras.

Entradas:

- accession do projeto;
- metadados ENA;
- YAML de parse por projeto;
- metadados de autor, se houver.

Saidas:

- `025-parse/010-raw_metadata/<PROJECT>.tsv`
- `025-parse/015-intermediate_folder/<PROJECT>_base.csv`
- `025-parse/015-intermediate_folder/<PROJECT>_enriched.csv`
- `025-parse/020-metadata_parsers/Allprojects/<PROJECT>_parsed.csv`
- `025-parse/030-metadata_final/AllProjects_metadata_new.csv`

## 030-qc-fastq: controle de qualidade

Objetivo:

Avaliar, limpar e organizar os reads.

Fluxo:

1. FastQC nos FASTQs brutos.
2. Trim Galore para remocao de adaptadores e bases de baixa qualidade.
3. FastQC nos reads aparados.
4. Merge de runs pertencentes a mesma amostra biologica.
5. FastQC nos FASTQs unidos por amostra.
6. MultiQC para um relatorio consolidado.

Entradas:

- FASTQs brutos;
- metadados finais;
- plano de QC gerado pelo pipeline.

Saidas:

- relatorios FastQC;
- FASTQs aparados;
- FASTQs unidos por amostra;
- relatorio MultiQC;
- `030-qc-fastq/work/<PROJECT>_qc_plan.csv`

## 040-alignment: quantificacao ou alinhamento

Objetivo:

Transformar reads limpos em medidas de expressao.

Modo Salmon:

- entrada: FASTQs limpos unidos por amostra;
- referencia: indice Salmon;
- saida: `quant.sf` por amostra.

Modo STAR:

- entrada: FASTQs limpos unidos por amostra;
- referencia: indice STAR com anotacao;
- saidas: BAM ordenado, `ReadsPerGene.out.tab` e logs STAR.

Saidas tipicas:

- `040-alignment/work/<PROJECT>_salmon_plan.csv`
- `${QUANT_DIR}/<PROJECT>/<sample_id>/quant.sf`
- `040-alignment/work/<PROJECT>_star_plan.csv`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/ReadsPerGene.out.tab`
- `${STAR_QUANT_DIR}/<PROJECT>/<sample_id>/Aligned.sortedByCoord.out.bam`

## 050-quantification: importacao e matrizes

Objetivo:

Importar resultados de Salmon ou STAR e gerar matrizes por gene.

Modo Salmon:

- importa `quant.sf` com tximport;
- gera matriz de counts;
- gera matriz TPM.

Modo STAR:

- importa `ReadsPerGene.out.tab`;
- gera matriz de counts;
- gera matriz CPM.

Saidas:

- `050-quantification/counts_matrix.tsv` ou `.tsv.gz`
- `050-quantification/tpm_matrix.tsv` ou `.tsv.gz`
- `050-quantification/star_cpm_matrix.tsv` ou `.tsv.gz`
- `050-quantification/quant_samples.tsv` ou `.tsv.gz`
- matrizes especificas por projeto.

## 055-batch-correction: lote opcional

Objetivo:

Avaliar e corrigir efeito de lote, quando apropriado.

Entradas:

- matriz de counts;
- tabela de amostras;
- coluna de lote, por padrao `dataset`;
- covariaveis biologicas a preservar, quando configuradas.

Saidas:

- `055-batch-correction/all_projects/assessment/`
- `055-batch-correction/all_projects/counts_batch_corrected.tsv`
- `055-batch-correction/all_projects/batch_correction_samples.tsv`
- PCA antes/depois da correcao, quando gerado.

Observacao:

A correcao de lote deve ser usada com cuidado quando lote e variavel biologica
estao confundidos.

## 060-deg-analysis: expressao diferencial

Objetivo:

Rodar DESeq2 para identificar genes diferencialmente expressos.

Entradas:

- matriz de counts;
- tabela de amostras;
- variaveis de teste, como `condition`, `stage`, `sex`, `tissue` ou
  `infection_mode`;
- covariaveis de desenho, se configuradas.

Saidas:

- `060-deg-analysis/work/deg_plan.csv`
- `deg_summary.tsv`
- `DEGs_all_results.tsv`
- `DEGs_significant.tsv`
- arquivos por contraste em `contrasts/`
- contagens normalizadas;
- graficos em `plots/`.

## 090-search-gene: relatorio de genes candidatos

Objetivo:

Gerar um relatorio exploratorio para genes ou grupos de genes de interesse.

Entradas:

- lista `090-search-gene/genes.txt`;
- matriz de expressao TPM ou CPM;
- tabela de amostras;
- anotacao GFF3/GTF;
- resultados de DESeq2.

Saidas:

- `090-search-gene/results/gene_set_report.html`
- `gene_catalog.tsv`
- `gene_expression_summary.tsv`
- `expression_long.tsv`
- `expression_summary_by_context.tsv`
- `deg_hits.tsv`
- graficos globais, por grupo e por gene.

