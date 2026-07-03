# Materiais para NotebookLM

Este diretorio contem textos preparados para alimentar uma LLM de apoio, como
NotebookLM, sobre o pipeline de RNA-seq deste repositorio.

## Como usar

Carregue estes arquivos como fontes no NotebookLM:

1. `00_resumo_executivo.md`
2. `01_pipeline_visao_geral.md`
3. `02_metaqc_e_metadados.md`
4. `03_etapas_entradas_saidas.md`
5. `04_resultados_e_interpretacao.md`
6. `05_glossario_e_perguntas.md`

Opcionalmente, carregue tambem os READMEs originais do repositorio:

- `README.md`
- `025-parse/README.md`
- `030-qc-fastq/README.md`
- `040-alignment/README.md`
- `050-quantification/README.md`
- `055-batch-correction/README.md`
- `060-deg-analysis/README.md`
- `090-search-gene/README.md`

## Perguntas boas para fazer ao NotebookLM

- Explique o pipeline em linguagem simples.
- Quais sao as entradas obrigatorias do pipeline?
- O que o MetaQC faz neste projeto?
- Qual a diferenca entre o modo Salmon e o modo STAR?
- Quais resultados finais o pipeline pode gerar?
- Como interpretar counts, TPM e CPM?
- Quais arquivos eu devo mostrar em uma apresentacao?
- O que pode dar errado na etapa de metadados?
- Como a correcao de lote entra no fluxo?
- Como o relatorio de genes candidatos integra os resultados?

