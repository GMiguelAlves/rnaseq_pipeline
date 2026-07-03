# Resumo executivo do pipeline

Este repositorio implementa um pipeline modular de RNA-seq para transformar
reads brutos em resultados biologicos interpretaveis. Ele pode rodar em
ambiente local ou em Slurm/HPC.

O fluxo principal comeca com arquivos FASTQ de projetos ENA/SRA e uma tabela
de metadados das amostras. Em seguida, o pipeline faz controle de qualidade,
limpeza dos reads, quantificacao da expressao genica, montagem de matrizes,
analise diferencial e, opcionalmente, correcao de lote e relatorio de genes
candidatos.

## Ideia central

O pipeline responde a pergunta:

> Como ir de reads RNA-seq brutos para matrizes de expressao, genes
> diferencialmente expressos e relatorios interpretaveis por amostra, grupo ou
> gene?

## Fluxo em uma frase

FASTQ + metadados + referencia -> controle de qualidade -> Salmon ou STAR ->
matrizes de expressao -> DESeq2 -> graficos, tabelas e relatorio de genes.

## Entradas principais

- Projetos ENA/SRA.
- Arquivos FASTQ brutos.
- Metadados das amostras.
- Genoma, transcriptoma e anotacao GFF3/GTF.
- Arquivo de configuracao do pipeline.

## Saidas principais

- Relatorios FastQC e MultiQC.
- FASTQs limpos e unidos por amostra.
- Arquivos `quant.sf`, quando o metodo e Salmon.
- Arquivos `ReadsPerGene.out.tab` e BAM, quando o metodo e STAR.
- Matriz de contagens por gene.
- Matriz TPM, no modo Salmon.
- Matriz CPM, no modo STAR.
- Tabela de amostras usada na quantificacao.
- Resultados de expressao diferencial com DESeq2.
- Graficos como PCA, heatmap, MA plot e volcano plot.
- Relatorio HTML de genes candidatos.

## O que torna o pipeline modular

Cada etapa tem um diretorio numerado e scripts proprios:

- `010-reference`: preparo da referencia.
- `020-data-download`: download de FASTQs.
- `025-parse`: metadados e MetaQC.
- `030-qc-fastq`: qualidade, trimming e MultiQC.
- `040-alignment`: Salmon ou STAR.
- `050-quantification`: importacao e matrizes.
- `055-batch-correction`: avaliacao/correcao de lote opcional.
- `060-deg-analysis`: analise diferencial.
- `090-search-gene`: relatorio de genes candidatos.

## Mensagem curta para apresentacao

Este pipeline organiza uma analise completa de RNA-seq, desde dados publicos
brutos ate resultados biologicos. Ele padroniza metadados, avalia a qualidade
dos reads, quantifica expressao por Salmon ou STAR, gera matrizes integradas,
executa DESeq2 e produz relatorios com tabelas e graficos.

