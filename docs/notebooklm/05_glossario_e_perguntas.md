# Glossario e perguntas guia

Este documento ajuda uma LLM a responder perguntas sobre o pipeline com
linguagem consistente.

## Glossario

FASTQ:

Arquivo com sequencias de reads e qualidade por base.

Read:

Pequena sequencia produzida pelo sequenciamento.

Metadata:

Tabela que descreve amostras, runs, condicoes experimentais e covariaveis.

MetaQC:

Ferramenta usada aqui para validar, enriquecer, padronizar e unir metadados.

FastQC:

Ferramenta de controle de qualidade de reads.

Trim Galore:

Ferramenta para remover adaptadores e bases de baixa qualidade.

MultiQC:

Ferramenta que junta relatorios de QC em um unico relatorio.

Salmon:

Ferramenta de quantificacao de expressao baseada em transcriptoma.

STAR:

Alinhador de RNA-seq contra genoma. Neste pipeline pode gerar BAM e contagens
por gene.

tximport:

Ferramenta R usada para importar resultados Salmon e sumarizar para genes.

Counts:

Contagens por gene. Sao a entrada principal para DESeq2.

TPM:

Transcripts per million. Medida normalizada de abundancia usada para
visualizacao e comparacao exploratoria.

CPM:

Counts per million. Medida normalizada derivada de contagens.

Batch:

Lote tecnico ou agrupamento que pode introduzir variacao nao biologica.

DESeq2:

Ferramenta R/Bioconductor para analise de expressao diferencial.

DEG:

Gene diferencialmente expresso.

PCA:

Analise de componentes principais, usada para visualizar estrutura global das
amostras.

Volcano plot:

Grafico que combina magnitude de mudanca e significancia estatistica.

Heatmap:

Mapa de calor de expressao genica.

## Perguntas e respostas curtas

### O que este pipeline faz?

Ele transforma dados RNA-seq brutos em resultados interpretaveis, incluindo QC,
matrizes de expressao, DEGs e relatorio de genes candidatos.

### Quais sao as entradas minimas?

Projetos ENA/SRA ou FASTQs, metadados padronizaveis, referencia biologica e
configuracao do pipeline.

### Qual e a funcao do MetaQC?

Padronizar metadados. Ele baixa/valida metadados, aplica regras YAML, pode
enriquecer com dados de autor e gera uma tabela final combinada.

### Qual e a diferenca entre Salmon e STAR?

Salmon quantifica contra o transcriptoma e gera `quant.sf`. STAR alinha reads
ao genoma e gera BAM e contagens por gene.

### Qual matriz deve ser usada no DESeq2?

A matriz de counts.

### Para que servem TPM e CPM?

TPM e CPM sao mais adequados para visualizacao, exploracao e relatorios de
expressao.

### Quando usar correcao de lote?

Quando ha evidencia de que amostras se agrupam por lote tecnico, e quando o
lote nao esta completamente confundido com a variavel biologica de interesse.

### O que e o relatorio de genes candidatos?

Um HTML que integra expressao, metadados, anotacao e DEGs para genes ou grupos
de genes definidos pelo usuario.

## Perguntas que o NotebookLM deve conseguir responder

- Explique o pipeline para uma banca de mestrado.
- Explique o pipeline para alguem que nao conhece bioinformatica.
- Liste as entradas e saidas de cada etapa.
- Explique como o MetaQC transforma os metadados.
- Diga quais arquivos finais sao mais importantes.
- Compare Salmon e STAR neste pipeline.
- Explique por que counts sao usados no DESeq2.
- Explique por que TPM/CPM sao usados nos relatorios.
- Aponte riscos de correcao de lote.
- Sugira um roteiro de apresentacao de 5 minutos.
- Sugira uma legenda para a figura metodologica do pipeline.
- Sugira um texto curto para materiais e metodos.

## Texto curto para materiais e metodos

Os dados de RNA-seq foram processados por um pipeline modular executavel em
ambiente local ou Slurm/HPC. Inicialmente, arquivos de referencia e metadados
foram preparados e padronizados. Os arquivos FASTQ foram avaliados com FastQC,
aparados com Trim Galore, reunidos por amostra quando necessario e sumarizados
com MultiQC. A expressao genica foi quantificada por Salmon ou STAR, seguida da
geracao de matrizes de counts e de expressao normalizada. As matrizes foram
usadas para analise diferencial com DESeq2 e para a geracao de graficos,
tabelas e relatorios exploratorios de genes candidatos.

