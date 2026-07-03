# Resultados possiveis e interpretacao

O pipeline gera resultados em diferentes niveis: qualidade dos dados, arquivos
intermediarios processados, matrizes de expressao, analise diferencial e
relatorios interpretativos.

## Relatorios de qualidade

FastQC avalia a qualidade dos reads antes e depois do trimming. MultiQC junta
os relatorios em uma visao unica.

Perguntas que esses resultados ajudam a responder:

- A qualidade das bases e aceitavel?
- Existem adaptadores ou sequencias super-representadas?
- O trimming melhorou os reads?
- Alguma amostra parece problematica?

Resultados esperados:

- relatorios FastQC por arquivo;
- relatorio MultiQC por projeto.

## Quantificacao

No modo Salmon, cada amostra gera um `quant.sf`. Esse arquivo contem estimativas
de abundancia por transcrito. Depois, o pipeline usa tximport para sumarizar os
dados em nivel de gene.

No modo STAR, cada amostra gera alinhamentos no genoma e uma tabela
`ReadsPerGene.out.tab`, que contem contagens por gene.

## Matrizes de expressao

As matrizes principais sao:

- `counts_matrix`: contagens por gene e amostra;
- `tpm_matrix`: expressao normalizada em TPM, quando o modo e Salmon;
- `star_cpm_matrix`: expressao normalizada em CPM, quando o modo e STAR;
- `quant_samples`: tabela de amostras associada as matrizes.

## Counts, TPM e CPM

Counts:

- usados para DESeq2;
- representam contagens brutas ou estimadas por gene;
- preservam informacao necessaria para modelagem estatistica.

TPM:

- usado para comparar abundancia relativa dentro e entre amostras;
- vem do modo Salmon;
- adequado para visualizacao e relatorio de genes.

CPM:

- counts per million;
- usado como medida normalizada no modo STAR;
- adequado para visualizacao e relatorio.

Regra pratica:

- Para expressao diferencial com DESeq2, use counts.
- Para visualizacao de expressao, use TPM ou CPM.

## Correcao de lote

A etapa `055-batch-correction` pode avaliar e corrigir efeitos de lote. O
batch padrao e `dataset`, mas pode ser outra coluna configurada.

Resultados possiveis:

- matriz corrigida;
- tabela de amostras usada na correcao;
- relatorio JSON;
- PCA antes/depois;
- metricas de associacao entre PCs e lote.

Interpretacao:

Uma correcao bem-sucedida reduz agrupamentos por lote sem apagar sinal
biologico real. Se lote e condicao biologica estiverem confundidos, a correcao
pode remover sinal biologico relevante.

## DESeq2 e expressao diferencial

A etapa `060-deg-analysis` gera contrastes para variaveis biologicas presentes
nos metadados. Exemplos de variaveis:

- `condition`
- `stage`
- `sex`
- `tissue`
- `infection_mode`

Resultados principais:

- `deg_summary.tsv`: resumo das analises e contrastes;
- `DEGs_all_results.tsv`: todos os resultados testados;
- `DEGs_significant.tsv`: genes significativos;
- `contrasts/DEG_<contrast>.tsv`: resultado por contraste;
- `normalized_counts_<variable>.tsv`: contagens normalizadas;
- graficos PCA;
- heatmaps;
- volcano plots.

Colunas comuns em resultados DEG:

- `gene_id`: identificador do gene;
- `log2FoldChange`: magnitude e direcao da diferenca;
- `pvalue`: valor p;
- `padj`: valor p ajustado;
- `contrast`: comparacao testada;
- `biotype` ou anotacoes, quando disponiveis.

Interpretacao basica:

- `log2FoldChange > 0`: gene mais expresso no grupo de referencia do contraste
  definido pelo resultado.
- `log2FoldChange < 0`: gene menos expresso nesse grupo.
- `padj` baixo indica maior evidencia estatistica apos correcao de multiplos
  testes.

## Graficos DEG

PCA:

Mostra agrupamento global das amostras. Ajuda a ver separacao por condicao,
fase, tecido, sexo ou lote.

Volcano plot:

Mostra magnitude da mudanca versus significancia. Genes mais relevantes tendem
a aparecer longe do centro e no topo.

MA plot:

Mostra mudanca de expressao em funcao da abundancia media.

Heatmap:

Mostra padroes de expressao de genes selecionados entre amostras.

## Relatorio de genes candidatos

A etapa `090-search-gene` produz um HTML interativo que junta:

- matriz de expressao;
- metadados;
- anotacao GFF3/GTF;
- resultados DEG;
- genes ou grupos definidos em `genes.txt`.

Saidas importantes:

- `gene_set_report.html`
- `gene_catalog.tsv`
- `gene_expression_summary.tsv`
- `expression_long.tsv`
- `expression_summary_by_context.tsv`
- `deg_hits.tsv`
- graficos por grupo;
- graficos por gene.

O relatorio ajuda a responder:

- O gene candidato esta presente na matriz?
- Em quais amostras ou contextos ele e mais expresso?
- Ele aparece como diferencialmente expresso?
- Qual e sua anotacao, biotipo e localizacao genomica?
- Como genes de um mesmo grupo se comportam juntos?

## Resultados mais importantes para apresentar

Para uma apresentacao curta, priorize:

1. Diagrama do pipeline.
2. MultiQC resumido.
3. PCA geral.
4. Matriz de expressao ou heatmap.
5. Volcano plot de um contraste relevante.
6. Tabela de genes diferencialmente expressos.
7. Relatorio de genes candidatos, se houver foco em familias ou genes
   especificos.

