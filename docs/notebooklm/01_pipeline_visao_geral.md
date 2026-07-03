# Visao geral do pipeline RNA-seq

O pipeline foi desenhado para ser generico e reprodutivel. Ele nao depende de
um unico organismo: o usuario define organismo, projetos, referencias,
metadados e metodos no arquivo de configuracao.

## Como o pipeline e executado

O script principal e:

```bash
bash rnaseq_pipeline.sh --all
```

Para rodar localmente, sem Slurm:

```bash
bash rnaseq_pipeline.sh --all --local
```

Para apenas visualizar os comandos que seriam executados:

```bash
bash rnaseq_pipeline.sh --all --dry-run
```

O pipeline tambem aceita etapas individuais:

```bash
bash rnaseq_pipeline.sh --step metadata --step qc
```

As etapas reconhecidas incluem:

- `reference`
- `download`
- `metadata`
- `qc`
- `salmon`
- `star`
- `tximport`
- `batch`
- `deg`
- `report`

## Configuracao

O usuario normalmente edita `config/user_settings.sh`, criado a partir de
`config/user_settings_template.sh`.

Configuracoes importantes:

- `PIPELINE_PROJECTS`: projetos ENA/SRA a processar.
- `ORGANISM_NAME`: nome do organismo.
- `SCRATCH_ROOT`: area de trabalho para arquivos grandes.
- `PIPELINE_EXECUTOR`: `slurm` ou `local`.
- `QUANT_METHOD`: `salmon` ou `star`.
- Caminhos para genoma, transcriptoma e anotacao GFF3/GTF.

## Dois modos de quantificacao

O pipeline pode usar Salmon ou STAR.

No modo Salmon:

- usa o transcriptoma como referencia;
- gera um indice Salmon;
- quantifica abundancia de transcritos;
- produz `quant.sf`;
- gera matriz TPM alem da matriz de counts.

No modo STAR:

- usa o genoma anotado como referencia;
- gera um indice STAR com GTF/GFF;
- alinha reads ao genoma;
- produz BAM ordenado;
- produz `ReadsPerGene.out.tab`;
- gera matriz CPM alem da matriz de counts.

## Ordem logica das etapas

1. Preparar referencia.
2. Baixar FASTQs.
3. Baixar, validar e padronizar metadados.
4. Fazer QC dos reads.
5. Quantificar ou alinhar reads limpos.
6. Importar quantificacoes e gerar matrizes.
7. Opcionalmente avaliar/corrigir lote.
8. Rodar DESeq2.
9. Opcionalmente gerar relatorio de genes candidatos.

## Dependencias entre etapas

A referencia e necessaria para Salmon ou STAR. Os FASTQs e metadados sao
necessarios para o controle de qualidade. O QC gera os arquivos limpos usados
pela quantificacao. A quantificacao gera os arquivos que alimentam as matrizes.
As matrizes alimentam a correcao de lote e o DESeq2. O relatorio de genes usa
matrizes, metadados, anotacao e resultados DEG.

## O papel do Slurm

Em modo Slurm, o pipeline submete jobs com dependencias `afterok`. Isso permite
que etapas independentes rodem em paralelo e que etapas posteriores so iniciem
quando as anteriores forem concluidas com sucesso.

Em modo local, o pipeline executa as mesmas etapas em ordem no computador atual.

