# MetaQC e metadados

Os metadados sao uma parte central do pipeline. Eles ligam os arquivos FASTQ
as informacoes biologicas e experimentais das amostras. Sem metadados bem
padronizados, as etapas de quantificacao, correcao de lote, DESeq2 e relatorio
de genes ficam ambiguas ou incorretas.

## O que o MetaQC faz neste pipeline

Neste repositorio, o comando `metaqc` e usado em tres momentos:

1. `metaqc validate`
2. `metaqc enrich`
3. `metaqc parse`
4. `metaqc merge`

Na pratica, a etapa `025-parse` baixa metadados do ENA, valida campos basicos,
opcionalmente enriquece os metadados com informacoes de autor, aplica regras de
parse definidas em YAML e une os projetos em uma tabela final.

## Etapa validate

O script `scripts/025-parse/run_metaqc.sh` baixa uma tabela TSV do ENA usando o
endpoint de filereport. O resultado bruto e salvo em:

```text
025-parse/010-raw_metadata/<PROJECT>.tsv
```

Depois, o pipeline roda:

```bash
metaqc validate <PROJECT>.tsv --output <PROJECT>_base.csv
```

A saida principal e:

```text
025-parse/015-intermediate_folder/<PROJECT>_base.csv
```

Essa tabela base preserva colunas essenciais, como `sample_id`,
`run_accession` e `study_accession`, dependendo da configuracao
`METAQC_KEEP_COLUMNS`.

## Etapa enrich

A etapa `enrich` e opcional. Ela e usada quando existe:

```text
025-parse/020-metadata_parsers/<PROJECT>/configs/<PROJECT>_enrich.yaml
025-parse/020-metadata_parsers/<PROJECT>/author_metadata.tsv
```

Nesse caso, o pipeline combina metadados baixados do ENA com metadados
fornecidos pelo autor ou curador do projeto. A saida e:

```text
025-parse/015-intermediate_folder/<PROJECT>_enriched.csv
```

Essa etapa e util quando o ENA nao contem informacoes biologicas suficientes,
como fase de vida, sexo, tecido, condicao experimental ou lote.

## Etapa parse

O parse transforma informacoes textuais, nomes de amostras e colunas brutas em
campos padronizados.

Cada projeto deve ter um YAML:

```text
025-parse/020-metadata_parsers/<PROJECT>/configs/<PROJECT>.yaml
```

O YAML define:

- valores padrao;
- quais colunas preservar;
- regras de regex para extrair `stage`, `tissue`, `sex`, `condition` e outros
  campos;
- como construir `sample_id`;
- quais colunas finais devem aparecer na tabela.

A saida por projeto e:

```text
025-parse/020-metadata_parsers/Allprojects/<PROJECT>_parsed.csv
```

## Etapa merge

Depois que todos os projetos foram processados, o pipeline une os arquivos
`*_parsed.csv` com:

```bash
metaqc merge <METADATA_PARSED_DIR> --output <METADATA_FINAL>
```

As saidas finais sao:

```text
025-parse/030-metadata_final/AllProjects_metadata.csv
025-parse/030-metadata_final/AllProjects_metadata_new.csv
```

## Colunas obrigatorias

A tabela final deve conter pelo menos:

- `dataset`
- `sample_id`
- `run_accession`

Essas colunas conectam projeto, amostra biologica e corrida de sequenciamento.

## Colunas recomendadas

Para analises biologicas e estatisticas, sao recomendadas:

- `condition`
- `stage`
- `tissue`
- `sex`
- `batch`
- `replicate`

Essas colunas podem ser usadas para contrastes no DESeq2, correcao de lote,
agrupamentos em graficos e filtros no relatorio de genes.

## Por que metadados sao criticos

Metadados mal formatados podem causar:

- amostras sem FASTQ correspondente;
- runs diferentes tratados como amostras independentes;
- contrastes biologicos incorretos;
- confundimento entre lote e condicao biologica;
- graficos dificeis de interpretar;
- relatorio de genes com agrupamentos errados.

## Exemplo de contrato minimo

Uma tabela minima de metadados segue esta ideia:

```text
dataset  sample_id   run_accession  condition  stage    tissue   sex      batch  replicate
PRJXXXX  SAMPLE_001  RUN000001      control    unknown  unknown  unknown  B1     1
```

## Como explicar o MetaQC em uma frase

MetaQC e a camada que transforma metadados publicos e heterogeneos em uma
tabela padronizada, validada e pronta para guiar o restante da analise RNA-seq.

