#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

has_flag <- function(flag) flag %in% args

log_info <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")

split_csv <- function(x) {
  x <- trimws(x)
  if (is.na(x) || x == "") return(character())
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

sanitize <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "unknown"
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- "unknown"
  x
}

resolve_table_file <- function(path) {
  if (path == "" || file.exists(path)) return(path)
  if (grepl("\\.gz$", path)) {
    plain <- sub("\\.gz$", "", path)
    if (file.exists(plain)) return(plain)
  } else {
    gz <- paste0(path, ".gz")
    if (file.exists(gz)) return(gz)
  }
  path
}

table_suffix <- Sys.getenv("PIPELINE_TABLE_SUFFIX", unset = "")

add_table_suffix <- function(path) {
  if (table_suffix != "" && grepl("\\.tsv$", path)) return(paste0(path, table_suffix))
  path
}

write_tsv2 <- function(df, path) readr::write_tsv(df, add_table_suffix(path), na = "")

normalize_id <- function(x, prefix = "") {
  x <- as.character(x)
  if (prefix != "") x <- gsub(paste0("^", prefix, ":"), "", x)
  x <- gsub("\\.[0-9]+$", "", x)
  x
}

metadata_file <- resolve_table_file(get_arg("--metadata", Sys.getenv("METADATA_FINAL_NEW", unset = Sys.getenv("METADATA_FINAL", unset = ""))))
quant_root <- get_arg("--quant-root", Sys.getenv("QUANT_DIR", unset = ""))
tx2gene_file <- resolve_table_file(get_arg("--tx2gene", Sys.getenv("TX2GENE_FILE", unset = file.path(Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"), Sys.getenv("TX2GENE_NAME", unset = "tx2gene.tsv")))))
output_root <- get_arg("--output-root", Sys.getenv("DTU_DIR", unset = "."))
projects_arg <- get_arg("--projects", "auto")
test_variables <- split_csv(get_arg("--test-variables", Sys.getenv("DTU_TEST_VARIABLES", unset = "condition,stage,sex,tissue,infection_mode")))
min_replicates <- as.integer(get_arg("--min-replicates", Sys.getenv("DTU_MIN_REPLICATES", unset = "2")))
min_gene_count <- as.numeric(get_arg("--min-gene-count", Sys.getenv("DTU_MIN_GENE_COUNT", unset = "10")))
min_transcripts <- as.integer(get_arg("--min-transcripts-per-gene", Sys.getenv("DTU_MIN_TRANSCRIPTS_PER_GENE", unset = "2")))
include_all <- has_flag("--include-all")
allow_missing <- has_flag("--allow-missing")

required_paths <- c(metadata = metadata_file, quant_root = quant_root, tx2gene = tx2gene_file)
missing_args <- names(required_paths)[required_paths == ""]
if (length(missing_args) > 0) {
  stop("[ERRO] Argumentos/caminhos ausentes: ", paste(missing_args, collapse = ", "))
}
if (!file.exists(metadata_file)) stop("[ERRO] Metadata nao encontrado: ", metadata_file)
if (!dir.exists(quant_root)) stop("[ERRO] QUANT_DIR nao encontrado: ", quant_root)
if (!file.exists(tx2gene_file)) stop("[ERRO] tx2gene nao encontrado: ", tx2gene_file)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

metadata <- readr::read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
required_cols <- c("dataset", "sample_id")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("[ERRO] Metadata sem colunas obrigatorias: ", paste(missing_cols, collapse = ", "))
}

metadata <- metadata %>%
  filter(!is.na(dataset), dataset != "", !is.na(sample_id), sample_id != "") %>%
  distinct(dataset, sample_id, .keep_all = TRUE)

if (nrow(metadata) == 0) stop("[ERRO] Nenhuma amostra valida no metadata.")

if (projects_arg == "auto") {
  projects <- sort(unique(metadata$dataset))
} else {
  projects <- split_csv(projects_arg)
}

tx2gene_raw <- readr::read_tsv(tx2gene_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
tx_col <- intersect(c("transcript_id", "TXNAME", "tx", "transcript"), colnames(tx2gene_raw))[1]
gene_col <- intersect(c("gene_id", "GENEID", "gene"), colnames(tx2gene_raw))[1]
if (is.na(tx_col) || is.na(gene_col)) {
  stop("[ERRO] tx2gene precisa ter colunas transcript_id/gene_id ou TXNAME/GENEID.")
}

tx2gene <- tx2gene_raw %>%
  transmute(
    transcript_id = normalize_id(.data[[tx_col]], "transcript"),
    gene_id = normalize_id(gsub("^gene:", "", .data[[gene_col]]), "")
  ) %>%
  filter(!is.na(transcript_id), transcript_id != "", !is.na(gene_id), gene_id != "") %>%
  distinct()

if (nrow(tx2gene) == 0) stop("[ERRO] tx2gene sem relacoes transcript-gene validas.")

empty_results <- tibble(
  analysis_id = character(),
  scope = character(),
  project = character(),
  variable = character(),
  gene_id = character(),
  transcript_id = character(),
  n_samples = integer(),
  n_levels = integer(),
  min_level_n = integer(),
  max_level = character(),
  min_level = character(),
  max_mean_usage = double(),
  min_mean_usage = double(),
  delta_usage = double(),
  pvalue = double(),
  padj = double(),
  method = character()
)

read_quant <- function(row) {
  quant_file <- file.path(quant_root, row$dataset, row$sample_id, "quant.sf")
  if (!file.exists(quant_file)) {
    if (allow_missing) {
      warning("[WARN] Ignorando quant.sf ausente: ", quant_file)
      return(tibble())
    }
    stop("[ERRO] quant.sf ausente: ", quant_file)
  }
  readr::read_tsv(
    quant_file,
    show_col_types = FALSE,
    col_types = cols(
      Name = col_character(),
      TPM = col_double(),
      NumReads = col_double(),
      .default = col_skip()
    )
  ) %>%
    transmute(
      dataset = row$dataset,
      sample_id = row$sample_id,
      import_id = row$import_id,
      transcript_id = normalize_id(Name, "transcript"),
      transcript_count = NumReads,
      transcript_tpm = TPM
    ) %>%
    inner_join(tx2gene, by = "transcript_id")
}

run_variable_test <- function(usage_tbl, variable) {
  if (!variable %in% colnames(usage_tbl)) return(empty_results)

  test_tbl <- usage_tbl %>%
    mutate(level = as.character(.data[[variable]])) %>%
    filter(!is.na(level), level != "", !tolower(level) %in% c("na", "nan", "none", "unknown"))

  level_counts <- test_tbl %>%
    distinct(import_id, level) %>%
    count(level, name = "n")

  valid_levels <- level_counts %>% filter(n >= min_replicates) %>% pull(level)
  if (length(valid_levels) < 2) return(empty_results)

  test_tbl <- test_tbl %>% filter(level %in% valid_levels)

  test_tbl %>%
    group_by(gene_id, transcript_id) %>%
    group_modify(function(df, key) {
      level_summary <- df %>%
        group_by(level) %>%
        summarise(mean_usage = mean(usage, na.rm = TRUE), n = n_distinct(import_id), .groups = "drop") %>%
        arrange(desc(mean_usage))
      if (nrow(level_summary) < 2 || min(level_summary$n, na.rm = TRUE) < min_replicates) {
        return(tibble())
      }
      pvalue <- tryCatch(
        stats::kruskal.test(usage ~ level, data = df)$p.value,
        error = function(e) NA_real_
      )
      tibble(
        variable = variable,
        n_samples = n_distinct(df$import_id),
        n_levels = nrow(level_summary),
        min_level_n = min(level_summary$n, na.rm = TRUE),
        max_level = level_summary$level[1],
        min_level = level_summary$level[nrow(level_summary)],
        max_mean_usage = level_summary$mean_usage[1],
        min_mean_usage = level_summary$mean_usage[nrow(level_summary)],
        delta_usage = max_mean_usage - min_mean_usage,
        pvalue = pvalue
      )
    }) %>%
    ungroup() %>%
    mutate(padj = p.adjust(pvalue, method = "BH"))
}

run_scope <- function(scope, project = "all_projects") {
  if (scope == "project") {
    sample_meta <- metadata %>% filter(dataset == project)
    analysis_id <- sanitize(project)
    out_dir <- file.path(output_root, project)
  } else {
    sample_meta <- metadata
    analysis_id <- "all_projects"
    out_dir <- file.path(output_root, "all_projects")
  }

  sample_meta <- sample_meta %>%
    arrange(dataset, sample_id) %>%
    mutate(import_id = if (.env$scope == "project") sample_id else paste(dataset, sample_id, sep = "__"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(sample_meta) == 0) {
    warning("[WARN] Sem amostras para escopo ", analysis_id)
    return(invisible(NULL))
  }

  log_info(paste0("DTU scope: ", analysis_id, " (", nrow(sample_meta), " amostras)"))

  quant_long <- purrr::map_dfr(seq_len(nrow(sample_meta)), function(i) read_quant(sample_meta[i, , drop = FALSE]))
  if (nrow(quant_long) == 0) {
    warning("[WARN] Nenhum quant.sf importado para ", analysis_id)
    write_tsv2(empty_results, file.path(out_dir, "dtu_results.tsv"))
    write_tsv2(empty_results, file.path(out_dir, "dtu_significant.tsv"))
    write_tsv2(tibble(), file.path(out_dir, "transcript_usage_long.tsv"))
    write_tsv2(tibble(analysis_id = analysis_id, scope = scope, project = project, n_samples = nrow(sample_meta), n_genes_tested = 0, n_transcripts_tested = 0, n_tests = 0, n_significant = 0), file.path(out_dir, "dtu_summary.tsv"))
    return(invisible(NULL))
  }

  gene_counts <- quant_long %>%
    group_by(import_id, gene_id) %>%
    summarise(gene_total_count = sum(transcript_count, na.rm = TRUE), .groups = "drop")

  eligible_genes <- gene_counts %>%
    group_by(gene_id) %>%
    summarise(total_gene_count = sum(gene_total_count, na.rm = TRUE), .groups = "drop") %>%
    inner_join(tx2gene %>% count(gene_id, name = "n_transcripts"), by = "gene_id") %>%
    filter(total_gene_count >= min_gene_count, n_transcripts >= min_transcripts)

  usage_long <- quant_long %>%
    inner_join(gene_counts, by = c("import_id", "gene_id")) %>%
    semi_join(eligible_genes, by = "gene_id") %>%
    mutate(usage = ifelse(gene_total_count > 0, transcript_count / gene_total_count, NA_real_)) %>%
    left_join(sample_meta, by = c("dataset", "sample_id", "import_id")) %>%
    mutate(analysis_id = .env$analysis_id, scope = .env$scope, project = .env$project) %>%
    select(analysis_id, scope, project, dataset, sample_id, import_id, gene_id, transcript_id, transcript_count, transcript_tpm, gene_total_count, usage, everything())

  write_tsv2(usage_long, file.path(out_dir, "transcript_usage_long.tsv"))

  test_results <- if (length(test_variables) > 0) {
    purrr::map_dfr(test_variables, function(variable) run_variable_test(usage_long, variable))
  } else {
    empty_results
  }

  results <- if (nrow(test_results) > 0) {
    test_results %>%
      mutate(
        analysis_id = .env$analysis_id,
        scope = .env$scope,
        project = .env$project,
        method = "kruskal_transcript_usage_screen"
      ) %>%
      select(analysis_id, scope, project, variable, gene_id, transcript_id, n_samples, n_levels, min_level_n, max_level, min_level, max_mean_usage, min_mean_usage, delta_usage, pvalue, padj, method) %>%
      arrange(padj, pvalue, desc(delta_usage))
  } else {
    empty_results
  }

  significant <- results %>%
    filter(!is.na(padj), padj < 0.05, !is.na(delta_usage), delta_usage >= 0.10)

  summary <- tibble(
    analysis_id = analysis_id,
    scope = scope,
    project = project,
    n_samples = nrow(sample_meta),
    n_genes_tested = n_distinct(results$gene_id),
    n_transcripts_tested = n_distinct(results$transcript_id),
    n_tests = nrow(results),
    n_significant = nrow(significant),
    min_replicates = min_replicates,
    min_gene_count = min_gene_count,
    min_transcripts_per_gene = min_transcripts,
    tested_variables = paste(test_variables, collapse = ",")
  )

  write_tsv2(results, file.path(out_dir, "dtu_results.tsv"))
  write_tsv2(significant, file.path(out_dir, "dtu_significant.tsv"))
  write_tsv2(summary, file.path(out_dir, "dtu_summary.tsv"))
  log_info(paste0("DTU concluido: ", analysis_id, " (", nrow(significant), " transcritos significativos)"))
}

for (project in projects) {
  run_scope("project", project)
}

if (include_all) {
  run_scope("all_projects", "all_projects")
}

log_info("Etapa DTU concluida.")
