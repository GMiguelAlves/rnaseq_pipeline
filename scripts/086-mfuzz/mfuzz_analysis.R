#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
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

metadata_file <- resolve_table_file(get_arg("--metadata", Sys.getenv("METADATA_FINAL_NEW", unset = Sys.getenv("METADATA_FINAL", unset = ""))))
expression_file <- resolve_table_file(get_arg("--expression", Sys.getenv("EXPRESSION_MATRIX_FILE", unset = "")))
samples_file <- resolve_table_file(get_arg("--samples", Sys.getenv("QUANT_SAMPLES_FILE", unset = "")))
quantification_dir <- get_arg("--quantification-dir", Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"))
method <- tolower(get_arg("--method", Sys.getenv("QUANT_METHOD", unset = "salmon")))
output_root <- get_arg("--output-root", Sys.getenv("MFUZZ_DIR", unset = "."))
projects_arg <- get_arg("--projects", "auto")
include_all <- has_flag("--include-all")
allow_missing <- has_flag("--allow-missing")
time_variable <- get_arg("--time-variable", Sys.getenv("MFUZZ_TIME_VARIABLE", unset = "stage"))
time_levels <- split_csv(get_arg("--time-levels", Sys.getenv("MFUZZ_TIME_LEVELS", unset = "")))
group_columns <- split_csv(get_arg("--group-columns", Sys.getenv("MFUZZ_GROUP_COLUMNS", unset = "")))
clusters_requested <- as.integer(get_arg("--clusters", Sys.getenv("MFUZZ_CLUSTERS", unset = "6")))
m_requested <- as.numeric(get_arg("--m", Sys.getenv("MFUZZ_M", unset = "0")))
min_samples <- as.integer(get_arg("--min-samples", Sys.getenv("MFUZZ_MIN_SAMPLES", unset = "6")))
min_genes <- as.integer(get_arg("--min-genes", Sys.getenv("MFUZZ_MIN_GENES", unset = "50")))
min_expression <- as.numeric(get_arg("--min-expression", Sys.getenv("MFUZZ_MIN_EXPRESSION", unset = "1")))
min_fraction <- as.numeric(get_arg("--min-fraction", Sys.getenv("MFUZZ_MIN_FRACTION", unset = "0.20")))

missing_packages <- c("Mfuzz", "Biobase", "e1071")[!vapply(c("Mfuzz", "Biobase", "e1071"), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "[ERRO] Pacotes R ausentes no ambiente ativo: ", paste(missing_packages, collapse = ", "),
    "\n[ERRO] O arquivo envs/r-analysis.yml foi atualizado, mas ambientes Conda existentes nao mudam automaticamente.",
    "\n[ERRO] No cluster, rode: conda env update -n ", Sys.getenv("MFUZZ_ENV", unset = "r-analysis"), " -f envs/r-analysis.yml --prune",
    "\n[ERRO] Ou instale direto: conda install -n ", Sys.getenv("MFUZZ_ENV", unset = "r-analysis"), " -c conda-forge -c bioconda bioconductor-mfuzz bioconductor-biobase r-e1071"
  )
}

# Some Mfuzz versions call Biobase/e1071 functions such as exprs() and
# cmeans() by name. Attaching them makes those functions available during
# Mfuzz execution.
library(Biobase)
library(e1071)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (ncol(df) < 2) stop("[ERRO] Matriz invalida: ", path)
  colnames(df)[1] <- "gene_id"
  df %>%
    dplyr::mutate(dplyr::across(-gene_id, ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::filter(!is.na(gene_id), gene_id != "") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

read_samples <- function(path, sample_names) {
  if (path == "" || !file.exists(path)) {
    return(tibble(import_id = sample_names, sample_id = sample_names, dataset = "unknown"))
  }
  samples <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (!"import_id" %in% colnames(samples)) {
    if (all(c("dataset", "sample_id") %in% colnames(samples)) && all(sample_names %in% paste(samples$dataset, samples$sample_id, sep = "__"))) {
      samples$import_id <- paste(samples$dataset, samples$sample_id, sep = "__")
    } else if ("sample_id" %in% colnames(samples)) {
      samples$import_id <- samples$sample_id
    } else {
      stop("[ERRO] Tabela de amostras precisa de import_id ou sample_id: ", path)
    }
  }
  samples %>% dplyr::distinct(import_id, .keep_all = TRUE)
}

merge_metadata <- function(samples) {
  if (metadata_file == "" || !file.exists(metadata_file)) return(samples)
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (!all(c("dataset", "sample_id") %in% colnames(metadata))) return(samples)
  metadata$import_id_combined <- paste(metadata$dataset, metadata$sample_id, sep = "__")
  key <- if (all(samples$import_id %in% metadata$import_id_combined)) "import_id_combined" else "sample_id"
  extra <- metadata[match(samples$import_id, metadata[[key]]), , drop = FALSE]
  add_cols <- setdiff(colnames(extra), colnames(samples))
  dplyr::bind_cols(samples, extra[, add_cols, drop = FALSE])
}

project_expression_name <- function(project) {
  suffix <- if (method == "star") "star_cpm_matrix" else "tpm_matrix"
  resolve_table_file(file.path(quantification_dir, paste0(project, "_", suffix, ".tsv", table_suffix)))
}

project_samples_name <- function(project) {
  resolve_table_file(file.path(quantification_dir, paste0(project, "_quant_samples.tsv", table_suffix)))
}

load_scope <- function(scope, project = "all_projects") {
  if (scope == "project") {
    matrix_path <- project_expression_name(project)
    sample_path <- project_samples_name(project)
  } else {
    matrix_path <- expression_file
    sample_path <- samples_file
  }
  if (!file.exists(matrix_path) || !file.exists(sample_path)) {
    msg <- paste("[WARN] Arquivos ausentes para", project, ":", matrix_path, sample_path)
    if (allow_missing) {
      warning(msg)
      return(NULL)
    }
    stop(msg)
  }
  expr <- read_matrix(matrix_path)
  sample_names <- setdiff(colnames(expr), "gene_id")
  samples <- read_samples(sample_path, sample_names) %>% merge_metadata()
  key <- if (all(sample_names %in% samples$import_id)) "import_id" else if ("sample_id" %in% colnames(samples) && all(sample_names %in% samples$sample_id)) "sample_id" else ""
  if (key == "") stop("[ERRO] Amostras da matriz nao batem com tabela de amostras: ", matrix_path)
  samples <- samples[match(sample_names, samples[[key]]), , drop = FALSE]
  samples$import_id <- sample_names
  list(expr = expr, samples = samples, matrix_path = matrix_path, sample_path = sample_path)
}

write_skip <- function(out_dir, analysis_id, scope, project, reason, n_samples = 0, n_genes = 0, n_contexts = 0) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- tibble(
    analysis_id = analysis_id,
    scope = scope,
    project = project,
    status = "skipped",
    reason = reason,
    n_samples = n_samples,
    n_genes_input = n_genes,
    n_genes_used = 0,
    n_contexts = n_contexts,
    clusters = 0,
    m = NA_real_,
    time_variable = time_variable,
    group_columns = paste(group_columns, collapse = ",")
  )
  write_tsv2(summary, file.path(out_dir, "mfuzz_summary.tsv"))
  write_tsv2(tibble(gene_id = character(), cluster = integer(), membership = numeric()), file.path(out_dir, "mfuzz_membership.tsv"))
  write_tsv2(tibble(), file.path(out_dir, "cluster_centers.tsv"))
  write_tsv2(tibble(), file.path(out_dir, "expression_context_matrix.tsv"))
}

make_contexts <- function(samples) {
  if (!time_variable %in% colnames(samples)) {
    stop("[ERRO] Variavel temporal ausente na metadata: ", time_variable)
  }
  samples$time_value <- as.character(samples[[time_variable]])
  samples$time_value[is.na(samples$time_value) | samples$time_value == ""] <- "unknown"
  if (length(time_levels) > 0) {
    samples$time_value <- factor(samples$time_value, levels = unique(c(time_levels, setdiff(unique(samples$time_value), time_levels))))
  } else {
    samples$time_value <- factor(samples$time_value, levels = unique(samples$time_value))
  }
  group_present <- intersect(group_columns, colnames(samples))
  if (length(group_present) > 0) {
    for (col in group_present) {
      samples[[col]][is.na(samples[[col]]) | samples[[col]] == ""] <- "unknown"
    }
    group_text <- do.call(paste, c(samples[group_present], sep = " | "))
    samples$context <- paste(group_text, samples$time_value, sep = " | ")
  } else {
    samples$context <- as.character(samples$time_value)
  }
  samples$context <- factor(samples$context, levels = unique(samples$context))
  samples
}

plot_cluster_centers <- function(center_long, out_file) {
  if (nrow(center_long) == 0) return(FALSE)
  p <- ggplot(center_long, aes(x = context, y = center, group = cluster, color = cluster)) +
    geom_hline(yintercept = 0, color = "gray85") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Mfuzz cluster centers", x = "Context", y = "Standardized expression")
  ggsave(out_file, p, width = max(8, min(16, length(unique(center_long$context)) * 0.35 + 4)), height = max(5, min(14, length(unique(center_long$cluster)) * 1.4 + 2)), dpi = 300)
  TRUE
}

run_scope <- function(scope, project = "all_projects") {
  analysis_id <- if (scope == "project") sanitize(project) else "all_projects"
  out_dir <- file.path(output_root, analysis_id)
  dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  loaded <- load_scope(scope, project)
  if (is.null(loaded)) return(invisible(NULL))

  expr <- loaded$expr
  samples <- make_contexts(loaded$samples)
  sample_cols <- setdiff(colnames(expr), "gene_id")
  expr_mat <- as.matrix(expr[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- expr$gene_id
  storage.mode(expr_mat) <- "numeric"

  n_samples <- ncol(expr_mat)
  n_genes_input <- nrow(expr_mat)
  if (n_samples < min_samples) {
    write_skip(out_dir, analysis_id, scope, project, paste0("n_samples < ", min_samples), n_samples, n_genes_input, 0)
    return(invisible(NULL))
  }

  contexts <- levels(samples$context)
  if (length(contexts) < 2) {
    write_skip(out_dir, analysis_id, scope, project, "fewer than 2 time/context levels", n_samples, n_genes_input, length(contexts))
    return(invisible(NULL))
  }

  min_present <- max(1, ceiling(min_fraction * n_samples))
  keep <- rowSums(expr_mat >= min_expression, na.rm = TRUE) >= min_present
  keep <- keep & apply(expr_mat, 1, function(x) stats::var(x, na.rm = TRUE) > 0)
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (nrow(expr_mat) < min_genes) {
    write_skip(out_dir, analysis_id, scope, project, paste0("n_genes_used < ", min_genes), n_samples, n_genes_input, length(contexts))
    return(invisible(NULL))
  }

  expr_long <- as_tibble(expr_mat, rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "import_id", values_to = "expression") %>%
    left_join(samples %>% select(import_id, context), by = "import_id") %>%
    group_by(gene_id, context) %>%
    summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")

  context_matrix <- expr_long %>%
    mutate(context = factor(context, levels = contexts)) %>%
    arrange(gene_id, context) %>%
    pivot_wider(names_from = context, values_from = expression, values_fill = 0)
  context_mat <- as.matrix(context_matrix[, -1, drop = FALSE])
  rownames(context_mat) <- context_matrix$gene_id
  storage.mode(context_mat) <- "numeric"
  context_mat <- context_mat[apply(context_mat, 1, function(x) stats::var(x, na.rm = TRUE) > 0), , drop = FALSE]

  if (nrow(context_mat) < min_genes) {
    write_skip(out_dir, analysis_id, scope, project, "context matrix left too few variable genes", n_samples, n_genes_input, ncol(context_mat))
    return(invisible(NULL))
  }

  write_tsv2(as_tibble(context_mat, rownames = "gene_id"), file.path(out_dir, "expression_context_matrix.tsv"))

  eset <- Biobase::ExpressionSet(assayData = context_mat)
  eset <- Mfuzz::standardise(eset)
  m_value <- m_requested
  if (is.na(m_value) || m_value <= 0) {
    m_value <- Mfuzz::mestimate(eset)
  }
  clusters <- max(2, min(clusters_requested, nrow(context_mat) - 1))
  cl <- Mfuzz::mfuzz(eset, c = clusters, m = m_value)

  membership <- as.data.frame(cl$membership)
  colnames(membership) <- paste0("cluster_", seq_len(ncol(membership)))
  membership_tbl <- as_tibble(membership, rownames = "gene_id") %>%
    mutate(
      cluster = max.col(as.matrix(dplyr::select(., starts_with("cluster_"))), ties.method = "first"),
      membership = apply(as.matrix(dplyr::select(., starts_with("cluster_"))), 1, max)
    ) %>%
    relocate(cluster, membership, .after = gene_id) %>%
    arrange(cluster, desc(membership), gene_id)
  write_tsv2(membership_tbl, file.path(out_dir, "mfuzz_membership.tsv"))

  centers <- as.data.frame(cl$centers)
  colnames(centers) <- colnames(context_mat)
  center_long <- as_tibble(centers, rownames = "cluster") %>%
    mutate(cluster = paste0("cluster_", cluster)) %>%
    pivot_longer(-cluster, names_to = "context", values_to = "center") %>%
    mutate(context = factor(context, levels = colnames(context_mat))) %>%
    arrange(cluster, context)
  write_tsv2(center_long, file.path(out_dir, "cluster_centers.tsv"))
  plot_cluster_centers(center_long, file.path(out_dir, "plots", "cluster_centers.png"))

  cluster_summary <- membership_tbl %>%
    group_by(cluster) %>%
    summarise(
      n_genes = n(),
      median_membership = median(membership, na.rm = TRUE),
      high_membership_genes = sum(membership >= 0.7, na.rm = TRUE),
      .groups = "drop"
    )
  write_tsv2(cluster_summary, file.path(out_dir, "cluster_summary.tsv"))

  summary <- tibble(
    analysis_id = analysis_id,
    scope = scope,
    project = project,
    status = "completed",
    reason = "",
    n_samples = n_samples,
    n_genes_input = n_genes_input,
    n_genes_used = nrow(context_mat),
    n_contexts = ncol(context_mat),
    clusters = clusters,
    m = m_value,
    time_variable = time_variable,
    group_columns = paste(group_columns, collapse = ","),
    matrix_file = loaded$matrix_path,
    sample_file = loaded$sample_path
  )
  write_tsv2(summary, file.path(out_dir, "mfuzz_summary.tsv"))
  log_info(paste0("Mfuzz concluido: ", analysis_id, " (", clusters, " clusters)"))
}

if (projects_arg == "auto") {
  if (metadata_file != "" && file.exists(metadata_file)) {
    metadata_projects <- readr::read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
      filter(!is.na(dataset), dataset != "") %>%
      pull(dataset) %>%
      unique() %>%
      sort()
    projects <- metadata_projects
  } else {
    projects <- character()
  }
} else {
  projects <- split_csv(projects_arg)
}

for (project in projects) run_scope("project", project)
if (include_all) run_scope("all_projects", "all_projects")
log_info("Etapa Mfuzz concluida.")
