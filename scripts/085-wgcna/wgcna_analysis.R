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
output_root <- get_arg("--output-root", Sys.getenv("WGCNA_DIR", unset = "."))
projects_arg <- get_arg("--projects", "auto")
include_all <- has_flag("--include-all")
allow_missing <- has_flag("--allow-missing")
trait_columns <- split_csv(get_arg("--trait-columns", Sys.getenv("WGCNA_TRAIT_COLUMNS", unset = "condition,stage,sex,tissue,batch,dataset")))
min_samples <- as.integer(get_arg("--min-samples", Sys.getenv("WGCNA_MIN_SAMPLES", unset = "12")))
min_genes <- as.integer(get_arg("--min-genes", Sys.getenv("WGCNA_MIN_GENES", unset = "500")))
min_expression <- as.numeric(get_arg("--min-expression", Sys.getenv("WGCNA_MIN_EXPRESSION", unset = "1")))
min_fraction <- as.numeric(get_arg("--min-fraction", Sys.getenv("WGCNA_MIN_FRACTION", unset = "0.20")))
requested_power <- as.integer(get_arg("--power", Sys.getenv("WGCNA_POWER", unset = "0")))
fit_cutoff <- as.numeric(get_arg("--fit-cutoff", Sys.getenv("WGCNA_FIT_CUTOFF", unset = "0.80")))
network_type <- get_arg("--network-type", Sys.getenv("WGCNA_NETWORK_TYPE", unset = "signed"))
network_type <- gsub("-", " ", network_type, fixed = TRUE)
cor_method <- get_arg("--cor-method", Sys.getenv("WGCNA_COR_METHOD", unset = "pearson"))
min_module_size <- as.integer(get_arg("--min-module-size", Sys.getenv("WGCNA_MIN_MODULE_SIZE", unset = "30")))
merge_cut_height <- as.numeric(get_arg("--merge-cut-height", Sys.getenv("WGCNA_MERGE_CUT_HEIGHT", unset = "0.25")))
top_hubs <- as.integer(get_arg("--top-hubs", Sys.getenv("WGCNA_TOP_HUBS", unset = "30")))
max_block_size <- as.integer(get_arg("--max-block-size", Sys.getenv("WGCNA_MAX_BLOCK_SIZE", unset = "20000")))
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = Sys.getenv("THREADS", unset = "4")))

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("[ERRO] Pacote WGCNA nao encontrado no ambiente R ativo.")
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
tryCatch(WGCNA::allowWGCNAThreads(nThreads = threads), error = function(e) invisible(FALSE))

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

write_skip <- function(out_dir, analysis_id, scope, project, reason, n_samples = 0, n_genes = 0) {
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
    selected_power = NA_integer_,
    n_modules = 0
  )
  write_tsv2(summary, file.path(out_dir, "wgcna_summary.tsv"))
  write_tsv2(tibble(gene_id = character(), module_color = character(), module_label = integer()), file.path(out_dir, "module_assignments.tsv"))
  write_tsv2(tibble(), file.path(out_dir, "module_eigengenes.tsv"))
  write_tsv2(tibble(), file.path(out_dir, "module_trait_correlations.tsv"))
  write_tsv2(tibble(), file.path(out_dir, "hub_genes.tsv"))
}

build_trait_matrix <- function(samples, traits) {
  mats <- list()
  for (trait in traits) {
    if (!trait %in% colnames(samples)) next
    x <- samples[[trait]]
    if (all(is.na(x) | x == "")) next
    suppressWarnings(x_num <- as.numeric(x))
    if (sum(!is.na(x_num)) == length(x_num) && dplyr::n_distinct(x_num, na.rm = TRUE) > 1) {
      mats[[trait]] <- matrix(x_num, ncol = 1, dimnames = list(NULL, trait))
    } else {
      f <- factor(ifelse(is.na(x) | x == "", "unknown", as.character(x)))
      if (nlevels(f) < 2 || nlevels(f) > 30) next
      mm <- stats::model.matrix(~ 0 + f)
      colnames(mm) <- paste(trait, sanitize(levels(f)), sep = "__")
      mats[[trait]] <- mm
    }
  }
  if (length(mats) == 0) return(NULL)
  out <- do.call(cbind, mats)
  rownames(out) <- samples$import_id
  out
}

plot_soft_threshold <- function(fit, selected_power, out_file) {
  if (is.null(fit) || nrow(fit) == 0) return(FALSE)
  df <- fit %>% dplyr::mutate(selected = Power == selected_power)
  p <- ggplot(df, aes(x = Power, y = SFT.R.sq)) +
    geom_line(color = "#34699a") +
    geom_point(aes(color = selected), size = 2.5) +
    geom_hline(yintercept = fit_cutoff, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c(`FALSE` = "#34699a", `TRUE` = "#c63f3f"), guide = "none") +
    theme_bw(base_size = 10) +
    labs(title = "WGCNA soft-threshold selection", x = "Power", y = "Scale-free topology fit")
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
  TRUE
}

plot_module_sizes <- function(assignments, out_file) {
  df <- assignments %>%
    count(module_color, name = "n_genes") %>%
    arrange(desc(n_genes))
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = reorder(module_color, n_genes), y = n_genes, fill = module_color)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(title = "WGCNA module sizes", x = "Module", y = "Genes")
  ggsave(out_file, p, width = 7, height = max(4, min(12, nrow(df) * 0.25 + 2)), dpi = 300)
  TRUE
}

plot_trait_heatmap <- function(trait_long, out_file) {
  if (nrow(trait_long) == 0) return(FALSE)
  df <- trait_long %>% mutate(label = sprintf("%.2f\np=%.2g", correlation, pvalue))
  p <- ggplot(df, aes(x = trait, y = module, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 2.5) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030", limits = c(-1, 1), na.value = "gray90") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "WGCNA module-trait correlations", x = "Trait", y = "Module", fill = "cor")
  ggsave(out_file, p, width = max(7, min(16, length(unique(df$trait)) * 0.35 + 4)), height = max(4, min(14, length(unique(df$module)) * 0.28 + 3)), dpi = 300)
  TRUE
}

run_scope <- function(scope, project = "all_projects") {
  analysis_id <- if (scope == "project") sanitize(project) else "all_projects"
  out_dir <- file.path(output_root, analysis_id)
  dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  loaded <- load_scope(scope, project)
  if (is.null(loaded)) return(invisible(NULL))

  expr <- loaded$expr
  samples <- loaded$samples
  sample_cols <- setdiff(colnames(expr), "gene_id")
  expr_mat <- as.matrix(expr[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- expr$gene_id
  storage.mode(expr_mat) <- "numeric"

  n_samples <- ncol(expr_mat)
  n_genes_input <- nrow(expr_mat)
  if (n_samples < min_samples) {
    write_skip(out_dir, analysis_id, scope, project, paste0("n_samples < ", min_samples), n_samples, n_genes_input)
    return(invisible(NULL))
  }

  min_present <- max(1, ceiling(min_fraction * n_samples))
  keep <- rowSums(expr_mat >= min_expression, na.rm = TRUE) >= min_present
  keep <- keep & apply(expr_mat, 1, function(x) stats::var(x, na.rm = TRUE) > 0)
  expr_mat <- expr_mat[keep, , drop = FALSE]
  if (nrow(expr_mat) < min_genes) {
    write_skip(out_dir, analysis_id, scope, project, paste0("n_genes_used < ", min_genes), n_samples, n_genes_input)
    return(invisible(NULL))
  }

  dat_expr <- t(log2(expr_mat + 1))
  dat_expr <- as.data.frame(dat_expr)
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  if (!gsg$allOK) {
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    samples <- samples[gsg$goodSamples, , drop = FALSE]
  }

  if (nrow(dat_expr) < min_samples || ncol(dat_expr) < min_genes) {
    write_skip(out_dir, analysis_id, scope, project, "goodSamplesGenes left too few samples/genes", nrow(dat_expr), ncol(dat_expr))
    return(invisible(NULL))
  }

  powers <- c(1:10, seq(12, 20, 2))
  cor_fnc <- if (cor_method == "bicor") "bicor" else "cor"
  sft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = network_type,
    corFnc = cor_fnc,
    corOptions = list(use = "p"),
    verbose = 0
  )
  fit <- as_tibble(sft$fitIndices)
  selected_power <- requested_power
  if (is.na(selected_power) || selected_power <= 0) {
    candidates <- fit$Power[fit$SFT.R.sq >= fit_cutoff]
    selected_power <- if (length(candidates) > 0) candidates[1] else fit$Power[which.max(fit$SFT.R.sq)]
  }
  if (is.na(selected_power) || selected_power <= 0) selected_power <- 6
  write_tsv2(fit, file.path(out_dir, "soft_threshold.tsv"))
  plot_soft_threshold(fit, selected_power, file.path(out_dir, "plots", "soft_threshold.png"))

  tom_type <- if (network_type == "unsigned") "unsigned" else "signed"
  net <- WGCNA::blockwiseModules(
    dat_expr,
    power = selected_power,
    networkType = network_type,
    TOMType = tom_type,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    maxBlockSize = max_block_size,
    verbose = 0
  )

  module_colors <- WGCNA::labels2colors(net$colors)
  names(module_colors) <- colnames(dat_expr)
  assignments <- tibble(
    gene_id = names(module_colors),
    module_color = unname(module_colors),
    module_label = as.integer(net$colors)
  ) %>% arrange(module_color, gene_id)
  write_tsv2(assignments, file.path(out_dir, "module_assignments.tsv"))
  plot_module_sizes(assignments, file.path(out_dir, "plots", "module_sizes.png"))

  mes <- WGCNA::moduleEigengenes(dat_expr, module_colors)$eigengenes
  mes <- WGCNA::orderMEs(mes)
  eigengenes <- as_tibble(mes) %>%
    mutate(import_id = rownames(mes), .before = 1)
  write_tsv2(eigengenes, file.path(out_dir, "module_eigengenes.tsv"))

  trait_matrix <- build_trait_matrix(samples, trait_columns)
  trait_long <- tibble()
  if (!is.null(trait_matrix) && ncol(trait_matrix) > 0) {
    trait_matrix <- trait_matrix[rownames(mes), , drop = FALSE]
    cor_mat <- WGCNA::cor(mes, trait_matrix, use = "p")
    p_mat <- WGCNA::corPvalueStudent(cor_mat, nrow(dat_expr))
    trait_long <- as.data.frame(cor_mat) %>%
      rownames_to_column("module") %>%
      pivot_longer(-module, names_to = "trait", values_to = "correlation") %>%
      left_join(
        as.data.frame(p_mat) %>%
          rownames_to_column("module") %>%
          pivot_longer(-module, names_to = "trait", values_to = "pvalue"),
        by = c("module", "trait")
      ) %>%
      mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      arrange(padj, desc(abs(correlation)))
  }
  write_tsv2(trait_long, file.path(out_dir, "module_trait_correlations.tsv"))
  plot_trait_heatmap(trait_long, file.path(out_dir, "plots", "module_trait_heatmap.png"))

  hub_rows <- lapply(setdiff(sort(unique(module_colors)), "grey"), function(module) {
    genes <- names(module_colors)[module_colors == module]
    me_name <- paste0("ME", module)
    if (!me_name %in% colnames(mes) || length(genes) == 0) return(NULL)
    kme <- as.numeric(WGCNA::cor(dat_expr[, genes, drop = FALSE], mes[[me_name]], use = "p"))
    tibble(
      gene_id = genes,
      module_color = module,
      module_eigengene = me_name,
      kME = kme,
      abs_kME = abs(kme)
    ) %>%
      arrange(desc(abs_kME)) %>%
      mutate(module_rank = row_number()) %>%
      filter(module_rank <= top_hubs)
  })
  hubs <- bind_rows(hub_rows)
  if (nrow(hubs) == 0) {
    hubs <- tibble(
      gene_id = character(),
      module_color = character(),
      module_eigengene = character(),
      kME = numeric(),
      abs_kME = numeric(),
      module_rank = integer()
    )
  } else {
    hubs <- hubs %>% arrange(module_color, module_rank)
  }
  write_tsv2(hubs, file.path(out_dir, "hub_genes.tsv"))

  summary <- tibble(
    analysis_id = analysis_id,
    scope = scope,
    project = project,
    status = "completed",
    reason = "",
    n_samples = nrow(dat_expr),
    n_genes_input = n_genes_input,
    n_genes_used = ncol(dat_expr),
    selected_power = selected_power,
    n_modules = n_distinct(assignments$module_color[assignments$module_color != "grey"]),
    matrix_file = loaded$matrix_path,
    sample_file = loaded$sample_path
  )
  write_tsv2(summary, file.path(out_dir, "wgcna_summary.tsv"))
  log_info(paste0("WGCNA concluido: ", analysis_id, " (", summary$n_modules, " modulos)"))
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
log_info("Etapa WGCNA concluida.")
