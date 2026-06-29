#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

log_info <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")

genes_file <- get_arg("--genes", "genes.txt")
tpm_file <- get_arg("--tpm", Sys.getenv("EXPRESSION_MATRIX_FILE", unset = file.path(Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"), Sys.getenv("SALMON_TPM_MATRIX_NAME", unset = "tpm_matrix.tsv"))))
expression_unit <- get_arg("--expression-unit", Sys.getenv("EXPRESSION_UNIT", unset = "TPM"))
samples_file <- get_arg("--samples", Sys.getenv("QUANT_SAMPLES_FILE", unset = file.path(Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"), Sys.getenv("QUANT_SAMPLES_NAME", unset = "quant_samples.tsv"))))
metadata_file <- get_arg("--metadata", Sys.getenv("METADATA_FINAL_NEW", unset = Sys.getenv("METADATA_FINAL", unset = "")))
deg_root <- get_arg("--deg-root", Sys.getenv("DEG_DIR", unset = "../060-deg-analysis"))
dtu_root <- get_arg("--dtu-root", Sys.getenv("DTU_DIR", unset = "../070-dtu-analysis"))
splicing_root <- get_arg("--splicing-root", Sys.getenv("SPLICING_DIR", unset = "../080-splicing"))
batch_root <- get_arg("--batch-root", Sys.getenv("BATCH_DIR", unset = "../055-batch-correction"))
corrected_expression_file <- get_arg("--corrected-expression", Sys.getenv("BATCH_CORRECTED_COUNTS_FILE", unset = ""))
corrected_expression_unit <- get_arg("--corrected-expression-unit", Sys.getenv("BATCH_CORRECTED_EXPRESSION_UNIT", unset = "batch-corrected counts"))
wgcna_root <- get_arg("--wgcna-root", Sys.getenv("WGCNA_DIR", unset = "../085-wgcna"))
mfuzz_root <- get_arg("--mfuzz-root", Sys.getenv("MFUZZ_DIR", unset = "../086-mfuzz"))
gff_file <- get_arg("--gff", Sys.getenv("GENE_REPORT_ANNOTATION_FILE", unset = Sys.getenv("REF_GFF3", unset = "")))
out_dir <- get_arg("--output-dir", file.path(Sys.getenv("GENE_REPORT_DIR", unset = "."), "results"))
report_title <- get_arg("--title", "Relatorio exploratorio de genes")
if (is.na(expression_unit) || expression_unit == "") expression_unit <- "TPM"
expression_log_label <- paste0("log2(", expression_unit, "+1)")
expression_mean_log_label <- paste0("Media log2(", expression_unit, "+1)")

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

tpm_file <- resolve_table_file(tpm_file)
samples_file <- resolve_table_file(samples_file)
metadata_file <- resolve_table_file(metadata_file)

if (!file.exists(genes_file)) stop("[ERRO] genes.txt nao encontrado: ", genes_file)
if (!file.exists(tpm_file)) stop("[ERRO] Matriz de expressao nao encontrada: ", tpm_file)
if (!file.exists(samples_file)) warning("[WARN] Tabela de amostras nao encontrada; inferindo metadata minima pelos nomes das colunas da matriz de expressao: ", samples_file)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "genes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "groups"), recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "unknown"
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- "unknown"
  x
}

table_suffix <- Sys.getenv("PIPELINE_TABLE_SUFFIX", unset = "")

add_table_suffix <- function(path) {
  if (table_suffix != "" && grepl("\\.tsv$", path)) {
    return(paste0(path, table_suffix))
  }
  path
}

write_tsv2 <- function(df, path) readr::write_tsv(df, add_table_suffix(path), na = "")

find_batch_corrected_expression_file <- function(root) {
  if (root == "" || !dir.exists(root)) return("")
  candidates <- c(
    file.path(root, "all_projects", paste0("counts_batch_corrected.tsv", table_suffix)),
    file.path(root, "all_projects", "counts_batch_corrected.tsv"),
    file.path(root, "all_projects", "counts_batch_corrected.tsv.gz")
  )
  candidates <- unique(vapply(candidates, resolve_table_file, character(1)))
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) "" else found[[1]]
}

if (corrected_expression_file == "") {
  corrected_expression_file <- find_batch_corrected_expression_file(batch_root)
} else {
  corrected_expression_file <- resolve_table_file(corrected_expression_file)
  if (!file.exists(corrected_expression_file)) {
    warning("[WARN] Matriz corrigida informada mas nao encontrada: ", corrected_expression_file)
    corrected_expression_file <- ""
  }
}

safe_div <- function(x, y) {
  ifelse(is.na(y) | y == 0, NA_real_, x / y)
}

split_env_csv <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  value <- trimws(value)
  if (value == "") return(character())
  trimws(unlist(strsplit(value, ",")))
}

life_stage_levels <- split_env_csv("LIFE_STAGE_LEVELS", "unknown")
if (!"unknown" %in% life_stage_levels) life_stage_levels <- c(life_stage_levels, "unknown")
stage_synonym_map <- split_env_csv("STAGE_SYNONYM_MAP", "")
organism_specific_reports <- Sys.getenv("ORGANISM_SPECIFIC_REPORTS", unset = "0") %in% c("1", "true", "TRUE", "yes", "YES")
stage_tau_threshold <- suppressWarnings(as.numeric(Sys.getenv("GENE_REPORT_STAGE_TAU_THRESHOLD", unset = "0.60")))
if (is.na(stage_tau_threshold)) stage_tau_threshold <- 0.60
stage_min_expression <- suppressWarnings(as.numeric(Sys.getenv("GENE_REPORT_STAGE_MIN_EXPRESSION", unset = "1")))
if (is.na(stage_min_expression)) stage_min_expression <- 1
mfuzz_membership_threshold <- suppressWarnings(as.numeric(Sys.getenv("GENE_REPORT_MFUZZ_MEMBERSHIP_THRESHOLD", unset = "0.70")))
if (is.na(mfuzz_membership_threshold)) mfuzz_membership_threshold <- 0.70

normalize_stage_detail <- function(stage) {
  x <- tolower(trimws(as.character(stage)))
  x <- gsub("[[:space:]_-]+", "_", x)
  x[x %in% c("", "na", "nan", "none", "unknown", "not_available")] <- "unknown"
  if (length(stage_synonym_map) > 0) {
    for (rule in stage_synonym_map) {
      parts <- strsplit(rule, "=", fixed = TRUE)[[1]]
      if (length(parts) == 2 && nzchar(parts[1])) {
        x <- gsub(parts[1], parts[2], x)
      }
    }
  }
  x
}

classify_life_stage <- function(stage_detail) {
  x <- as.character(stage_detail)
  out <- rep("unknown", length(x))
  for (level in setdiff(life_stage_levels, "unknown")) {
    out[grepl(paste0("^", level, "($|_)"), x)] <- level
  }
  out
}

extract_stage_day <- function(stage_detail) {
  x <- as.character(stage_detail)
  day <- stringr::str_match(x, "(?:^|_)([0-9]+(?:\\.[0-9]+)?)(?:_)?d(?:$|_)")[, 2]
  suppressWarnings(as.numeric(day))
}

order_stage_details <- function(stage_detail) {
  details <- unique(as.character(stage_detail))
  stage_df <- tibble::tibble(
    stage = details,
    stage_class = classify_life_stage(details),
    stage_day = extract_stage_day(details)
  ) %>%
    dplyr::mutate(
      stage_class = factor(stage_class, levels = life_stage_levels),
      stage_day_sort = ifelse(is.na(stage_day), Inf, stage_day)
    ) %>%
    dplyr::arrange(stage_class, stage_day_sort, stage)
  stage_df$stage
}

clean_annotation_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tryCatch(utils::URLdecode(x), error = function(e) x)
  x <- gsub("[\t\r\n]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

shorten_annotation_label <- function(x, max_chars = 80) {
  x <- clean_annotation_text(x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 3), "...")
  x
}

is_uninformative_gene_name <- function(gene_name, gene_id) {
  gene_name <- clean_annotation_text(gene_name)
  gene_id <- as.character(gene_id)
  gene_name == "" |
    gene_name == gene_id |
    grepl("^gene:", gene_name)
}

make_gene_display_label <- function(gene_name, gene_id, description = "") {
  gene_name <- as.character(gene_name)
  gene_id <- as.character(gene_id)
  description <- as.character(description)
  label <- clean_annotation_text(gene_name)
  desc <- clean_annotation_text(description)
  use_desc <- is_uninformative_gene_name(label, gene_id) & desc != ""
  label[use_desc] <- desc[use_desc]
  label <- shorten_annotation_label(label)
  label[is.na(label) | label == ""] <- gene_id[is.na(label) | label == ""]
  ifelse(label == gene_id, gene_id, paste(label, gene_id, sep = " | "))
}

plot_or_skip <- function(label, plot_fun) {
  ok <- tryCatch(plot_fun(), error = function(e) {
    warning(label, ": ", e$message)
    FALSE
  })
  isTRUE(ok)
}

parse_gene_groups <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  rows <- list()
  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) next
    if (!grepl(":", line, fixed = TRUE)) {
      warning("Linha ignorada em genes.txt sem ':': ", line)
      next
    }
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    group <- trimws(parts[1])
    genes <- trimws(unlist(strsplit(paste(parts[-1], collapse = ":"), "[,;]")))
    genes <- genes[genes != ""]
    if (length(genes) == 0) next
    rows[[length(rows) + 1]] <- data.frame(group = group, query = genes, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) stop("[ERRO] Nenhum gene encontrado em ", path)
  dplyr::bind_rows(rows) %>% dplyr::distinct(group, query, .keep_all = TRUE)
}

read_matrix <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (ncol(df) < 2) stop("[ERRO] Matriz invalida: ", path)
  colnames(df)[1] <- "gene_id"
  df %>% dplyr::mutate(dplyr::across(-gene_id, ~ suppressWarnings(as.numeric(.x))))
}

read_samples <- function(path, sample_names) {
  if (!file.exists(path)) {
    return(tibble::tibble(
      import_id = sample_names,
      sample_id = sample_names,
      dataset = "unknown"
    ))
  }
  samples <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (!"import_id" %in% colnames(samples)) {
    if (all(c("dataset", "sample_id") %in% colnames(samples))) {
      combined <- paste(samples$dataset, samples$sample_id, sep = "__")
      samples$import_id <- if (all(sample_names %in% combined)) combined else samples$sample_id
    } else if ("sample_id" %in% colnames(samples)) {
      samples$import_id <- samples$sample_id
    } else {
      stop("[ERRO] Tabela de amostras precisa de import_id ou sample_id.")
    }
  }
  samples <- samples %>% dplyr::distinct(import_id, .keep_all = TRUE)
  missing <- setdiff(sample_names, samples$import_id)
  if (length(missing) > 0) stop("[ERRO] Amostras sem metadata: ", paste(head(missing, 20), collapse = ", "))
  samples[match(sample_names, samples$import_id), , drop = FALSE]
}

empty_annotations <- function() {
  tibble::tibble(
    gene_id = character(),
    gene_name = character(),
    biotype = character(),
    description = character(),
    chromosome = character(),
    gene_start = integer(),
    gene_end = integer(),
    strand = character(),
    location = character()
  )
}

load_annotations <- function(gff_file) {
  if (gff_file == "" || !file.exists(gff_file)) {
    return(empty_annotations())
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    warning("Pacote rtracklayer nao encontrado; seguindo sem anotacao GFF3.")
    return(empty_annotations())
  }
  log_info("Lendo anotacao GFF/GTF...")
  gff <- rtracklayer::import(gff_file)
  genes <- as.data.frame(gff[gff$type == "gene"])
  if (nrow(genes) == 0) {
    return(empty_annotations())
  }
  pick_col <- function(df, names, default = NA_character_) {
    found <- intersect(names, colnames(df))
    if (length(found) == 0) return(rep(default, nrow(df)))
    as.character(df[[found[1]]])
  }
  gene_id <- pick_col(genes, c("ID", "gene_id"))
  gene_id <- gsub("^gene:", "", gene_id)
  gene_id <- gsub("\\.[0-9]+$", "", gene_id)
  gene_name <- pick_col(genes, c("gene_name", "symbol", "gene", "Name", "locus_tag"))
  gene_name[is.na(gene_name) | gene_name == ""] <- gene_id[is.na(gene_name) | gene_name == ""]
  biotype <- pick_col(genes, c("biotype", "gene_biotype", "type"), "Unknown")
  biotype[is.na(biotype) | biotype == ""] <- "Unknown"
  description <- clean_annotation_text(pick_col(genes, c("description", "product", "Note", "note"), ""))
  chromosome <- clean_annotation_text(pick_col(genes, c("seqnames", "seqid", "chromosome", "chr"), ""))
  gene_start <- suppressWarnings(as.integer(pick_col(genes, c("start"), NA_character_)))
  gene_end <- suppressWarnings(as.integer(pick_col(genes, c("end"), NA_character_)))
  strand <- clean_annotation_text(pick_col(genes, c("strand"), ""))
  strand[is.na(strand) | strand == "*" | strand == "."] <- ""
  location <- ifelse(
    chromosome != "" & !is.na(gene_start) & !is.na(gene_end),
    paste0(chromosome, ":", gene_start, "-", gene_end, ifelse(strand != "", paste0("(", strand, ")"), "")),
    ""
  )
  tibble::tibble(
    gene_id = gene_id,
    gene_name = gene_name,
    biotype = biotype,
    description = description,
    chromosome = chromosome,
    gene_start = gene_start,
    gene_end = gene_end,
    strand = strand,
    location = location
  ) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}

build_gene_catalog <- function(gene_groups, tpm, annotations) {
  tpm_genes <- tpm$gene_id
  ann <- annotations
  gene_groups %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      matched_gene_id = dplyr::case_when(
        query %in% tpm_genes ~ query,
        query %in% ann$gene_id ~ query,
        query %in% ann$gene_name ~ ann$gene_id[match(query, ann$gene_name)],
        TRUE ~ query
      ),
      match_type = dplyr::case_when(
        query %in% tpm_genes ~ "gene_id",
        query %in% ann$gene_id ~ "annotation_gene_id",
        query %in% ann$gene_name ~ "gene_name",
        TRUE ~ "unmatched"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(ann, by = c("matched_gene_id" = "gene_id")) %>%
    dplyr::mutate(
      gene_name = ifelse(is.na(gene_name) | gene_name == "", matched_gene_id, gene_name),
      biotype = ifelse(is.na(biotype) | biotype == "", "Unknown", biotype),
      description = clean_annotation_text(ifelse(is.na(description), "", description)),
      chromosome = ifelse(is.na(chromosome), "", chromosome),
      gene_start = suppressWarnings(as.integer(gene_start)),
      gene_end = suppressWarnings(as.integer(gene_end)),
      strand = ifelse(is.na(strand), "", strand),
      location = ifelse(is.na(location), "", location),
      found_in_tpm = matched_gene_id %in% tpm_genes,
      found_in_expression_matrix = found_in_tpm,
      gene_display_label = make_gene_display_label(gene_name, matched_gene_id, description),
      query_display = ifelse(query == matched_gene_id, gene_display_label, paste(query, "->", gene_display_label))
    ) %>%
    dplyr::distinct(group, query, matched_gene_id, .keep_all = TRUE)
}

gene_display_lookup <- function(gene_catalog) {
  gene_catalog %>%
    dplyr::select(matched_gene_id, gene_name, gene_display_label, description, group) %>%
    dplyr::mutate(
      matched_gene_id = as.character(matched_gene_id),
      gene_name = ifelse(is.na(gene_name) | gene_name == "", matched_gene_id, as.character(gene_name)),
      description = clean_annotation_text(ifelse(is.na(description), "", description)),
      gene_display_label = ifelse(is.na(gene_display_label) | gene_display_label == "", make_gene_display_label(gene_name, matched_gene_id, description), gene_display_label),
      group = ifelse(is.na(group) | group == "", "unknown", as.character(group))
    ) %>%
    dplyr::distinct() %>%
    dplyr::group_by(matched_gene_id) %>%
    dplyr::summarise(
      gene_name = paste(sort(unique(gene_name)), collapse = "; "),
      gene_display_label = paste(sort(unique(gene_display_label)), collapse = "; "),
      description = paste(sort(unique(description[description != ""])), collapse = "; "),
      group = paste(sort(unique(group)), collapse = "; "),
      .groups = "drop"
    )
}

annotate_deg_hits <- function(deg_hits, gene_catalog) {
  lookup <- gene_display_lookup(gene_catalog)
  deg_hits %>%
    dplyr::select(-dplyr::any_of(c("group", "gene_name", "gene_display_label", "description"))) %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::relocate(dplyr::any_of(c("group", "gene_name", "gene_display_label", "description")), .after = gene_id)
}

left_join_gene_catalog <- function(x, y, by) {
  if ("relationship" %in% names(formals(dplyr::left_join))) {
    dplyr::left_join(x, y, by = by, relationship = "many-to-many")
  } else {
    suppressWarnings(dplyr::left_join(x, y, by = by))
  }
}

load_deg_hits <- function(deg_root, gene_catalog) {
  empty_deg <- tibble::tibble(
    gene_id = character(),
    contrast = character(),
    source_file = character(),
    result_dir = character(),
    deg_project = character(),
    deg_mode = character(),
    contrast_label = character(),
    padj_num = numeric(),
    log2FoldChange_num = numeric(),
    neg_log10_padj = numeric(),
    significant = logical()
  )
  files <- list.files(deg_root, pattern = "DEGs(_all)?_results\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(empty_deg)
  rows <- lapply(files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    if (!"contrast" %in% colnames(df)) df$contrast <- tools::file_path_sans_ext(basename(path))
    rel <- gsub("\\\\", "/", sub(paste0("^", normalizePath(deg_root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    result_dir <- dirname(rel)
    df %>%
      dplyr::filter(gene_id %in% gene_catalog$matched_gene_id) %>%
      dplyr::mutate(
        source_file = rel,
        result_dir = result_dir,
        deg_project = ifelse(result_dir %in% c(".", ""), "unknown", sub("/.*$", "", result_dir)),
        deg_mode = ifelse(grepl("/", result_dir), sub("^.*/", "", result_dir), "unknown"),
        contrast_label = paste(result_dir, contrast, sep = " | "),
        padj_num = suppressWarnings(as.numeric(padj)),
        log2FoldChange_num = suppressWarnings(as.numeric(log2FoldChange)),
        neg_log10_padj = ifelse(!is.na(padj_num) & padj_num > 0, -log10(padj_num), NA_real_),
        significant = !is.na(padj_num) & padj_num < 0.05 & abs(log2FoldChange_num) >= 1
      )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) empty_deg else out
}

annotate_dtu_hits <- function(dtu_hits, gene_catalog) {
  lookup <- gene_display_lookup(gene_catalog)
  dtu_hits %>%
    dplyr::select(-dplyr::any_of(c("group", "gene_name", "gene_display_label", "description"))) %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::relocate(dplyr::any_of(c("group", "gene_name", "gene_display_label", "description")), .after = gene_id)
}

load_dtu_hits <- function(dtu_root, gene_catalog = NULL) {
  empty_dtu <- tibble::tibble(
    gene_id = character(),
    transcript_id = character(),
    variable = character(),
    source_file = character(),
    dtu_project = character(),
    dtu_scope = character(),
    max_level = character(),
    min_level = character(),
    delta_usage = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    method = character()
  )
  if (dtu_root == "" || !dir.exists(dtu_root)) return(empty_dtu)
  files <- list.files(dtu_root, pattern = "dtu_significant\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(empty_dtu)
  rows <- lapply(files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    for (nm in c("transcript_id", "variable", "max_level", "min_level", "delta_usage", "pvalue", "padj", "method")) {
      if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
    }
    rel <- gsub("\\\\", "/", sub(paste0("^", normalizePath(dtu_root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    result_dir <- dirname(rel)
    dtu_project_value <- if ("project" %in% colnames(df)) df$project else ifelse(result_dir %in% c(".", ""), "unknown", sub("/.*$", "", result_dir))
    dtu_scope_value <- if ("scope" %in% colnames(df)) df$scope else ifelse(dtu_project_value == "all_projects", "all_projects", "project")
    # Keep all significant DTU hits. Gene-level sections filter against the
    # report catalog later, so the overview can still reflect the global DTU run.
    df %>%
      dplyr::mutate(
        source_file = rel,
        dtu_project = dtu_project_value,
        dtu_scope = dtu_scope_value,
        delta_usage = suppressWarnings(as.numeric(delta_usage)),
        pvalue = suppressWarnings(as.numeric(pvalue)),
        padj = suppressWarnings(as.numeric(padj))
      )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) empty_dtu else out
}

annotate_splicing_hits <- function(splicing_hits, gene_catalog) {
  lookup <- gene_display_lookup(gene_catalog)
  splicing_hits %>%
    dplyr::select(-dplyr::any_of(c("group", "gene_name", "gene_display_label", "description"))) %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::relocate(dplyr::any_of(c("group", "gene_name", "gene_display_label", "description")), .after = gene_id)
}

load_splicing_hits <- function(splicing_root, gene_catalog) {
  empty_splicing <- tibble::tibble(
    gene_id = character(),
    event_type = character(),
    source_file = character(),
    splicing_project = character(),
    splicing_variable = character(),
    splicing_contrast = character(),
    FDR = numeric(),
    IncLevelDifference = numeric()
  )
  if (splicing_root == "" || !dir.exists(splicing_root)) return(empty_splicing)
  files <- list.files(splicing_root, pattern = "significant_events\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(empty_splicing)

  catalog_ids <- gene_catalog$matched_gene_id
  catalog_names <- gene_catalog %>%
    dplyr::filter(!is.na(gene_name), gene_name != "") %>%
    dplyr::select(matched_gene_id, gene_name) %>%
    dplyr::distinct()

  normalize_gene_id <- function(x) {
    x <- as.character(x)
    x <- gsub("^gene:", "", x)
    x <- gsub("\\.[0-9]+$", "", x)
    x
  }

  rows <- lapply(files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    gene_col <- intersect(c("GeneID", "gene_id", "geneID"), colnames(df))[1]
    symbol_col <- intersect(c("geneSymbol", "gene_name", "symbol"), colnames(df))[1]
    if (is.na(gene_col) && is.na(symbol_col)) return(NULL)
    if (!"FDR" %in% colnames(df)) df$FDR <- NA_character_
    if (!"IncLevelDifference" %in% colnames(df)) df$IncLevelDifference <- NA_character_
    if (!"event_type" %in% colnames(df)) df$event_type <- "unknown"

    event_gene_id <- if (!is.na(gene_col)) normalize_gene_id(df[[gene_col]]) else rep("", nrow(df))
    event_gene_symbol <- if (!is.na(symbol_col)) as.character(df[[symbol_col]]) else rep("", nrow(df))
    by_id <- event_gene_id %in% catalog_ids
    by_symbol_match <- match(event_gene_symbol, catalog_names$gene_name)
    by_symbol <- !is.na(by_symbol_match)
    matched_gene <- ifelse(by_id, event_gene_id, ifelse(by_symbol, catalog_names$matched_gene_id[by_symbol_match], NA_character_))

    rel <- gsub("\\\\", "/", sub(paste0("^", normalizePath(splicing_root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    parts <- strsplit(dirname(rel), "/", fixed = TRUE)[[1]]
    splicing_project <- ifelse(length(parts) >= 1, parts[1], "unknown")
    splicing_variable <- ifelse(length(parts) >= 2, parts[2], "unknown")
    splicing_contrast <- ifelse(length(parts) >= 3, parts[3], "unknown")

    df %>%
      dplyr::mutate(
        gene_id = matched_gene,
        source_file = rel,
        splicing_project = splicing_project,
        splicing_variable = splicing_variable,
        splicing_contrast = splicing_contrast,
        FDR = suppressWarnings(as.numeric(FDR)),
        IncLevelDifference = suppressWarnings(as.numeric(IncLevelDifference))
      ) %>%
      dplyr::filter(!is.na(gene_id), gene_id != "")
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) empty_splicing else out
}

annotate_wgcna_hits <- function(wgcna_hits, gene_catalog) {
  lookup <- gene_display_lookup(gene_catalog)
  wgcna_hits %>%
    dplyr::select(-dplyr::any_of(c("group", "gene_name", "gene_display_label", "description"))) %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::relocate(dplyr::any_of(c("group", "gene_name", "gene_display_label", "description")), .after = gene_id)
}

load_wgcna_hits <- function(wgcna_root, gene_catalog) {
  empty_wgcna <- tibble::tibble(
    gene_id = character(),
    source_file = character(),
    wgcna_project = character(),
    wgcna_scope = character(),
    module_color = character(),
    module_label = integer(),
    is_hub = logical(),
    kME = numeric(),
    abs_kME = numeric(),
    module_rank = integer()
  )
  if (wgcna_root == "" || !dir.exists(wgcna_root)) return(empty_wgcna)

  root_norm <- normalizePath(wgcna_root, winslash = "/", mustWork = FALSE)
  project_from_path <- function(path) {
    rel <- gsub("\\\\", "/", sub(paste0("^", root_norm, "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    parts <- strsplit(dirname(rel), "/", fixed = TRUE)[[1]]
    ifelse(length(parts) >= 1 && parts[1] != ".", parts[1], "unknown")
  }
  rel_from_path <- function(path) {
    gsub("\\\\", "/", sub(paste0("^", root_norm, "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
  }

  assignment_files <- list.files(wgcna_root, pattern = "module_assignments\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  assignment_rows <- lapply(assignment_files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    if (!"module_color" %in% colnames(df)) df$module_color <- NA_character_
    if (!"module_label" %in% colnames(df)) df$module_label <- NA_character_
    project <- project_from_path(path)
    df %>%
      dplyr::filter(gene_id %in% gene_catalog$matched_gene_id) %>%
      dplyr::transmute(
        gene_id = gene_id,
        source_file = rel_from_path(path),
        wgcna_project = project,
        wgcna_scope = ifelse(project == "all_projects", "all_projects", "project"),
        module_color = module_color,
        module_label = suppressWarnings(as.integer(module_label))
      )
  })
  assignments <- dplyr::bind_rows(assignment_rows)

  hub_files <- list.files(wgcna_root, pattern = "hub_genes\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  hub_rows <- lapply(hub_files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    for (nm in c("module_color", "kME", "abs_kME", "module_rank")) {
      if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
    }
    project <- project_from_path(path)
    df %>%
      dplyr::filter(gene_id %in% gene_catalog$matched_gene_id) %>%
      dplyr::transmute(
        gene_id = gene_id,
        wgcna_project = project,
        module_color = module_color,
        kME = suppressWarnings(as.numeric(kME)),
        abs_kME = suppressWarnings(as.numeric(abs_kME)),
        module_rank = suppressWarnings(as.integer(module_rank))
      )
  })
  hubs <- dplyr::bind_rows(hub_rows)

  if (nrow(assignments) == 0 && nrow(hubs) == 0) return(empty_wgcna)
  if (nrow(assignments) == 0) {
    out <- hubs %>%
      dplyr::mutate(
        source_file = "",
        wgcna_scope = ifelse(wgcna_project == "all_projects", "all_projects", "project"),
        module_label = NA_integer_,
        is_hub = TRUE
      ) %>%
      dplyr::select(gene_id, source_file, wgcna_project, wgcna_scope, module_color, module_label, is_hub, kME, abs_kME, module_rank)
    return(out)
  }

  out <- assignments %>%
    dplyr::left_join(hubs, by = c("gene_id", "wgcna_project", "module_color")) %>%
    dplyr::mutate(
      is_hub = !is.na(module_rank),
      kME = suppressWarnings(as.numeric(kME)),
      abs_kME = suppressWarnings(as.numeric(abs_kME)),
      module_rank = suppressWarnings(as.integer(module_rank))
    ) %>%
    dplyr::arrange(wgcna_project, module_color, module_rank, gene_id)
  if (nrow(out) == 0) empty_wgcna else out
}

annotate_mfuzz_hits <- function(mfuzz_hits, gene_catalog) {
  lookup <- gene_display_lookup(gene_catalog)
  mfuzz_hits %>%
    dplyr::select(-dplyr::any_of(c("group", "gene_name", "gene_display_label", "description"))) %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::relocate(dplyr::any_of(c("group", "gene_name", "gene_display_label", "description")), .after = gene_id)
}

load_mfuzz_hits <- function(mfuzz_root, gene_catalog) {
  empty_mfuzz <- tibble::tibble(
    gene_id = character(),
    source_file = character(),
    mfuzz_project = character(),
    mfuzz_scope = character(),
    cluster = integer(),
    membership = numeric()
  )
  if (mfuzz_root == "" || !dir.exists(mfuzz_root)) return(empty_mfuzz)
  root_norm <- normalizePath(mfuzz_root, winslash = "/", mustWork = FALSE)
  files <- list.files(mfuzz_root, pattern = "mfuzz_membership\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(empty_mfuzz)
  rows <- lapply(files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    if (!"cluster" %in% colnames(df)) df$cluster <- NA_character_
    if (!"membership" %in% colnames(df)) df$membership <- NA_character_
    rel <- gsub("\\\\", "/", sub(paste0("^", root_norm, "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    parts <- strsplit(dirname(rel), "/", fixed = TRUE)[[1]]
    project <- ifelse(length(parts) >= 1 && parts[1] != ".", parts[1], "unknown")
    df %>%
      dplyr::filter(gene_id %in% gene_catalog$matched_gene_id) %>%
      dplyr::mutate(
        source_file = rel,
        mfuzz_project = project,
        mfuzz_scope = ifelse(project == "all_projects", "all_projects", "project"),
        cluster = suppressWarnings(as.integer(cluster)),
        membership = suppressWarnings(as.numeric(membership))
      )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) empty_mfuzz else out
}

complete_sample_fields <- function(samples) {
  for (nm in c("dataset", "sample_id", "stage", "tissue", "sex", "condition", "batch")) {
    if (!nm %in% colnames(samples)) samples[[nm]] <- NA_character_
  }
  out <- samples %>%
    dplyr::mutate(
      dataset = ifelse(is.na(dataset) | dataset == "", "unknown", dataset),
      sample_id = ifelse(is.na(sample_id) | sample_id == "", import_id, sample_id),
      stage_raw = ifelse(is.na(stage) | stage == "", "unknown", as.character(stage)),
      stage = normalize_stage_detail(stage_raw),
      stage_class = classify_life_stage(stage),
      stage_class = ifelse(stage_class %in% life_stage_levels, stage_class, "unknown"),
      stage_class = factor(stage_class, levels = life_stage_levels),
      stage_day = extract_stage_day(stage),
      tissue = ifelse(is.na(tissue) | tissue == "", "unknown", tissue),
      sex = ifelse(is.na(sex) | sex == "", "unknown", sex),
      condition = ifelse(is.na(condition) | condition == "", "unknown", condition),
      batch = ifelse(is.na(batch) | batch == "", dataset, batch)
    )
  out$stage <- factor(out$stage, levels = order_stage_details(out$stage))
  out
}

make_expression_long <- function(tpm, samples, gene_catalog, expression_source = "raw") {
  selected <- tpm %>% dplyr::filter(gene_id %in% gene_catalog$matched_gene_id)
  selected %>%
    tidyr::pivot_longer(-gene_id, names_to = "import_id", values_to = "TPM") %>%
    dplyr::left_join(samples, by = "import_id") %>%
    left_join_gene_catalog(gene_catalog %>% dplyr::select(group, query, matched_gene_id, gene_name, gene_display_label, description, biotype, chromosome, gene_start, gene_end, strand, location), by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::mutate(
      TPM = as.numeric(TPM),
      log2TPM = log2(TPM + 1),
      gene_display_label = ifelse(is.na(gene_display_label) | gene_display_label == "", make_gene_display_label(gene_name, gene_id, description), gene_display_label),
      sample_label = paste(dataset, sample_id, sep = " | "),
      context_full = paste(dataset, batch, condition, stage, tissue, sex, sep = " | "),
      context_biology = paste(condition, stage, tissue, sex, sep = " | "),
      expression_source = expression_source
    ) %>%
    dplyr::group_by(group, gene_id) %>%
    dplyr::mutate(z_log2TPM = as.numeric(scale(log2TPM))) %>%
    dplyr::ungroup()
}

summarise_expression <- function(expr_long) {
  expr_long %>%
    dplyr::group_by(group, gene_id, gene_name, gene_display_label, dataset, batch, condition, stage_class, stage_day, stage, tissue, sex) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_TPM = mean(TPM, na.rm = TRUE),
      median_TPM = median(TPM, na.rm = TRUE),
      mean_log2TPM = mean(log2TPM, na.rm = TRUE),
      fraction_expressed = mean(TPM > 1, na.rm = TRUE),
      .groups = "drop"
    )
}

stage_specificity_from_stage_means <- function(df) {
  df <- df %>%
    dplyr::arrange(stage_class, stage_day, stage) %>%
    dplyr::filter(!is.na(stage), as.character(stage) != "", as.character(stage) != "unknown")

  if (nrow(df) == 0) {
    return(tibble::tibble(
      n_stages_tested = 0L,
      dominant_stage = NA_character_,
      max_stage_mean_expression = NA_real_,
      second_stage_mean_expression = NA_real_,
      stage_specificity_tau = NA_real_,
      stage_specificity_max_fraction = NA_real_,
      stage_specificity_fold_vs_second = NA_real_,
      expression_stage_specific = FALSE
    ))
  }

  values <- as.numeric(df$mean_expression)
  stages <- as.character(df$stage)
  max_idx <- which.max(values)[1]
  max_value <- values[[max_idx]]
  sorted_values <- sort(values, decreasing = TRUE)
  second_value <- if (length(sorted_values) >= 2) sorted_values[[2]] else NA_real_
  tau <- if (length(values) > 1 && is.finite(max_value) && max_value > 0) {
    sum(1 - (values / max_value), na.rm = TRUE) / (length(values) - 1)
  } else {
    NA_real_
  }
  max_fraction <- if (sum(values, na.rm = TRUE) > 0) max_value / sum(values, na.rm = TRUE) else NA_real_
  fold_vs_second <- if (!is.na(second_value) && second_value > 0) max_value / second_value else NA_real_

  tibble::tibble(
    n_stages_tested = length(values),
    dominant_stage = stages[[max_idx]],
    max_stage_mean_expression = max_value,
    second_stage_mean_expression = second_value,
    stage_specificity_tau = tau,
    stage_specificity_max_fraction = max_fraction,
    stage_specificity_fold_vs_second = fold_vs_second,
    expression_stage_specific = !is.na(tau) && tau >= stage_tau_threshold && max_value >= stage_min_expression
  )
}

summarise_stage_specificity <- function(expr_long, deg_hits, mfuzz_hits) {
  empty <- tibble::tibble(
    group = character(),
    gene_id = character(),
    gene_name = character(),
    gene_display_label = character(),
    n_stages_tested = integer(),
    dominant_stage = character(),
    max_stage_mean_expression = numeric(),
    second_stage_mean_expression = numeric(),
    stage_specificity_tau = numeric(),
    stage_specificity_max_fraction = numeric(),
    stage_specificity_fold_vs_second = numeric(),
    expression_stage_specific = logical(),
    n_stage_deg_significant = integer(),
    stage_deg_significant = logical(),
    min_stage_deg_padj = numeric(),
    max_abs_stage_deg_log2FC = numeric(),
    top_stage_deg_contrast = character(),
    max_mfuzz_membership = numeric(),
    mfuzz_stage_specific = logical(),
    top_mfuzz_cluster = character(),
    stage_specificity_candidate = logical(),
    stage_specificity_evidence = character()
  )

  expr_stage <- expr_long %>%
    dplyr::filter(!is.na(stage), as.character(stage) != "", as.character(stage) != "unknown") %>%
    dplyr::group_by(group, gene_id, gene_name, gene_display_label, stage_class, stage_day, stage) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      mean_expression = mean(TPM, na.rm = TRUE),
      mean_log2_expression = mean(log2TPM, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(expr_stage) == 0) return(empty)

  expression_specificity <- expr_stage %>%
    dplyr::group_by(group, gene_id, gene_name, gene_display_label) %>%
    dplyr::group_modify(~ stage_specificity_from_stage_means(.x)) %>%
    dplyr::ungroup()

  if (!"variable" %in% colnames(deg_hits)) deg_hits$variable <- ""
  stage_deg <- deg_hits %>%
    dplyr::mutate(
      variable = ifelse(is.na(variable), "", as.character(variable)),
      contrast = ifelse(is.na(contrast), "", as.character(contrast)),
      is_stage_contrast = variable == "stage" | grepl("^stage__", contrast) | grepl("[| ][ ]*stage__", contrast)
    ) %>%
    dplyr::filter(is_stage_contrast, significant) %>%
    dplyr::arrange(padj_num) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      n_stage_deg_significant = dplyr::n(),
      stage_deg_significant = TRUE,
      min_stage_deg_padj = suppressWarnings(min(padj_num, na.rm = TRUE)),
      max_abs_stage_deg_log2FC = suppressWarnings(max(abs(log2FoldChange_num), na.rm = TRUE)),
      top_stage_deg_contrast = dplyr::first(contrast),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      min_stage_deg_padj = ifelse(is.infinite(min_stage_deg_padj), NA_real_, min_stage_deg_padj),
      max_abs_stage_deg_log2FC = ifelse(is.infinite(max_abs_stage_deg_log2FC), NA_real_, max_abs_stage_deg_log2FC)
    )

  mfuzz_summary <- if (nrow(mfuzz_hits) > 0) {
    mfuzz_hits %>%
      dplyr::mutate(membership_num = suppressWarnings(as.numeric(membership))) %>%
      dplyr::arrange(dplyr::desc(membership_num)) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        max_mfuzz_membership = suppressWarnings(max(membership_num, na.rm = TRUE)),
        top_mfuzz_cluster = as.character(dplyr::first(cluster)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        max_mfuzz_membership = ifelse(is.infinite(max_mfuzz_membership), NA_real_, max_mfuzz_membership),
        mfuzz_stage_specific = !is.na(max_mfuzz_membership) & max_mfuzz_membership >= mfuzz_membership_threshold
      )
  } else {
    tibble::tibble(gene_id = character(), max_mfuzz_membership = numeric(), mfuzz_stage_specific = logical(), top_mfuzz_cluster = character())
  }

  expression_specificity %>%
    dplyr::left_join(stage_deg, by = "gene_id") %>%
    dplyr::left_join(mfuzz_summary, by = "gene_id") %>%
    dplyr::mutate(
      n_stage_deg_significant = ifelse(is.na(n_stage_deg_significant), 0L, n_stage_deg_significant),
      stage_deg_significant = ifelse(is.na(stage_deg_significant), FALSE, stage_deg_significant),
      mfuzz_stage_specific = ifelse(is.na(mfuzz_stage_specific), FALSE, mfuzz_stage_specific),
      stage_specificity_candidate = expression_stage_specific | stage_deg_significant | mfuzz_stage_specific
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      stage_specificity_evidence = paste(
        c(
          if (expression_stage_specific) "expression_tau" else character(),
          if (stage_deg_significant) "stage_deg" else character(),
          if (mfuzz_stage_specific) "mfuzz_membership" else character()
        ),
        collapse = "; "
      ),
      stage_specificity_evidence = ifelse(stage_specificity_evidence == "", "none", stage_specificity_evidence)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group, dplyr::desc(stage_specificity_candidate), dominant_stage, dplyr::desc(stage_specificity_tau), gene_display_label)
}

add_stage_specificity_summary <- function(gene_summary, stage_specificity) {
  cols <- c(
    "group", "gene_id", "dominant_stage", "n_stages_tested", "stage_specificity_tau",
    "stage_specificity_max_fraction", "stage_specificity_fold_vs_second",
    "expression_stage_specific", "n_stage_deg_significant", "stage_deg_significant",
    "min_stage_deg_padj", "max_abs_stage_deg_log2FC", "max_mfuzz_membership",
    "mfuzz_stage_specific", "top_mfuzz_cluster", "stage_specificity_candidate",
    "stage_specificity_evidence"
  )
  if (is.null(stage_specificity) || nrow(stage_specificity) == 0) return(gene_summary)
  gene_summary %>%
    dplyr::left_join(stage_specificity %>% dplyr::select(dplyr::any_of(cols)), by = c("group", "gene_id")) %>%
    dplyr::mutate(
      n_stages_tested = ifelse(is.na(n_stages_tested), 0L, n_stages_tested),
      expression_stage_specific = ifelse(is.na(expression_stage_specific), FALSE, expression_stage_specific),
      n_stage_deg_significant = ifelse(is.na(n_stage_deg_significant), 0L, n_stage_deg_significant),
      stage_deg_significant = ifelse(is.na(stage_deg_significant), FALSE, stage_deg_significant),
      mfuzz_stage_specific = ifelse(is.na(mfuzz_stage_specific), FALSE, mfuzz_stage_specific),
      stage_specificity_candidate = ifelse(is.na(stage_specificity_candidate), FALSE, stage_specificity_candidate),
      stage_specificity_evidence = ifelse(is.na(stage_specificity_evidence) | stage_specificity_evidence == "", "none", stage_specificity_evidence)
    )
}

summarise_gene_descriptives <- function(expr_long, expr_summary, deg_hits, gene_catalog) {
  expr_gene <- expr_long %>%
    dplyr::group_by(group, gene_id, gene_name, gene_display_label) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      mean_TPM = mean(TPM, na.rm = TRUE),
      median_TPM = median(TPM, na.rm = TRUE),
      max_TPM = max(TPM, na.rm = TRUE),
      fraction_samples_TPM_gt1 = mean(TPM > 1, na.rm = TRUE),
      n_datasets = dplyr::n_distinct(dataset),
      datasets = paste(sort(unique(dataset)), collapse = "; "),
      n_batches = dplyr::n_distinct(batch),
      n_tissues = dplyr::n_distinct(tissue),
      .groups = "drop"
    )
  dominant <- expr_summary %>%
    dplyr::mutate(context = paste(dataset, batch, condition, stage, tissue, sex, sep = " | ")) %>%
    dplyr::group_by(group, gene_id) %>%
    dplyr::arrange(dplyr::desc(mean_log2TPM), .by_group = TRUE) %>%
    dplyr::summarise(context_with_highest_expression = dplyr::first(context), .groups = "drop")
  deg_summary <- if (nrow(deg_hits) > 0) {
    deg_hits %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        n_deg_records = dplyr::n(),
        n_significant_contrasts = sum(significant, na.rm = TRUE),
        n_deg_projects = dplyr::n_distinct(deg_project),
        deg_projects = paste(sort(unique(deg_project)), collapse = "; "),
        n_deg_modes = dplyr::n_distinct(deg_mode),
        max_abs_log2FC = suppressWarnings(max(abs(log2FoldChange_num), na.rm = TRUE)),
        min_padj = suppressWarnings(min(padj_num, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        max_abs_log2FC = ifelse(is.infinite(max_abs_log2FC), NA_real_, max_abs_log2FC),
        min_padj = ifelse(is.infinite(min_padj), NA_real_, min_padj)
      )
  } else {
    tibble::tibble(
      gene_id = character(),
      n_deg_records = integer(),
      n_significant_contrasts = integer(),
      n_deg_projects = integer(),
      deg_projects = character(),
      n_deg_modes = integer(),
      max_abs_log2FC = numeric(),
      min_padj = numeric()
    )
  }
  gene_catalog %>%
    dplyr::select(group, query, query_display, matched_gene_id, gene_name, gene_display_label, biotype, description, chromosome, gene_start, gene_end, strand, location, found_in_tpm, found_in_expression_matrix) %>%
    dplyr::rename(gene_id = matched_gene_id) %>%
    dplyr::left_join(expr_gene, by = c("group", "gene_id", "gene_name", "gene_display_label")) %>%
    dplyr::left_join(dominant, by = c("group", "gene_id")) %>%
    dplyr::left_join(deg_summary, by = "gene_id") %>%
    dplyr::mutate(
      dplyr::across(c(n_deg_records, n_significant_contrasts, n_deg_projects, n_deg_modes), ~ ifelse(is.na(.x), 0, .x))
    ) %>%
    dplyr::mutate(
      datasets = ifelse(is.na(datasets), "", datasets),
      deg_projects = ifelse(is.na(deg_projects), "", deg_projects)
    ) %>%
    dplyr::arrange(group, gene_name, gene_id)
}

add_optional_gene_summaries <- function(gene_summary, dtu_hits, splicing_hits, wgcna_hits, mfuzz_hits) {
  dtu_summary <- if (nrow(dtu_hits) > 0) {
    dtu_hits %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        n_dtu_records = dplyr::n(),
        n_dtu_transcripts = dplyr::n_distinct(transcript_id),
        dtu_projects = paste(sort(unique(dtu_project)), collapse = "; "),
        min_dtu_padj = suppressWarnings(min(padj, na.rm = TRUE)),
        max_delta_usage = suppressWarnings(max(delta_usage, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        min_dtu_padj = ifelse(is.infinite(min_dtu_padj), NA_real_, min_dtu_padj),
        max_delta_usage = ifelse(is.infinite(max_delta_usage), NA_real_, max_delta_usage)
      )
  } else {
    tibble::tibble(
      gene_id = character(),
      n_dtu_records = integer(),
      n_dtu_transcripts = integer(),
      dtu_projects = character(),
      min_dtu_padj = numeric(),
      max_delta_usage = numeric()
    )
  }

  splicing_summary <- if (nrow(splicing_hits) > 0) {
    splicing_hits %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        n_splicing_records = dplyr::n(),
        n_splicing_event_types = dplyr::n_distinct(event_type),
        splicing_projects = paste(sort(unique(splicing_project)), collapse = "; "),
        min_splicing_fdr = suppressWarnings(min(FDR, na.rm = TRUE)),
        max_abs_inc_level_difference = suppressWarnings(max(abs(IncLevelDifference), na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        min_splicing_fdr = ifelse(is.infinite(min_splicing_fdr), NA_real_, min_splicing_fdr),
        max_abs_inc_level_difference = ifelse(is.infinite(max_abs_inc_level_difference), NA_real_, max_abs_inc_level_difference)
      )
  } else {
    tibble::tibble(
      gene_id = character(),
      n_splicing_records = integer(),
      n_splicing_event_types = integer(),
      splicing_projects = character(),
      min_splicing_fdr = numeric(),
      max_abs_inc_level_difference = numeric()
    )
  }

  wgcna_summary <- if (nrow(wgcna_hits) > 0) {
    wgcna_hits %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        n_wgcna_records = dplyr::n(),
        n_wgcna_modules = dplyr::n_distinct(paste(wgcna_project, module_color, sep = "::")),
        n_wgcna_hub_records = sum(is_hub, na.rm = TRUE),
        wgcna_projects = paste(sort(unique(wgcna_project)), collapse = "; "),
        wgcna_modules = paste(sort(unique(module_color[!is.na(module_color) & module_color != ""])), collapse = "; "),
        max_abs_kME = suppressWarnings(max(abs_kME, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(max_abs_kME = ifelse(is.infinite(max_abs_kME), NA_real_, max_abs_kME))
  } else {
    tibble::tibble(
      gene_id = character(),
      n_wgcna_records = integer(),
      n_wgcna_modules = integer(),
      n_wgcna_hub_records = integer(),
      wgcna_projects = character(),
      wgcna_modules = character(),
      max_abs_kME = numeric()
    )
  }

  mfuzz_summary <- if (nrow(mfuzz_hits) > 0) {
    mfuzz_hits %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        n_mfuzz_records = dplyr::n(),
        n_mfuzz_clusters = dplyr::n_distinct(paste(mfuzz_project, cluster, sep = "::")),
        mfuzz_projects = paste(sort(unique(mfuzz_project)), collapse = "; "),
        mfuzz_clusters = paste(sort(unique(cluster[!is.na(cluster)])), collapse = "; "),
        max_mfuzz_membership = suppressWarnings(max(membership, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(max_mfuzz_membership = ifelse(is.infinite(max_mfuzz_membership), NA_real_, max_mfuzz_membership))
  } else {
    tibble::tibble(
      gene_id = character(),
      n_mfuzz_records = integer(),
      n_mfuzz_clusters = integer(),
      mfuzz_projects = character(),
      mfuzz_clusters = character(),
      max_mfuzz_membership = numeric()
    )
  }

  gene_summary %>%
    dplyr::left_join(dtu_summary, by = "gene_id") %>%
    dplyr::left_join(splicing_summary, by = "gene_id") %>%
    dplyr::left_join(wgcna_summary, by = "gene_id") %>%
    dplyr::left_join(mfuzz_summary, by = "gene_id") %>%
    dplyr::mutate(
      dplyr::across(c(n_dtu_records, n_dtu_transcripts, n_splicing_records, n_splicing_event_types, n_wgcna_records, n_wgcna_modules, n_wgcna_hub_records, n_mfuzz_records, n_mfuzz_clusters), ~ ifelse(is.na(.x), 0, .x)),
      dtu_projects = ifelse(is.na(dtu_projects), "", dtu_projects),
      splicing_projects = ifelse(is.na(splicing_projects), "", splicing_projects),
      wgcna_projects = ifelse(is.na(wgcna_projects), "", wgcna_projects),
      wgcna_modules = ifelse(is.na(wgcna_modules), "", wgcna_modules),
      mfuzz_projects = ifelse(is.na(mfuzz_projects), "", mfuzz_projects),
      mfuzz_clusters = ifelse(is.na(mfuzz_clusters), "", mfuzz_clusters)
    )
}

heatmap_scale_mode <- function(mat) {
  if (nrow(mat) > 1 && ncol(mat) > 1) "row" else "none"
}

heatmap_has_signal <- function(mat) {
  values <- as.numeric(mat)
  values <- values[is.finite(values)]
  length(unique(values)) > 1
}

plot_expression_heatmap <- function(expr_summary, outfile, title = "Expressao media por contexto") {
  mat_df <- expr_summary %>%
    dplyr::arrange(dataset, batch, condition, stage_class, stage_day, stage, tissue, sex, group, gene_display_label) %>%
    dplyr::mutate(
      label = paste(group, gene_display_label, sep = " | "),
      context = paste(dataset, batch, condition, stage, tissue, sex, sep = " | ")
    ) %>%
    dplyr::group_by(label, context) %>%
    dplyr::summarise(mean_log2TPM = mean(mean_log2TPM, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = context, values_from = mean_log2TPM, values_fill = 0)
  if (nrow(mat_df) == 0 || ncol(mat_df) < 2) return(FALSE)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$label
  if (!heatmap_has_signal(mat)) return(FALSE)
  pheatmap::pheatmap(mat, scale = heatmap_scale_mode(mat), border_color = NA,
                     cluster_rows = nrow(mat) > 1,
                     cluster_cols = ncol(mat) > 1,
                     fontsize_row = 7, fontsize_col = 6,
                     main = title, filename = outfile,
                     width = 14, height = max(5, min(18, nrow(mat) * 0.32 + 3)))
  TRUE
}

plot_expression_dotplot <- function(expr_summary, outfile, title = "Expressao media e fracao expressa") {
  df <- expr_summary %>%
    dplyr::arrange(dataset, batch, condition, stage_class, stage_day, stage, tissue, sex, group, gene_display_label) %>%
    dplyr::mutate(
      context = paste(dataset, batch, condition, stage, tissue, sex, sep = " | "),
      gene_label = paste(group, gene_display_label, sep = " | ")
    )
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = context, y = gene_label)) +
    geom_point(aes(size = fraction_expressed, color = mean_log2TPM), alpha = 0.85) +
    scale_color_viridis_c(option = "C") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1), panel.grid.major.y = element_line(color = "gray92")) +
    labs(title = title, x = "Projeto | batch | condicao | estagio | tecido | sexo", y = "Grupo | gene | ID",
         color = expression_mean_log_label, size = paste0("Frac. ", expression_unit, ">1"))
  ggsave(outfile, p, width = 15, height = max(5, min(18, length(unique(df$gene_label)) * 0.32 + 3)), dpi = 300)
  TRUE
}

plot_tissue_sex_heatmap <- function(expr_summary, outfile, title = "Padroes por tecido e sexo") {
  mat_df <- expr_summary %>%
    dplyr::arrange(tissue, sex, stage_class, stage_day, stage, group, gene_display_label) %>%
    dplyr::mutate(context = paste(tissue, sex, sep = " | "),
                  label = paste(group, gene_display_label, sep = " | ")) %>%
    dplyr::group_by(label, context) %>%
    dplyr::summarise(mean_log2TPM = mean(mean_log2TPM, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = context, values_from = mean_log2TPM, values_fill = 0)
  if (nrow(mat_df) == 0 || ncol(mat_df) < 2) return(FALSE)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$label
  if (!heatmap_has_signal(mat)) return(FALSE)
  pheatmap::pheatmap(mat, scale = heatmap_scale_mode(mat), border_color = NA,
                     cluster_rows = nrow(mat) > 1,
                     cluster_cols = ncol(mat) > 1,
                     fontsize_row = 7, fontsize_col = 8,
                     main = title, filename = outfile,
                     width = 11, height = max(5, min(18, nrow(mat) * 0.32 + 3)))
  TRUE
}

plot_batch_project_boxplot <- function(expr_long, outfile, title = "Expressao por projeto e batch") {
  if (nrow(expr_long) == 0) return(FALSE)
  df <- expr_long %>% dplyr::mutate(gene_label = paste(group, gene_display_label, sep = " | "))
  p <- ggplot(df, aes(x = batch, y = log2TPM, fill = dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = dataset), width = 0.18, alpha = 0.35, size = 1) +
    facet_wrap(~ gene_label, scales = "free_y") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "Batch", y = expression_log_label, fill = "Projeto", color = "Projeto")
  ggsave(outfile, p, width = 14, height = max(6, min(18, length(unique(df$gene_label)) * 1.15 + 3)), dpi = 300)
  TRUE
}

plot_group_sample_heatmap <- function(expr_long, outfile, title = "Amostras individuais") {
  mat_df <- expr_long %>%
    dplyr::mutate(
      gene_label = gene_display_label,
      sample_label = paste(dataset, batch, condition, stage, tissue, sex, sample_id, sep = " | ")
    ) %>%
    dplyr::group_by(gene_label, sample_label) %>%
    dplyr::summarise(log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_label, values_from = log2TPM, values_fill = 0)
  if (nrow(mat_df) == 0 || ncol(mat_df) < 2) return(FALSE)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$gene_label
  if (!heatmap_has_signal(mat)) return(FALSE)
  pheatmap::pheatmap(mat, scale = heatmap_scale_mode(mat), border_color = NA,
                     cluster_rows = nrow(mat) > 1,
                     cluster_cols = ncol(mat) > 1,
                     fontsize_row = 7, fontsize_col = 5,
                     main = title, filename = outfile,
                     width = 16, height = max(5, min(16, nrow(mat) * 0.35 + 3)))
  TRUE
}

expression_matrix_by_gene <- function(expr_long, include_group = TRUE) {
  mat_df <- expr_long %>%
    dplyr::mutate(gene_label = if (include_group) paste(group, gene_display_label, sep = " | ") else gene_display_label) %>%
    dplyr::group_by(gene_label, import_id) %>%
    dplyr::summarise(log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = import_id, values_from = log2TPM, values_fill = 0)
  if (nrow(mat_df) == 0 || ncol(mat_df) < 2) return(NULL)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$gene_label
  mat
}

sample_annotation_for_matrix <- function(expr_long, sample_ids) {
  ann <- expr_long %>%
    dplyr::distinct(import_id, dataset, batch, condition, stage, tissue, sex) %>%
    dplyr::filter(import_id %in% sample_ids)
  ann <- ann[match(sample_ids, ann$import_id), , drop = FALSE]
  ann <- as.data.frame(ann[, c("dataset", "batch", "condition", "stage", "tissue", "sex"), drop = FALSE])
  rownames(ann) <- sample_ids
  ann
}

plot_annotated_sample_heatmap <- function(expr_long, outfile, title = "Heatmap gene x amostra anotado") {
  mat <- expression_matrix_by_gene(expr_long, include_group = TRUE)
  if (is.null(mat) || nrow(mat) < 1 || ncol(mat) < 2) return(FALSE)
  if (!heatmap_has_signal(mat)) return(FALSE)
  ann_col <- sample_annotation_for_matrix(expr_long, colnames(mat))
  pheatmap::pheatmap(mat, scale = heatmap_scale_mode(mat), border_color = NA,
                     cluster_rows = nrow(mat) > 1,
                     cluster_cols = ncol(mat) > 1,
                     annotation_col = ann_col,
                     show_colnames = FALSE,
                     fontsize_row = 7,
                     main = title,
                     filename = outfile,
                     width = 15, height = max(5, min(18, nrow(mat) * 0.32 + 4)))
  TRUE
}

plot_gene_correlation <- function(expr_long, outfile, title = "Correlacao entre genes") {
  mat <- expression_matrix_by_gene(expr_long, include_group = TRUE)
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 3) return(FALSE)
  keep <- apply(mat, 1, stats::sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) return(FALSE)
  cor_mat <- stats::cor(t(mat), use = "pairwise.complete.obs", method = "spearman")
  pheatmap::pheatmap(cor_mat, border_color = NA,
                     color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
                     breaks = seq(-1, 1, length.out = 102),
                     fontsize_row = 7, fontsize_col = 7,
                     main = title,
                     filename = outfile,
                     width = max(6, min(16, nrow(cor_mat) * 0.28 + 4)),
                     height = max(6, min(16, nrow(cor_mat) * 0.28 + 4)))
  TRUE
}

sample_scores_long <- function(expr_long, method = c("pca", "mds")) {
  method <- match.arg(method)
  mat <- expression_matrix_by_gene(expr_long, include_group = TRUE)
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 3) return(tibble::tibble())
  keep <- apply(mat, 1, stats::sd, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) return(tibble::tibble())
  sample_mat <- t(mat)
  if (method == "pca") {
    pc <- stats::prcomp(sample_mat, center = TRUE, scale. = TRUE)
    coords <- as.data.frame(pc$x[, 1:2, drop = FALSE])
    names(coords) <- c("Dim1", "Dim2")
    variance <- round(100 * (pc$sdev^2 / sum(pc$sdev^2))[1:2], 1)
    axis_labels <- c(paste0("PC1 (", variance[1], "%)"), paste0("PC2 (", variance[2], "%)"))
  } else {
    d <- stats::dist(sample_mat)
    coords <- as.data.frame(stats::cmdscale(d, k = 2))
    names(coords) <- c("Dim1", "Dim2")
    axis_labels <- c("MDS1", "MDS2")
  }
  coords$import_id <- rownames(sample_mat)
  ann <- expr_long %>% dplyr::distinct(import_id, dataset, batch, condition, stage, tissue, sex)
  coords <- coords %>% dplyr::left_join(ann, by = "import_id")
  vars <- c("dataset", "batch", "condition", "stage", "tissue", "sex")
  out <- dplyr::bind_rows(lapply(vars, function(v) {
    coords %>%
      dplyr::mutate(variable = v, value = as.character(.data[[v]])) %>%
      dplyr::select(import_id, Dim1, Dim2, variable, value)
  }))
  attr(out, "axis_labels") <- axis_labels
  out
}

plot_sample_ordination <- function(expr_long, outfile, method = c("pca", "mds"), title = "Ordenacao de amostras") {
  method <- match.arg(method)
  df <- sample_scores_long(expr_long, method = method)
  if (nrow(df) == 0) return(FALSE)
  axis_labels <- attr(df, "axis_labels")
  p <- ggplot(df, aes(x = Dim1, y = Dim2, color = value)) +
    geom_point(size = 2.2, alpha = 0.9) +
    facet_wrap(~ variable, scales = "free") +
    theme_bw(base_size = 10) +
    labs(title = title, x = axis_labels[1], y = axis_labels[2], color = "Valor")
  ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  TRUE
}

plot_ovary_testis_panel <- function(expr_summary, outfile, title = "Ovario versus testiculo") {
  if (!organism_specific_reports) return(FALSE)
  df <- expr_summary %>%
    dplyr::filter(tissue %in% c("ovary", "testis")) %>%
    dplyr::group_by(group, gene_display_label, dataset, stage, sex, tissue) %>%
    dplyr::summarise(mean_log2TPM = mean(mean_log2TPM, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = tissue, values_from = mean_log2TPM)
  if (nrow(df) == 0 || !all(c("ovary", "testis") %in% colnames(df))) return(FALSE)
  df <- df %>%
    dplyr::mutate(
      ovary = dplyr::coalesce(ovary, 0),
      testis = dplyr::coalesce(testis, 0),
      ovary_minus_testis = ovary - testis,
      gene_label = paste(group, gene_display_label, sep = " | ")
    )
  p <- ggplot(df, aes(x = testis, y = ovary, color = ovary_minus_testis)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55") +
    geom_point(alpha = 0.85, size = 2) +
    geom_text(aes(label = gene_display_label), check_overlap = TRUE, size = 2.4, vjust = -0.7) +
    facet_grid(dataset ~ sex) +
    scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    theme_bw(base_size = 9) +
    labs(title = title, x = paste0("Testiculo: media ", expression_log_label), y = paste0("Ovario: media ", expression_log_label), color = "Ovario - testiculo")
  ggsave(outfile, p, width = 12, height = max(5, min(14, length(unique(df$dataset)) * 2.2 + 3)), dpi = 300)
  TRUE
}

plot_group_aggregate_profile <- function(expr_long, outfile, title = "Perfil agregado por grupo") {
  df <- expr_long %>%
    dplyr::group_by(group, dataset, batch, condition, stage_class, stage_day, stage, tissue, sex) %>%
    dplyr::summarise(mean_z_log2TPM = mean(z_log2TPM, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(group, dataset, batch, condition, stage_class, stage_day, stage, tissue, sex)
  if (nrow(df) == 0) return(FALSE)
  line_df <- df %>%
    dplyr::group_by(group, dataset, batch, condition, tissue, sex) %>%
    dplyr::filter(dplyr::n_distinct(stage) > 1) %>%
    dplyr::ungroup()
  p <- ggplot(df, aes(x = stage, y = mean_z_log2TPM, color = condition, shape = sex,
                      group = interaction(dataset, batch, condition, tissue, sex))) +
    geom_hline(yintercept = 0, color = "gray75", linewidth = 0.3) +
    geom_point(size = 2) +
    facet_grid(group + dataset ~ tissue, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "Estagio detalhado", y = paste0("Media z-score ", expression_log_label), color = "Condicao", shape = "Sexo")
  if (nrow(line_df) > 0) p <- p + geom_line(data = line_df, alpha = 0.7)
  ggsave(outfile, p, width = 14, height = max(6, min(18, length(unique(df$group)) * length(unique(df$dataset)) * 1.8 + 3)), dpi = 300)
  TRUE
}

plot_deg_direction_summary <- function(deg_hits, gene_catalog, outfile, title = "Direcao DEG por contraste") {
  if (nrow(deg_hits) == 0) return(FALSE)
  lookup <- gene_display_lookup(gene_catalog)
  df <- deg_hits %>%
    dplyr::left_join(lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::mutate(
      gene_label = paste(group, gene_display_label, sep = " | "),
      direction = dplyr::case_when(
        significant & log2FoldChange_num > 0 ~ "up",
        significant & log2FoldChange_num < 0 ~ "down",
        TRUE ~ "not_sig"
      )
    )
  if (nrow(df) == 0) return(FALSE)
  if (length(unique(df$contrast_label)) > 90) {
    keep <- df %>%
      dplyr::group_by(contrast_label) %>%
      dplyr::summarise(best = suppressWarnings(min(padj_num, na.rm = TRUE)), .groups = "drop") %>%
      dplyr::arrange(best) %>%
      utils::head(90) %>%
      dplyr::pull(contrast_label)
    df <- df %>% dplyr::filter(contrast_label %in% keep)
  }
  p <- ggplot(df, aes(x = contrast_label, y = gene_label, fill = direction)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("up" = "#b2182b", "down" = "#2166ac", "not_sig" = "gray88")) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) +
    labs(title = title, x = "Projeto/modo | contraste", y = "Gene", fill = "Direcao")
  ggsave(outfile, p, width = 15, height = max(5, min(18, length(unique(df$gene_label)) * 0.35 + 3)), dpi = 300)
  TRUE
}

plot_deg_heatmap <- function(deg_hits, gene_catalog, outfile, title = "log2FC em contrastes DEG") {
  if (nrow(deg_hits) == 0) return(FALSE)
  gene_lookup <- gene_display_lookup(gene_catalog)
  df <- deg_hits %>%
    dplyr::group_by(gene_id, contrast_label) %>%
    dplyr::summarise(log2FoldChange_num = mean(log2FoldChange_num, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = contrast_label, values_from = log2FoldChange_num, values_fill = 0) %>%
    dplyr::left_join(gene_lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::mutate(label = paste(group, gene_display_label, sep = " | "))
  if (nrow(df) == 0 || ncol(df) <= 4) return(FALSE)
  mat <- as.matrix(df[, setdiff(colnames(df), c("gene_id", "gene_name", "gene_display_label", "description", "group", "label")), drop = FALSE])
  rownames(mat) <- df$label
  if (!heatmap_has_signal(mat)) return(FALSE)
  pheatmap::pheatmap(mat, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
                     cluster_rows = nrow(mat) > 1,
                     cluster_cols = ncol(mat) > 1,
                     border_color = NA, fontsize_row = 7, fontsize_col = 6,
                     main = title, filename = outfile,
                     width = 14, height = max(5, min(16, nrow(mat) * 0.32 + 3)))
  TRUE
}

plot_deg_context_tile <- function(deg_hits, gene_catalog, outfile, title = "Presenca DEG por contraste/projeto") {
  if (nrow(deg_hits) == 0) return(FALSE)
  gene_lookup <- gene_display_lookup(gene_catalog)
  df <- deg_hits %>%
    dplyr::left_join(gene_lookup, by = c("gene_id" = "matched_gene_id")) %>%
    dplyr::mutate(
      gene_label = paste(group, gene_display_label, sep = " | "),
      sig_label = ifelse(significant, "significativo", "nao_significativo")
    )
  if (nrow(df) == 0) return(FALSE)
  if (length(unique(df$contrast_label)) > 90) {
    keep <- df %>%
      dplyr::group_by(contrast_label) %>%
      dplyr::summarise(best = suppressWarnings(min(padj_num, na.rm = TRUE)), .groups = "drop") %>%
      dplyr::arrange(best) %>%
      utils::head(90) %>%
      dplyr::pull(contrast_label)
    df <- df %>% dplyr::filter(contrast_label %in% keep)
  }
  p <- ggplot(df, aes(x = contrast_label, y = gene_label, fill = log2FoldChange_num, alpha = sig_label)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", na.value = "gray90") +
    scale_alpha_manual(values = c("significativo" = 1, "nao_significativo" = 0.35)) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) +
    labs(title = title, x = "Projeto/modo | contraste", y = "Gene", fill = "log2FC", alpha = "")
  ggsave(outfile, p, width = 15, height = max(5, min(18, length(unique(df$gene_label)) * 0.35 + 3)), dpi = 300)
  TRUE
}

plot_gene_expression_boxplot <- function(df, outfile, label) {
  p <- ggplot(df, aes(x = interaction(tissue, sex, drop = TRUE), y = log2TPM, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = batch), width = 0.18, alpha = 0.55, size = 1.5) +
    facet_grid(dataset ~ stage, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = label, x = "Tecido.sexo", y = expression_log_label, fill = "Condicao", color = "Batch")
  ggsave(outfile, p, width = 13, height = 8, dpi = 300)
  TRUE
}

plot_gene_batch_boxplot <- function(df, outfile, label) {
  p <- ggplot(df, aes(x = batch, y = log2TPM, fill = dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = condition), width = 0.18, alpha = 0.55, size = 1.5) +
    facet_grid(tissue ~ sex, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Batch/projeto:", label), x = "Batch", y = expression_log_label, fill = "Projeto", color = "Condicao")
  ggsave(outfile, p, width = 13, height = 8, dpi = 300)
  TRUE
}

plot_gene_profile_line <- function(df, outfile, label) {
  profile_df <- df %>%
    dplyr::group_by(dataset, batch, condition, stage_class, stage_day, stage, tissue, sex) %>%
    dplyr::summarise(mean_log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop")
  line_df <- profile_df %>%
    dplyr::group_by(dataset, batch, condition, tissue, sex) %>%
    dplyr::filter(dplyr::n_distinct(stage) > 1) %>%
    dplyr::ungroup()
  p <- ggplot(profile_df, aes(x = stage, y = mean_log2TPM, color = condition, shape = sex,
                              group = interaction(dataset, batch, condition, tissue, sex))) +
    geom_point(size = 2) +
    facet_grid(dataset + batch ~ tissue, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Perfil medio:", label), x = "Estagio", y = expression_mean_log_label, color = "Condicao", shape = "Sexo")
  if (nrow(line_df) > 0) p <- p + geom_line(data = line_df, alpha = 0.75)
  ggsave(outfile, p, width = 14, height = 9, dpi = 300)
  TRUE
}

plot_gene_sample_tile <- function(df, outfile, label) {
  tile_df <- df %>%
    dplyr::mutate(sample_context = paste(dataset, batch, condition, stage, tissue, sex, sample_id, sep = " | ")) %>%
    dplyr::arrange(dataset, batch, condition, stage_class, stage_day, stage, tissue, sex, sample_id)
  p <- ggplot(tile_df, aes(x = sample_context, y = gene_display_label, fill = log2TPM)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "C") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) +
    labs(title = paste("Expressao por amostra:", label), x = "Amostra", y = "", fill = expression_log_label)
  ggsave(outfile, p, width = 15, height = 3.8, dpi = 300)
  TRUE
}

plot_gene_deg_lollipop <- function(deg_df, outfile, label) {
  if (nrow(deg_df) == 0) return(FALSE)
  df <- deg_df %>%
    dplyr::mutate(contrast_display = paste(deg_project, deg_mode, contrast, sep = " | ")) %>%
    dplyr::arrange(log2FoldChange_num)
  if (nrow(df) > 80) {
    df <- df %>%
      dplyr::arrange(padj_num) %>%
      utils::head(80) %>%
      dplyr::arrange(log2FoldChange_num)
  }
  p <- ggplot(df, aes(x = log2FoldChange_num, y = reorder(contrast_display, log2FoldChange_num), color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
    geom_segment(aes(x = 0, xend = log2FoldChange_num, yend = contrast_display), linewidth = 0.35) +
    geom_point(aes(size = neg_log10_padj), alpha = 0.85) +
    scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = "#b2182b")) +
    theme_bw(base_size = 8) +
    labs(title = paste("DEG:", label), x = "log2FC", y = "Contraste", color = "Significativo", size = "-log10(padj)")
  ggsave(outfile, p, width = 11, height = max(5, min(18, nrow(df) * 0.25 + 3)), dpi = 300)
  TRUE
}

plot_gene_deg_scatter <- function(deg_df, outfile, label) {
  if (nrow(deg_df) == 0) return(FALSE)
  p <- ggplot(deg_df, aes(x = log2FoldChange_num, y = neg_log10_padj, color = significant, shape = deg_project)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
    geom_point(size = 2.4, alpha = 0.85) +
    facet_wrap(~ deg_mode) +
    scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = "#b2182b")) +
    theme_bw(base_size = 10) +
    labs(title = paste("Contrastes DEG:", label), x = "log2FC", y = "-log10(padj)", color = "Significativo", shape = "Projeto")
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
  TRUE
}

plot_png_name <- function(base, suffix = "") paste0(base, suffix, ".png")
project_plot_suffix <- function(project) paste0("_project_", sanitize(project))
project_plot_path <- function(src, project) {
  sub("\\.png$", paste0(project_plot_suffix(project), ".png"), src, ignore.case = TRUE)
}

with_expression_unit <- function(unit, expr) {
  old_expression_unit <- expression_unit
  old_expression_log_label <- expression_log_label
  old_expression_mean_log_label <- expression_mean_log_label
  expression_unit <<- unit
  expression_log_label <<- paste0("log2(", unit, "+1)")
  expression_mean_log_label <<- paste0("Media log2(", unit, "+1)")
  on.exit({
    expression_unit <<- old_expression_unit
    expression_log_label <<- old_expression_log_label
    expression_mean_log_label <<- old_expression_mean_log_label
  }, add = TRUE)
  force(expr)
}

plot_group_outputs <- function(expr_long, expr_summary, deg_hits, gene_catalog, out_dir, suffix = "", title_suffix = "", include_deg = TRUE) {
  groups <- unique(expr_long$group)
  for (grp in groups) {
    group_dir <- file.path(out_dir, "groups", sanitize(grp))
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    expr_g <- expr_long %>% dplyr::filter(group == grp)
    summary_g <- expr_summary %>% dplyr::filter(group == grp)
    genes_g <- unique(expr_g$gene_id)
    deg_g <- deg_hits %>% dplyr::filter(gene_id %in% genes_g)
    catalog_g <- gene_catalog %>% dplyr::filter(group == grp)
    plot_or_skip(paste("group heatmap", grp, suffix), function() plot_expression_heatmap(summary_g, file.path(group_dir, plot_png_name("expression_heatmap", suffix)), paste("Grupo:", grp, title_suffix)))
    plot_or_skip(paste("group dotplot", grp, suffix), function() plot_expression_dotplot(summary_g, file.path(group_dir, plot_png_name("expression_dotplot", suffix)), paste("Grupo:", grp, title_suffix)))
    plot_or_skip(paste("group sample heatmap", grp, suffix), function() plot_group_sample_heatmap(expr_g, file.path(group_dir, plot_png_name("sample_heatmap", suffix)), paste("Amostras -", grp, title_suffix)))
    plot_or_skip(paste("group annotated sample heatmap", grp, suffix), function() plot_annotated_sample_heatmap(expr_g, file.path(group_dir, plot_png_name("sample_heatmap_annotated", suffix)), paste("Amostras anotadas -", grp, title_suffix)))
    plot_or_skip(paste("group gene correlation", grp, suffix), function() plot_gene_correlation(expr_g, file.path(group_dir, plot_png_name("gene_correlation", suffix)), paste("Correlacao entre genes -", grp, title_suffix)))
    plot_or_skip(paste("group PCA", grp, suffix), function() plot_sample_ordination(expr_g, file.path(group_dir, plot_png_name("sample_pca", suffix)), method = "pca", title = paste("PCA -", grp, title_suffix)))
    plot_or_skip(paste("group MDS", grp, suffix), function() plot_sample_ordination(expr_g, file.path(group_dir, plot_png_name("sample_mds", suffix)), method = "mds", title = paste("MDS -", grp, title_suffix)))
    plot_or_skip(paste("group aggregate profile", grp, suffix), function() plot_group_aggregate_profile(expr_g, file.path(group_dir, plot_png_name("aggregate_profile", suffix)), paste("Perfil agregado -", grp, title_suffix)))
    plot_or_skip(paste("group ovary/testis", grp, suffix), function() plot_ovary_testis_panel(summary_g, file.path(group_dir, plot_png_name("ovary_testis_panel", suffix)), paste("Ovario/testiculo -", grp, title_suffix)))
    plot_or_skip(paste("group batch", grp, suffix), function() plot_batch_project_boxplot(expr_g, file.path(group_dir, plot_png_name("batch_project_boxplot", suffix)), paste("Batch/projeto -", grp, title_suffix)))
    if (include_deg) {
      plot_or_skip(paste("group DEG heatmap", grp), function() plot_deg_heatmap(deg_g, catalog_g, file.path(group_dir, "deg_log2fc_heatmap.png"), paste("DEG -", grp)))
      plot_or_skip(paste("group DEG tile", grp), function() plot_deg_context_tile(deg_g, catalog_g, file.path(group_dir, "deg_context_tile.png"), paste("DEG por contraste -", grp)))
      plot_or_skip(paste("group DEG direction", grp), function() plot_deg_direction_summary(deg_g, catalog_g, file.path(group_dir, "deg_direction_summary.png"), paste("Direcao DEG -", grp)))
    }
  }
}

plot_gene_outputs <- function(expr_long, deg_hits, out_dir, suffix = "", title_suffix = "", include_deg = TRUE) {
  gene_keys <- expr_long %>% dplyr::distinct(group, gene_id, gene_name, gene_display_label, biotype)
  for (i in seq_len(nrow(gene_keys))) {
    key <- gene_keys[i, ]
    gene_dir <- file.path(out_dir, "genes", sanitize(key$group), sanitize(key$gene_id))
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
    df <- expr_long %>% dplyr::filter(group == key$group, gene_id == key$gene_id)
    deg_df <- deg_hits %>% dplyr::filter(gene_id == key$gene_id)
    label <- paste(key$group, key$gene_display_label, sep = " | ")
    corrected_label <- trimws(paste(label, title_suffix))
    plot_or_skip(paste("gene expression", key$gene_id, suffix), function() plot_gene_expression_boxplot(df, file.path(gene_dir, plot_png_name("expression_tissue_sex_condition", suffix)), corrected_label))
    plot_or_skip(paste("gene batch", key$gene_id, suffix), function() plot_gene_batch_boxplot(df, file.path(gene_dir, plot_png_name("expression_batch_project", suffix)), corrected_label))
    plot_or_skip(paste("gene profile", key$gene_id, suffix), function() plot_gene_profile_line(df, file.path(gene_dir, plot_png_name("expression_stage_profile", suffix)), corrected_label))
    plot_or_skip(paste("gene sample tile", key$gene_id, suffix), function() plot_gene_sample_tile(df, file.path(gene_dir, plot_png_name("expression_sample_tile", suffix)), corrected_label))
    if (include_deg) {
      plot_or_skip(paste("gene DEG lollipop", key$gene_id), function() plot_gene_deg_lollipop(deg_df, file.path(gene_dir, "deg_lollipop.png"), label))
      plot_or_skip(paste("gene DEG scatter", key$gene_id), function() plot_gene_deg_scatter(deg_df, file.path(gene_dir, "deg_scatter.png"), label))
    }
  }
}

plot_expression_scope_outputs <- function(expr_long, expr_summary, plots, title_prefix) {
  dirs <- unique(dirname(unlist(plots, use.names = FALSE)))
  for (dir in dirs) dir.create(file.path(out_dir, dir), recursive = TRUE, showWarnings = FALSE)
  invisible(plot_or_skip(paste(title_prefix, "expression heatmap"), function() plot_expression_heatmap(expr_summary, file.path(out_dir, plots$heatmap), paste(title_prefix, "- expressao media"))))
  invisible(plot_or_skip(paste(title_prefix, "expression dotplot"), function() plot_expression_dotplot(expr_summary, file.path(out_dir, plots$dotplot), paste(title_prefix, "- expressao media e fracao expressa"))))
  invisible(plot_or_skip(paste(title_prefix, "annotated sample heatmap"), function() plot_annotated_sample_heatmap(expr_long, file.path(out_dir, plots$annotated_sample_heatmap), paste(title_prefix, "- amostras anotadas"))))
  invisible(plot_or_skip(paste(title_prefix, "gene correlation"), function() plot_gene_correlation(expr_long, file.path(out_dir, plots$gene_correlation), paste(title_prefix, "- correlacao entre genes"))))
  invisible(plot_or_skip(paste(title_prefix, "sample PCA"), function() plot_sample_ordination(expr_long, file.path(out_dir, plots$sample_pca), method = "pca", title = paste(title_prefix, "- PCA"))))
  invisible(plot_or_skip(paste(title_prefix, "sample MDS"), function() plot_sample_ordination(expr_long, file.path(out_dir, plots$sample_mds), method = "mds", title = paste(title_prefix, "- MDS"))))
  invisible(plot_or_skip(paste(title_prefix, "tissue/sex heatmap"), function() plot_tissue_sex_heatmap(expr_summary, file.path(out_dir, plots$tissue_sex_heatmap), paste(title_prefix, "- tecido e sexo"))))
  invisible(plot_or_skip(paste(title_prefix, "ovary/testis"), function() plot_ovary_testis_panel(expr_summary, file.path(out_dir, plots$ovary_testis), paste(title_prefix, "- ovario versus testiculo"))))
  invisible(plot_or_skip(paste(title_prefix, "group aggregate profile"), function() plot_group_aggregate_profile(expr_long, file.path(out_dir, plots$group_aggregate), paste(title_prefix, "- perfil agregado"))))
  invisible(plot_or_skip(paste(title_prefix, "batch/project"), function() plot_batch_project_boxplot(expr_long, file.path(out_dir, plots$batch_project), paste(title_prefix, "- batch/projeto"))))
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

table_to_html <- function(df, max_rows = 30) {
  if (is.null(df)) return("<p><em>Nenhum registro.</em></p>")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    if (ncol(df) == 0) return("<p><em>Nenhum registro.</em></p>")
    header <- paste0("<tr>", paste0("<th>", html_escape(colnames(df)), "</th>", collapse = ""), "</tr>")
    empty_row <- paste0("<tr><td colspan='", ncol(df), "'><em>Nenhum registro encontrado.</em></td></tr>")
    return(paste0("<div class='table-wrap'><table>", header, empty_row, "</table></div>"))
  }
  truncated <- nrow(df) > max_rows
  total_rows <- nrow(df)
  if (truncated) df <- df[seq_len(max_rows), , drop = FALSE]
  header <- paste0("<tr>", paste0("<th>", html_escape(colnames(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(row) {
    row_text <- paste(row, collapse = " ")
    paste0(
      "<tr class='searchable table-row' data-kind='table' data-search='", html_escape(tolower(row_text)), "'>",
      paste0("<td>", html_escape(row), "</td>", collapse = ""),
      "</tr>"
    )
  })
  note <- if (truncated) paste0("<p class='table-note'>Mostrando ", max_rows, " de ", total_rows, " registros. A tabela completa esta em <code>tables/</code>.</p>") else ""
  paste0("<div class='table-wrap'><table>", header, paste(rows, collapse = "\n"), "</table></div>", note)
}

figure_explanation <- function(caption) {
  caption_l <- tolower(caption)
  dplyr::case_when(
    grepl("pca", caption_l) ~ "Resume a variacao entre amostras usando os genes exibidos. Separacoes por projeto podem indicar batch ou diferencas biologicas fortes.",
    grepl("mds", caption_l) ~ "Mostra distancias entre amostras; pontos proximos tem perfis de expressao parecidos.",
    grepl("correlacao", caption_l) ~ "Indica se genes variam juntos entre amostras. Correlacoes altas sugerem perfis coordenados.",
    grepl("dotplot|fracao expressa", caption_l) ~ "Combina intensidade media de expressao com a proporcao de amostras em que o gene esta expresso.",
    grepl("heatmap gene x amostra|amostra", caption_l) ~ "Mostra expressao por amostra individual, util para checar outliers e consistencia entre replicatas.",
    grepl("heatmap|expressao media", caption_l) ~ "Resume expressao media por contexto; cores mais intensas indicam maior expressao relativa.",
    grepl("batch|projeto", caption_l) ~ "Ajuda a avaliar se o sinal acompanha projeto ou batch, o que exige cautela na interpretacao biologica.",
    grepl("deg|log2fc|contraste", caption_l) ~ "Resume efeitos diferenciais. log2FC indica direcao/tamanho do efeito e padj indica suporte estatistico.",
    grepl("perfil agregado|perfil medio", caption_l) ~ "Mostra tendencias medias do grupo ao longo dos contextos disponiveis.",
    TRUE ~ "Figura exploratoria para revisar padroes de expressao e consistencia dos resultados."
  )
}

img_tag <- function(src, caption, projects = character()) {
  if (!file.exists(file.path(out_dir, src))) return("")
  explanation <- figure_explanation(caption)
  search_text <- paste(caption, explanation, src)
  projects <- projects[projects != "" & !is.na(projects)]
  project_attr <- if (length(projects) > 0) paste0(" data-projects='", html_escape(paste(sort(unique(projects)), collapse = "|")), "'") else ""
  paste0(
    "<figure class='searchable report-figure' data-kind='figure'", project_attr, " data-search='", html_escape(tolower(search_text)), "'>",
    "<img src='", gsub("\\\\", "/", src), "' alt='", html_escape(caption), "'>",
    "<figcaption><strong>", html_escape(caption), "</strong><span>", html_escape(explanation), "</span></figcaption>",
    "</figure>"
  )
}

corrected_plot_path <- function(src) {
  sub("\\.png$", "_corrected.png", src, ignore.case = TRUE)
}

corrected_img_tag <- function(src, caption) {
  img_tag(corrected_plot_path(src), paste(caption, "- corrigido 055"), projects = "all_projects")
}

project_img_tag <- function(src, caption, project) {
  img_tag(project_plot_path(src, project), paste(caption, "-", project), projects = project)
}

figure_details <- function(title, body, count = NA_integer_, open = FALSE, note = "") {
  count_html <- if (!is.na(count)) paste0("<span>", count, "</span>") else ""
  if (body == "" || (!is.na(count) && count == 0)) body <- "<p><em>Nenhuma figura encontrada.</em></p>"
  paste0(
    "<details class='data-block figure-block' ", if (open) "open" else "", ">",
    "<summary>", html_escape(title), count_html, "</summary>",
    ifelse(note == "", "", paste0("<p class='description'>", html_escape(note), "</p>")),
    body,
    "</details>"
  )
}

count_existing_plots <- function(paths) {
  sum(file.exists(file.path(out_dir, unlist(paths, use.names = FALSE))))
}

expression_figure_grid <- function(plots, suffix = "", caption_suffix = "") {
  caption_suffix <- trimws(caption_suffix)
  cap <- function(x) trimws(paste(x, caption_suffix))
  paste0(
    "<div class='figure-grid'>",
    img_tag(if (suffix == "") plots$heatmap else corrected_plot_path(plots$heatmap), cap("Expressao media integrada por contexto")),
    img_tag(if (suffix == "") plots$dotplot else corrected_plot_path(plots$dotplot), cap("Expressao media e fracao expressa")),
    img_tag(if (suffix == "") plots$annotated_sample_heatmap else corrected_plot_path(plots$annotated_sample_heatmap), cap("Heatmap gene x amostra com anotacoes")),
    img_tag(if (suffix == "") plots$gene_correlation else corrected_plot_path(plots$gene_correlation), cap("Correlacao de expressao entre genes")),
    img_tag(if (suffix == "") plots$sample_pca else corrected_plot_path(plots$sample_pca), cap("PCA das amostras usando os genes de interesse")),
    img_tag(if (suffix == "") plots$sample_mds else corrected_plot_path(plots$sample_mds), cap("MDS das amostras usando os genes de interesse")),
    img_tag(if (suffix == "") plots$tissue_sex_heatmap else corrected_plot_path(plots$tissue_sex_heatmap), cap("Padroes por tecido e sexo")),
    img_tag(if (suffix == "") plots$ovary_testis else corrected_plot_path(plots$ovary_testis), cap("Comparacao ovario versus testiculo")),
    img_tag(if (suffix == "") plots$group_aggregate else corrected_plot_path(plots$group_aggregate), cap("Perfil agregado por grupo")),
    img_tag(if (suffix == "") plots$batch_project else corrected_plot_path(plots$batch_project), cap("Distribuicao de expressao por batch e projeto")),
    "</div>"
  )
}

deg_figure_grid <- function(plots) {
  paste0(
    "<div class='figure-grid'>",
    img_tag(plots$deg_heatmap, "log2FC integrado nos contrastes DEG"),
    img_tag(plots$deg_tile, "Sinais DEG por contraste/projeto"),
    img_tag(plots$deg_direction, "Direcao DEG por contraste/projeto"),
    "</div>"
  )
}

project_plot_paths <- function(project) {
  project_dir <- file.path("plots", "by_project", sanitize(project))
  list(
    heatmap = file.path(project_dir, "expression_heatmap.png"),
    dotplot = file.path(project_dir, "expression_dotplot.png"),
    annotated_sample_heatmap = file.path(project_dir, "sample_heatmap_annotated.png"),
    gene_correlation = file.path(project_dir, "gene_correlation.png"),
    sample_pca = file.path(project_dir, "sample_pca.png"),
    sample_mds = file.path(project_dir, "sample_mds.png"),
    tissue_sex_heatmap = file.path(project_dir, "tissue_sex_heatmap.png"),
    ovary_testis = file.path(project_dir, "ovary_testis_panel.png"),
    group_aggregate = file.path(project_dir, "aggregate_profile.png"),
    batch_project = file.path(project_dir, "batch_project_boxplot.png")
  )
}

project_figure_sections <- function(projects) {
  sections <- vapply(projects, function(project) {
    plots <- project_plot_paths(project)
    body <- expression_figure_grid(plots, caption_suffix = paste("-", project, "raw"))
    figure_details(
      paste("Projeto", project, "- raw"),
      body,
      count_existing_plots(plots),
      open = FALSE,
      note = "Visualizacoes restritas a um unico projeto; nao misturam amostras do all_projects nem matriz corrigida."
    )
  }, character(1))
  paste(sections, collapse = "\n")
}

empty_plot_message <- function(count) {
  if (!is.na(count) && count == 0) "<p><em>Nenhuma figura encontrada.</em></p>" else ""
}

scope_tabs_html <- function(global_plots, projects, n_corrected_global_figures) {
  tabs <- list()
  all_raw_count <- count_existing_plots(global_plots[setdiff(names(global_plots), c("deg_heatmap", "deg_tile", "deg_direction"))])
  tabs[[length(tabs) + 1]] <- list(
    id = "all_raw",
    label = "All raw",
    count = all_raw_count,
    note = "Integra todos os projetos usando a matriz principal do relatorio. Nao inclui figuras corrigidas.",
    body = if (all_raw_count == 0) empty_plot_message(all_raw_count) else expression_figure_grid(global_plots, caption_suffix = "- all_projects raw")
  )
  tabs[[length(tabs) + 1]] <- list(
    id = "all_corrected",
    label = "All corrected",
    count = n_corrected_global_figures,
    note = "Gerado apenas quando a etapa 055 produziu matriz corrigida all_projects; nao mistura projetos individuais.",
    body = if (n_corrected_global_figures == 0) empty_plot_message(n_corrected_global_figures) else expression_figure_grid(global_plots, suffix = "_corrected", caption_suffix = "- all_projects corrected 055")
  )
  for (project in projects) {
    plots <- project_plot_paths(project)
    project_count <- count_existing_plots(plots)
    tabs[[length(tabs) + 1]] <- list(
      id = paste0("project_", sanitize(project)),
      label = project,
      count = project_count,
      note = "Visualizacoes restritas a um unico projeto; nao misturam amostras do all_projects nem matriz corrigida.",
      body = if (project_count == 0) empty_plot_message(project_count) else expression_figure_grid(plots, caption_suffix = paste("-", project, "raw"))
    )
  }

  group_id <- "scope_tabs"
  buttons <- paste(vapply(seq_along(tabs), function(i) {
    tab <- tabs[[i]]
    target <- paste(group_id, tab$id, sep = "_")
    paste0(
      "<button type='button' class='tab-button", if (i == 1) " active" else "", "' data-tab-group='", group_id, "' data-tab-target='", target, "'>",
      "<span>", html_escape(tab$label), "</span><small>", tab$count, "</small>",
      "</button>"
    )
  }, character(1)), collapse = "\n")

  panels <- paste(vapply(seq_along(tabs), function(i) {
    tab <- tabs[[i]]
    target <- paste(group_id, tab$id, sep = "_")
    paste0(
      "<section id='", target, "' class='tab-panel", if (i == 1) " active" else "", "' data-tab-panel='", group_id, "'>",
      "<p class='description'>", html_escape(tab$note), "</p>",
      tab$body,
      "</section>"
    )
  }, character(1)), collapse = "\n")

  paste0(
    "<div class='tabs' data-tabs='", group_id, "'>",
    "<div class='tab-list' role='tablist'>", buttons, "</div>",
    panels,
    "</div>"
  )
}

copy_analysis_plots <- function(root, analysis_name, caption_map = character()) {
  empty <- tibble::tibble(src = character(), caption = character(), analysis = character(), project = character(), source_file = character())
  if (root == "" || !dir.exists(root)) return(empty)
  root_norm <- normalizePath(root, winslash = "/", mustWork = FALSE)
  files <- list.files(root, pattern = "\\.(png|jpg|jpeg|svg)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) return(empty)

  rows <- lapply(files, function(path) {
    path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
    rel <- gsub("\\\\", "/", sub(paste0("^", root_norm, "/?"), "", path_norm))
    rel_parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    project <- ifelse(length(rel_parts) >= 1 && rel_parts[1] != ".", rel_parts[1], "unknown")
    base <- basename(path)
    stem <- tools::file_path_sans_ext(base)
    caption_base <- if (base %in% names(caption_map)) caption_map[[base]] else gsub("_", " ", stem, fixed = TRUE)
    target_rel <- file.path("plots", tolower(analysis_name), sanitize(project), base)
    target_abs <- file.path(out_dir, target_rel)
    dir.create(dirname(target_abs), recursive = TRUE, showWarnings = FALSE)
    ok <- suppressWarnings(file.copy(path, target_abs, overwrite = TRUE))
    if (!isTRUE(ok) && !file.exists(target_abs)) return(NULL)
    tibble::tibble(
      src = gsub("\\\\", "/", target_rel),
      caption = paste(analysis_name, project, "-", caption_base),
      analysis = analysis_name,
      project = project,
      source_file = rel
    )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(empty)
  out %>% dplyr::arrange(project, source_file)
}

figure_gallery_html <- function(plot_rows) {
  if (is.null(plot_rows) || nrow(plot_rows) == 0) return("<p><em>Nenhuma figura encontrada.</em></p>")
  paste0(
    "<div class='figure-grid'>",
    paste(vapply(seq_len(nrow(plot_rows)), function(i) img_tag(plot_rows$src[i], plot_rows$caption[i]), character(1)), collapse = "\n"),
    "</div>"
  )
}

details_block <- function(title, df, max_rows = 50, open = FALSE) {
  count <- if (is.null(df)) 0 else nrow(df)
  paste0(
    "<details class='data-block' ", if (open) "open" else "", ">",
    "<summary>", html_escape(title), "<span>", count, "</span></summary>",
    table_to_html(df, max_rows),
    "</details>"
  )
}

coalesce_text <- function(x, fallback = "") {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- fallback
  x
}

write_html_report <- function(path, title, catalog, gene_summary, stage_specificity, deg_hits, dtu_hits, splicing_hits, wgcna_hits, mfuzz_hits, global_plots, expression_unit = "TPM", corrected_expression_file = "", corrected_expression_unit = "batch-corrected counts") {
  if (is.na(expression_unit) || expression_unit == "") expression_unit <- "TPM"
  if (is.null(stage_specificity)) {
    stage_specificity <- tibble::tibble(
      group = character(),
      gene_id = character(),
      gene_display_label = character(),
      dominant_stage = character(),
      stage_specificity_candidate = logical()
    )
  }
  batch_plot_gallery <- copy_analysis_plots(
    batch_root,
    "BatchCorrection",
    c("batch_pca_before_after.png" = "PCA antes/depois da correcao de batch")
  )
  wgcna_plot_gallery <- copy_analysis_plots(
    wgcna_root,
    "WGCNA",
    c(
      "soft_threshold.png" = "selecao do soft-threshold",
      "module_sizes.png" = "tamanho dos modulos",
      "module_trait_heatmap.png" = "correlacao modulo-traco"
    )
  )
  mfuzz_plot_gallery <- copy_analysis_plots(
    mfuzz_root,
    "Mfuzz",
    c("cluster_centers.png" = "centroides dos clusters")
  )
  n_found_genes <- catalog %>%
    dplyr::filter(found_in_expression_matrix %in% TRUE) %>%
    dplyr::distinct(matched_gene_id) %>%
    nrow()
  n_annotated_genes <- catalog %>%
    dplyr::filter(location != "" | biotype != "Unknown" | description != "") %>%
    dplyr::distinct(matched_gene_id) %>%
    nrow()
  n_deg_sig_genes <- if (nrow(deg_hits) > 0) {
    deg_hits %>%
      dplyr::filter(significant %in% TRUE) %>%
      dplyr::distinct(gene_id) %>%
      nrow()
  } else {
    0
  }
  n_dtu_sig_genes <- if (nrow(dtu_hits) > 0) dtu_hits %>% dplyr::distinct(gene_id) %>% nrow() else 0
  n_dtu_sig_transcripts <- if (nrow(dtu_hits) > 0) dtu_hits %>% dplyr::distinct(gene_id, transcript_id) %>% nrow() else 0
  n_dtu_report_genes <- if (nrow(dtu_hits) > 0) {
    dtu_hits %>%
      dplyr::filter(gene_id %in% catalog$matched_gene_id) %>%
      dplyr::distinct(gene_id) %>%
      nrow()
  } else {
    0
  }
  n_splicing_sig_genes <- if (nrow(splicing_hits) > 0) splicing_hits %>% dplyr::distinct(gene_id) %>% nrow() else 0
  n_wgcna_genes <- if (nrow(wgcna_hits) > 0) wgcna_hits %>% dplyr::distinct(gene_id) %>% nrow() else 0
  n_wgcna_hub_genes <- if (nrow(wgcna_hits) > 0) wgcna_hits %>% dplyr::filter(is_hub %in% TRUE) %>% dplyr::distinct(gene_id) %>% nrow() else 0
  n_mfuzz_genes <- if (nrow(mfuzz_hits) > 0) mfuzz_hits %>% dplyr::distinct(gene_id) %>% nrow() else 0
  n_stage_specific_genes <- if (nrow(stage_specificity) > 0) {
    stage_specificity %>%
      dplyr::filter(stage_specificity_candidate %in% TRUE) %>%
      dplyr::distinct(gene_id) %>%
      nrow()
  } else {
    0
  }
  n_corrected_figures <- length(list.files(out_dir, pattern = "_corrected\\.png$", recursive = TRUE, full.names = TRUE))
  corrected_global_expression_paths <- vapply(global_plots[setdiff(names(global_plots), c("deg_heatmap", "deg_tile", "deg_direction"))], corrected_plot_path, character(1))
  n_corrected_global_figures <- count_existing_plots(as.list(corrected_global_expression_paths))
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M")

  split_projects <- function(x) {
    x <- paste(x, collapse = "; ")
    vals <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
    vals <- vals[vals != "" & !is.na(vals)]
    unique(vals)
  }
  projects_for_gene <- function(group, gene_id) {
    summary_projects <- gene_summary %>%
      dplyr::filter(.data$group == .env$group, .data$gene_id == .env$gene_id) %>%
      dplyr::select(dplyr::any_of(c("datasets", "deg_projects", "dtu_projects", "splicing_projects", "wgcna_projects", "mfuzz_projects")))
    deg_projects <- deg_hits %>%
      dplyr::filter(.data$gene_id == .env$gene_id) %>%
      dplyr::pull(deg_project)
    dtu_projects <- dtu_hits %>%
      dplyr::filter(.data$gene_id == .env$gene_id) %>%
      dplyr::pull(dtu_project)
    splicing_projects <- splicing_hits %>%
      dplyr::filter(.data$gene_id == .env$gene_id) %>%
      dplyr::pull(splicing_project)
    wgcna_projects <- wgcna_hits %>%
      dplyr::filter(.data$gene_id == .env$gene_id) %>%
      dplyr::pull(wgcna_project)
    mfuzz_projects <- mfuzz_hits %>%
      dplyr::filter(.data$gene_id == .env$gene_id) %>%
      dplyr::pull(mfuzz_project)
    sort(unique(c(split_projects(unlist(summary_projects)), deg_projects, dtu_projects, splicing_projects, wgcna_projects, mfuzz_projects)))
  }
  project_attr <- function(projects) {
    projects <- projects[projects != "" & !is.na(projects)]
    html_escape(paste(sort(unique(projects)), collapse = "|"))
  }
  all_projects <- sort(unique(c(
    split_projects(gene_summary$datasets),
    split_projects(gene_summary$deg_projects),
    split_projects(gene_summary$dtu_projects),
    split_projects(gene_summary$splicing_projects),
    split_projects(gene_summary$wgcna_projects),
    split_projects(gene_summary$mfuzz_projects),
    if (nrow(deg_hits) > 0) deg_hits$deg_project else character(),
    if (nrow(dtu_hits) > 0) dtu_hits$dtu_project else character(),
    if (nrow(splicing_hits) > 0) splicing_hits$splicing_project else character(),
    if (nrow(wgcna_hits) > 0) wgcna_hits$wgcna_project else character(),
    if (nrow(mfuzz_hits) > 0) mfuzz_hits$mfuzz_project else character()
  )))
  project_options <- paste(vapply(all_projects, function(project) {
    paste0("<option value='", html_escape(project), "'>", html_escape(project), "</option>")
  }, character(1)), collapse = "")
  dataset_projects <- sort(setdiff(split_projects(gene_summary$datasets), c("", "unknown", NA)))

  dtu_overview_table <- dtu_hits %>%
    dplyr::mutate(
      group = coalesce_text(group, "fora do genes.txt"),
      gene_display_label = coalesce_text(gene_display_label, gene_id)
    ) %>%
    dplyr::select(group, gene_id, gene_display_label, transcript_id, dtu_project, dtu_scope, variable, max_level, min_level, delta_usage, padj, method, source_file) %>%
    dplyr::arrange(padj, gene_id, transcript_id)

  wgcna_overview_table <- wgcna_hits %>%
    dplyr::mutate(gene_display_label = coalesce_text(gene_display_label, gene_id)) %>%
    dplyr::select(group, gene_id, gene_display_label, wgcna_project, wgcna_scope, module_color, module_label, is_hub, kME, abs_kME, module_rank, source_file) %>%
    dplyr::arrange(wgcna_project, module_color, module_rank, gene_id)

  mfuzz_overview_table <- mfuzz_hits %>%
    dplyr::mutate(gene_display_label = coalesce_text(gene_display_label, gene_id)) %>%
    dplyr::select(group, gene_id, gene_display_label, mfuzz_project, mfuzz_scope, cluster, membership, source_file) %>%
    dplyr::arrange(mfuzz_project, cluster, dplyr::desc(membership), gene_id)

  stage_specificity_table <- stage_specificity %>%
    dplyr::mutate(gene_display_label = coalesce_text(gene_display_label, gene_id)) %>%
    dplyr::select(
      group, gene_id, gene_display_label, dominant_stage, n_stages_tested,
      stage_specificity_candidate, stage_specificity_evidence,
      stage_specificity_tau, stage_specificity_max_fraction,
      stage_specificity_fold_vs_second, expression_stage_specific,
      n_stage_deg_significant, min_stage_deg_padj, max_abs_stage_deg_log2FC,
      max_mfuzz_membership, top_mfuzz_cluster
    ) %>%
    dplyr::arrange(dplyr::desc(stage_specificity_candidate), dominant_stage, dplyr::desc(stage_specificity_tau), group, gene_display_label)

  group_links <- paste(vapply(unique(catalog$group), function(grp) {
    paste0("<li><a href='#group_", sanitize(grp), "'>", html_escape(grp), "</a></li>")
  }, character(1)), collapse = "\n")

  gene_index <- paste(vapply(seq_len(nrow(catalog)), function(i) {
    row <- catalog[i, ]
    gene_projects <- projects_for_gene(row$group, row$matched_gene_id)
    search_text <- paste(row$group, row$query, row$matched_gene_id, row$gene_name, row$gene_display_label, row$biotype, row$description, row$chromosome, row$location, paste(gene_projects, collapse = " "))
    paste0(
      "<a class='searchable gene-chip' data-kind='gene' data-projects='", project_attr(gene_projects), "' data-search='", html_escape(tolower(search_text)), "' href='#gene_", sanitize(row$group), "_", sanitize(row$matched_gene_id), "'>",
      "<span>", html_escape(row$gene_display_label), "</span>",
      "<small>", html_escape(row$group), "</small>",
      "</a>"
    )
  }, character(1)), collapse = "\n")

  group_sections <- paste(vapply(unique(catalog$group), function(grp) {
    group_dir <- file.path("groups", sanitize(grp))
    group_catalog <- catalog %>%
      dplyr::filter(group == grp) %>%
      dplyr::select(group, query, query_display, matched_gene_id, gene_name, gene_display_label, biotype, chromosome, gene_start, gene_end, strand, location, found_in_expression_matrix)
    group_projects <- sort(unique(unlist(lapply(group_catalog$matched_gene_id, function(gene_id) projects_for_gene(grp, gene_id)))))
    group_search <- paste(group_catalog$group, group_catalog$query, group_catalog$matched_gene_id, group_catalog$gene_name, group_catalog$gene_display_label, group_catalog$biotype, group_catalog$chromosome, group_catalog$location, paste(group_projects, collapse = " "), collapse = " ")
    group_raw_figures <- paste0(
      img_tag(file.path(group_dir, "expression_heatmap.png"), "Expressao media por contexto biologico, projeto e batch", projects = "all_projects"),
      img_tag(file.path(group_dir, "expression_dotplot.png"), "Media de expressao e fracao expressa por contexto", projects = "all_projects"),
      img_tag(file.path(group_dir, "sample_heatmap.png"), "Expressao nas amostras individuais", projects = "all_projects"),
      img_tag(file.path(group_dir, "sample_heatmap_annotated.png"), "Heatmap gene x amostra com anotacoes de projeto, batch e biologia", projects = "all_projects"),
      img_tag(file.path(group_dir, "gene_correlation.png"), "Correlacao de expressao entre genes do grupo", projects = "all_projects"),
      img_tag(file.path(group_dir, "sample_pca.png"), "PCA das amostras usando apenas genes do grupo", projects = "all_projects"),
      img_tag(file.path(group_dir, "sample_mds.png"), "MDS das amostras usando apenas genes do grupo", projects = "all_projects"),
      img_tag(file.path(group_dir, "aggregate_profile.png"), "Perfil agregado medio do grupo", projects = "all_projects"),
      img_tag(file.path(group_dir, "ovary_testis_panel.png"), "Comparacao ovario versus testiculo", projects = "all_projects"),
      img_tag(file.path(group_dir, "batch_project_boxplot.png"), "Distribuicao de expressao por batch e projeto", projects = "all_projects"),
      img_tag(file.path(group_dir, "deg_log2fc_heatmap.png"), "log2FC dos genes do grupo nos contrastes DEG", projects = "all_projects"),
      img_tag(file.path(group_dir, "deg_context_tile.png"), "Consistencia dos sinais DEG por contraste/projeto", projects = "all_projects"),
      img_tag(file.path(group_dir, "deg_direction_summary.png"), "Direcao DEG por contraste/projeto", projects = "all_projects")
    )
    group_project_sections <- paste(vapply(group_projects[group_projects != "all_projects" & group_projects != "unknown"], function(project) {
      project_figures <- paste0(
        project_img_tag(file.path(group_dir, "expression_heatmap.png"), "Expressao media por contexto biologico, projeto e batch", project),
        project_img_tag(file.path(group_dir, "expression_dotplot.png"), "Media de expressao e fracao expressa por contexto", project),
        project_img_tag(file.path(group_dir, "sample_heatmap.png"), "Expressao nas amostras individuais", project),
        project_img_tag(file.path(group_dir, "sample_heatmap_annotated.png"), "Heatmap gene x amostra com anotacoes de projeto, batch e biologia", project),
        project_img_tag(file.path(group_dir, "gene_correlation.png"), "Correlacao de expressao entre genes do grupo", project),
        project_img_tag(file.path(group_dir, "sample_pca.png"), "PCA das amostras usando apenas genes do grupo", project),
        project_img_tag(file.path(group_dir, "sample_mds.png"), "MDS das amostras usando apenas genes do grupo", project),
        project_img_tag(file.path(group_dir, "aggregate_profile.png"), "Perfil agregado medio do grupo", project),
        project_img_tag(file.path(group_dir, "ovary_testis_panel.png"), "Comparacao ovario versus testiculo", project),
        project_img_tag(file.path(group_dir, "batch_project_boxplot.png"), "Distribuicao de expressao por batch e projeto", project),
        project_img_tag(file.path(group_dir, "deg_log2fc_heatmap.png"), "log2FC dos genes do grupo nos contrastes DEG", project),
        project_img_tag(file.path(group_dir, "deg_context_tile.png"), "Consistencia dos sinais DEG por contraste/projeto", project),
        project_img_tag(file.path(group_dir, "deg_direction_summary.png"), "Direcao DEG por contraste/projeto", project)
      )
      project_count <- if (project_figures == "") 0L else NA_integer_
      project_body <- if (project_figures == "") "" else paste0("<div class='figure-grid'>", project_figures, "</div>")
      figure_details(paste("Projeto", project, "raw - grupo"), project_body, count = project_count, open = FALSE, note = paste("Figuras filtradas para", project, "; nao incluem amostras de outros projetos/batches."))
    }, character(1)), collapse = "")
    group_corrected_figures <- paste0(
      corrected_img_tag(file.path(group_dir, "expression_heatmap.png"), "Expressao media por contexto biologico, projeto e batch"),
      corrected_img_tag(file.path(group_dir, "expression_dotplot.png"), "Media de expressao e fracao expressa por contexto"),
      corrected_img_tag(file.path(group_dir, "sample_heatmap.png"), "Expressao nas amostras individuais"),
      corrected_img_tag(file.path(group_dir, "sample_heatmap_annotated.png"), "Heatmap gene x amostra com anotacoes de projeto, batch e biologia"),
      corrected_img_tag(file.path(group_dir, "gene_correlation.png"), "Correlacao de expressao entre genes do grupo"),
      corrected_img_tag(file.path(group_dir, "sample_pca.png"), "PCA das amostras usando apenas genes do grupo"),
      corrected_img_tag(file.path(group_dir, "sample_mds.png"), "MDS das amostras usando apenas genes do grupo"),
      corrected_img_tag(file.path(group_dir, "aggregate_profile.png"), "Perfil agregado medio do grupo"),
      corrected_img_tag(file.path(group_dir, "ovary_testis_panel.png"), "Comparacao ovario versus testiculo"),
      corrected_img_tag(file.path(group_dir, "batch_project_boxplot.png"), "Distribuicao de expressao por batch e projeto")
    )
    paste0(
      "<section class='searchable group-section' data-kind='group' data-projects='", project_attr(group_projects), "' data-search='", html_escape(tolower(group_search)), "' id='group_", sanitize(grp), "'><h2>Grupo: ", html_escape(grp), "</h2>",
      details_block("Genes deste grupo", group_catalog, 100, TRUE),
      figure_details("All projects raw - grupo", paste0("<div class='figure-grid'>", group_raw_figures, "</div>"), open = TRUE, note = "Figuras geradas com a matriz principal do relatorio."),
      group_project_sections,
      figure_details("All projects corrected 055 - grupo", paste0("<div class='figure-grid'>", group_corrected_figures, "</div>"), open = FALSE, note = "Figuras geradas somente quando a matriz corrigida do 055 existe; ficam separadas para nao misturar escopos."),
      "</section>"
    )
  }, character(1)), collapse = "\n")

  gene_sections <- paste(vapply(seq_len(nrow(catalog)), function(i) {
    row <- catalog[i, ]
    gene_dir <- file.path("genes", sanitize(row$group), sanitize(row$matched_gene_id))
    deg_table <- deg_hits %>%
      dplyr::filter(gene_id == row$matched_gene_id) %>%
      dplyr::select(gene_display_label, deg_project, deg_mode, contrast, log2FoldChange_num, padj_num, significant) %>%
      dplyr::arrange(padj_num)
    dtu_table <- dtu_hits %>%
      dplyr::filter(gene_id == row$matched_gene_id) %>%
      dplyr::select(gene_display_label, dtu_project, dtu_scope, variable, transcript_id, max_level, min_level, delta_usage, padj, method) %>%
      dplyr::arrange(padj, dplyr::desc(delta_usage))
    splicing_table <- splicing_hits %>%
      dplyr::filter(gene_id == row$matched_gene_id) %>%
      dplyr::select(gene_display_label, splicing_project, splicing_variable, splicing_contrast, event_type, FDR, IncLevelDifference) %>%
      dplyr::arrange(FDR, dplyr::desc(abs(IncLevelDifference)))
    wgcna_table <- wgcna_hits %>%
      dplyr::filter(gene_id == row$matched_gene_id) %>%
      dplyr::select(gene_display_label, wgcna_project, wgcna_scope, module_color, module_label, is_hub, kME, abs_kME, module_rank) %>%
      dplyr::arrange(wgcna_project, module_color, module_rank)
    mfuzz_table <- mfuzz_hits %>%
      dplyr::filter(gene_id == row$matched_gene_id) %>%
      dplyr::select(gene_display_label, mfuzz_project, mfuzz_scope, cluster, membership) %>%
      dplyr::arrange(mfuzz_project, cluster, dplyr::desc(membership))
    stage_table <- stage_specificity %>%
      dplyr::filter(group == row$group, gene_id == row$matched_gene_id) %>%
      dplyr::select(
        gene_display_label, dominant_stage, n_stages_tested,
        stage_specificity_candidate, stage_specificity_evidence,
        stage_specificity_tau, stage_specificity_max_fraction,
        stage_specificity_fold_vs_second, expression_stage_specific,
        n_stage_deg_significant, min_stage_deg_padj, max_abs_stage_deg_log2FC,
        max_mfuzz_membership, top_mfuzz_cluster
      )
    gene_projects <- projects_for_gene(row$group, row$matched_gene_id)
    gene_search <- paste(row$group, row$query, row$matched_gene_id, row$gene_name, row$gene_display_label, row$biotype, row$description, row$chromosome, row$location, paste(gene_projects, collapse = " "), paste(stage_table$dominant_stage, collapse = " "), paste(stage_table$stage_specificity_evidence, collapse = " "), paste(deg_table$contrast, collapse = " "), paste(dtu_table$transcript_id, collapse = " "), paste(splicing_table$splicing_contrast, collapse = " "), paste(wgcna_table$module_color, collapse = " "), paste(mfuzz_table$cluster, collapse = " "))
    gene_raw_figures <- paste0(
      img_tag(file.path(gene_dir, "expression_tissue_sex_condition.png"), "Expressao por tecido, sexo, condicao, estagio e projeto", projects = "all_projects"),
      img_tag(file.path(gene_dir, "expression_batch_project.png"), "Expressao por batch/projeto", projects = "all_projects"),
      img_tag(file.path(gene_dir, "expression_stage_profile.png"), "Perfil medio por estagio, tecido, batch e condicao", projects = "all_projects"),
      img_tag(file.path(gene_dir, "expression_sample_tile.png"), "Expressao por amostra individual", projects = "all_projects"),
      img_tag(file.path(gene_dir, "deg_lollipop.png"), "Efeito DEG do gene nos contrastes disponiveis", projects = "all_projects"),
      img_tag(file.path(gene_dir, "deg_scatter.png"), "log2FC versus -log10(padj) nos contrastes DEG", projects = "all_projects")
    )
    gene_project_sections <- paste(vapply(gene_projects[gene_projects != "all_projects" & gene_projects != "unknown"], function(project) {
      project_figures <- paste0(
        project_img_tag(file.path(gene_dir, "expression_tissue_sex_condition.png"), "Expressao por tecido, sexo, condicao, estagio e projeto", project),
        project_img_tag(file.path(gene_dir, "expression_batch_project.png"), "Expressao por batch/projeto", project),
        project_img_tag(file.path(gene_dir, "expression_stage_profile.png"), "Perfil medio por estagio, tecido, batch e condicao", project),
        project_img_tag(file.path(gene_dir, "expression_sample_tile.png"), "Expressao por amostra individual", project),
        project_img_tag(file.path(gene_dir, "deg_lollipop.png"), "Efeito DEG do gene nos contrastes disponiveis", project),
        project_img_tag(file.path(gene_dir, "deg_scatter.png"), "log2FC versus -log10(padj) nos contrastes DEG", project)
      )
      project_count <- if (project_figures == "") 0L else NA_integer_
      project_body <- if (project_figures == "") "" else paste0("<div class='figure-grid'>", project_figures, "</div>")
      figure_details(paste("Projeto", project, "raw - gene"), project_body, count = project_count, open = FALSE, note = paste("Visualizacoes filtradas para", project, "; nao incluem amostras de outros projetos/batches."))
    }, character(1)), collapse = "")
    gene_corrected_figures <- paste0(
      corrected_img_tag(file.path(gene_dir, "expression_tissue_sex_condition.png"), "Expressao por tecido, sexo, condicao, estagio e projeto"),
      corrected_img_tag(file.path(gene_dir, "expression_batch_project.png"), "Expressao por batch/projeto"),
      corrected_img_tag(file.path(gene_dir, "expression_stage_profile.png"), "Perfil medio por estagio, tecido, batch e condicao"),
      corrected_img_tag(file.path(gene_dir, "expression_sample_tile.png"), "Expressao por amostra individual")
    )
    paste0(
      "<section class='searchable gene' data-kind='gene' data-projects='", project_attr(gene_projects), "' data-search='", html_escape(tolower(gene_search)), "' id='gene_", sanitize(row$group), "_", sanitize(row$matched_gene_id), "'>",
      "<h3>", html_escape(row$gene_display_label), "</h3>",
      "<div class='gene-meta'>",
      "<span><b>Grupo</b>", html_escape(row$group), "</span>",
      "<span><b>Query</b>", html_escape(row$query), "</span>",
      "<span><b>Biotipo</b>", html_escape(row$biotype), "</span>",
      "<span><b>Localizacao</b>", html_escape(row$location), "</span>",
      "<span><b>Na matriz ", html_escape(expression_unit), "</b>", html_escape(row$found_in_expression_matrix), "</span>",
      "</div>",
      ifelse(row$description != "", paste0("<p class='description'>", html_escape(row$description), "</p>"), ""),
      figure_details("All projects raw - gene", paste0("<div class='figure-grid'>", gene_raw_figures, "</div>"), open = TRUE, note = "Visualizacoes do gene na matriz principal."),
      gene_project_sections,
      figure_details("All projects corrected 055 - gene", paste0("<div class='figure-grid'>", gene_corrected_figures, "</div>"), open = FALSE, note = "Visualizacoes corrigidas ficam separadas da matriz raw."),
      details_block("Especificidade por estagio", stage_table, 50, nrow(stage_table) > 0),
      details_block("DEG do gene", deg_table, 50, nrow(deg_table) > 0),
      details_block("DTU do gene: gene/transcrito diferencial", dtu_table, 50, nrow(dtu_table) > 0),
      details_block("Splicing do gene", splicing_table, 50, nrow(splicing_table) > 0),
      details_block("WGCNA do gene", wgcna_table, 50, nrow(wgcna_table) > 0),
      details_block("Mfuzz do gene", mfuzz_table, 50, nrow(mfuzz_table) > 0),
      "</section>"
    )
  }, character(1)), collapse = "\n")

  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    paste0("<title>", html_escape(title), "</title>"),
    "<style>
      :root{--ink:#202833;--muted:#607080;--line:#d8e1ea;--soft:#f6f8fb;--accent:#255f85;--accent2:#8a5a27}
      *{box-sizing:border-box}
      html{scroll-behavior:smooth}
      body{font-family:Arial,sans-serif;max-width:1440px;margin:0 auto;padding:28px 24px 48px;line-height:1.45;color:var(--ink);background:#fff}
      nav{position:sticky;top:0;background:#fff;border-bottom:1px solid var(--line);padding:10px 0;margin:18px 0 18px;z-index:3;display:flex;gap:16px;flex-wrap:wrap}
      nav a{color:var(--accent);text-decoration:none;font-weight:700}
      h1{font-size:30px;line-height:1.15;margin:0 0 8px;color:#17324d}
      h2{color:#17324d;border-top:2px solid #e6edf3;padding-top:24px;margin-top:42px}
      h3{margin:22px 0 10px;color:#17324d}
      h4{margin:18px 0 8px;color:#17324d}
      p{max-width:980px}
      .cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:18px 0}
      .card{background:var(--soft);border:1px solid var(--line);border-radius:6px;padding:14px}
      .card .num{font-size:28px;font-weight:700;color:#17324d}
      .card .num.text{font-size:18px;line-height:1.2;overflow-wrap:anywhere}
      .toolbar{position:sticky;top:46px;background:#fff;border:1px solid #d8e1ea;border-radius:6px;padding:12px;margin:14px 0 24px 0;z-index:2;box-shadow:0 2px 10px rgba(20,45,70,.06)}
      .toolbar input{box-sizing:border-box;width:100%;font-size:16px;padding:10px 12px;border:1px solid #bdc9d6;border-radius:4px}
      .filters{display:flex;flex-wrap:wrap;gap:14px;margin-top:10px;font-size:13px;color:#34495e}
      .filters label{display:inline-flex;gap:6px;align-items:center}
      .search-count{font-size:13px;color:#5d6d7e;margin-top:8px}
      .gene-index{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:8px;margin:14px 0 24px 0}
      .gene-chip{display:block;border:1px solid #d8e1ea;border-radius:6px;padding:9px 10px;text-decoration:none;color:#17324d;background:#fbfcfd}
      .gene-chip span{display:block;font-weight:700;overflow-wrap:anywhere}.gene-chip small{display:block;color:#697b8c;margin-top:2px}
      .analysis-panel{border-top:1px solid var(--line);padding-top:16px;margin-top:22px}
      .analysis-panel>p{color:var(--muted);margin-top:0}
      .notice{border-left:4px solid var(--accent2);background:#fff8ef;padding:10px 12px;margin:12px 0 18px;color:#3f3428}
      .tabs{border:1px solid var(--line);border-radius:6px;margin:14px 0 24px;background:#fff}
      .tab-list{display:flex;flex-wrap:wrap;gap:6px;padding:8px;border-bottom:1px solid var(--line);background:#f8fafc}
      .tab-button{display:inline-flex;align-items:center;gap:8px;border:1px solid var(--line);background:#fff;color:#17324d;border-radius:6px;padding:8px 10px;font-weight:700;cursor:pointer}
      .tab-button small{min-width:24px;text-align:center;border-radius:999px;background:#e7eef5;color:#17324d;padding:1px 7px;font-size:11px}
      .tab-button.active{border-color:var(--accent);background:#eaf3f8;color:#12344d}
      .tab-panel{display:none;padding:12px}
      .tab-panel.active{display:block}
      .figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:18px;align-items:start;margin:12px 0 24px}
      .figure-grid figure{margin:0}
      .table-wrap{overflow-x:auto}
      table{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0 22px 0}
      th,td{border:1px solid #ddd;padding:5px;vertical-align:top} th{background:#f3f3f3}
      .table-note{font-size:12px;color:var(--muted);margin-top:-14px}
      code{background:#f6f6f6;padding:2px 4px}
      figure{margin:18px 0 30px 0} figcaption{font-size:13px;color:#555;margin-top:6px}
      figcaption strong{display:block;color:#2a3f53} figcaption span{display:block;margin-top:3px}
      img{display:block;width:100%;max-width:100%;border:1px solid #ddd;margin:8px 0 10px 0;background:#fff}
      .gene,.group-section{border-top:1px solid #e6e6e6;padding-top:12px}
      .gene-meta{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0 10px}
      .gene-meta span{border:1px solid var(--line);background:#fbfcfd;border-radius:6px;padding:6px 8px;font-size:12px;color:#263442}
      .gene-meta b{display:block;font-size:11px;color:var(--muted);font-weight:700}
      .description{color:var(--muted)}
      details.data-block{border:1px solid var(--line);border-radius:6px;margin:10px 0 14px;background:#fff}
      details.data-block summary{cursor:pointer;display:flex;justify-content:space-between;gap:12px;align-items:center;padding:10px 12px;font-weight:700;color:#17324d;background:#f8fafc}
      details.data-block summary span{min-width:28px;text-align:center;border-radius:999px;background:#e7eef5;color:#17324d;padding:2px 8px;font-size:12px}
      details.data-block .table-wrap{padding:0 10px}
      .hidden-by-search{display:none!important}
      ul{columns:2}
      @media(max-width:760px){
        body{padding:20px 12px 36px}
        .toolbar{top:82px}
        .figure-grid{grid-template-columns:1fr}
        ul{columns:1}
      }
    </style>",
    "<script>
      document.addEventListener('DOMContentLoaded', function(){
        var input = document.getElementById('geneSearch');
        var projectSelect = document.getElementById('projectFilter');
        var count = document.getElementById('searchCount');
        var toggles = Array.prototype.slice.call(document.querySelectorAll('[data-filter-kind]'));
        var items = Array.prototype.slice.call(document.querySelectorAll('.searchable'));
        function activeKinds(){
          return toggles.filter(function(t){return t.checked;}).map(function(t){return t.getAttribute('data-filter-kind');});
        }
        function applySearch(){
          var query = (input.value || '').trim().toLowerCase();
          var project = projectSelect ? projectSelect.value : 'all';
          var kinds = activeKinds();
          var visible = 0;
          items.forEach(function(el){
            var kind = el.getAttribute('data-kind') || '';
            var text = el.getAttribute('data-search') || el.textContent.toLowerCase();
            var projects = el.getAttribute('data-projects') || '';
            var kindOk = kind === '' || kinds.indexOf(kind) !== -1;
            var queryOk = query === '' || text.indexOf(query) !== -1;
            var projectOk = project === 'all' || projects === '' || projects.split('|').indexOf(project) !== -1;
            var show = kindOk && queryOk && projectOk;
            el.classList.toggle('hidden-by-search', !show);
            if (show && kind !== '') visible += 1;
          });
          items.forEach(function(el){
            if (!el.classList.contains('hidden-by-search')) {
              var parent = el.closest('.group-section.hidden-by-search,.gene.hidden-by-search');
              if (parent) parent.classList.remove('hidden-by-search');
            }
          });
          count.textContent = query === '' && project === 'all' ? 'Filtro inativo.' : visible + ' itens encontrados.';
        }
        input.addEventListener('input', applySearch);
        if (projectSelect) projectSelect.addEventListener('change', applySearch);
        toggles.forEach(function(t){t.addEventListener('change', applySearch);});
        document.querySelectorAll('.tab-button').forEach(function(btn){
          btn.addEventListener('click', function(){
            var group = btn.getAttribute('data-tab-group');
            var target = btn.getAttribute('data-tab-target');
            document.querySelectorAll('.tab-button[data-tab-group=\"' + group + '\"]').forEach(function(other){
              other.classList.toggle('active', other === btn);
            });
            document.querySelectorAll('.tab-panel[data-tab-panel=\"' + group + '\"]').forEach(function(panel){
              panel.classList.toggle('active', panel.id === target);
            });
          });
        });
        applySearch();
      });
    </script>",
    "</head><body>",
    paste0("<h1>", html_escape(title), "</h1>"),
    "<nav><a href='#overview'>Resumo</a><a href='#evidence'>Evidencias</a><a href='#groups'>Grupos</a><a href='#genes'>Genes</a><a href='#tables'>Tabelas</a></nav>",
    "<div class='toolbar' role='search'>",
    "<input id='geneSearch' type='search' placeholder='Buscar por gene, ID, grupo, biotipo, descricao ou contraste'>",
    "<div class='filters'>",
    "<label>Projeto <select id='projectFilter'><option value='all'>Todos</option>", project_options, "</select></label>",
    "<label><input type='checkbox' data-filter-kind='gene' checked>Genes</label>",
    "<label><input type='checkbox' data-filter-kind='group' checked>Grupos</label>",
    "<label><input type='checkbox' data-filter-kind='figure' checked>Figuras</label>",
    "<label><input type='checkbox' data-filter-kind='table' checked>Tabelas</label>",
    "</div>",
    "<div class='search-count' id='searchCount'>Filtro inativo.</div>",
    "</div>",
    "<section id='overview'>",
    "<div class='cards'>",
    paste0("<div class='card'><div class='num text'>", html_escape(expression_unit), "</div><div>matriz de expressao</div></div>"),
    paste0("<div class='card'><div class='num'>", nrow(catalog), "</div><div>entradas no genes.txt</div></div>"),
    paste0("<div class='card'><div class='num'>", n_found_genes, "</div><div>genes encontrados</div></div>"),
    paste0("<div class='card'><div class='num'>", n_annotated_genes, "</div><div>genes anotados</div></div>"),
    paste0("<div class='card'><div class='num'>", length(unique(catalog$group)), "</div><div>grupos</div></div>"),
    paste0("<div class='card'><div class='num'>", n_deg_sig_genes, "</div><div>genes com DEG significativo</div></div>"),
    paste0("<div class='card'><div class='num'>", n_dtu_sig_genes, "</div><div>genes com DTU significativo</div></div>"),
    paste0("<div class='card'><div class='num'>", n_dtu_sig_transcripts, "</div><div>transcritos com DTU significativo</div></div>"),
    paste0("<div class='card'><div class='num'>", n_dtu_report_genes, "</div><div>genes do relatorio com DTU</div></div>"),
    paste0("<div class='card'><div class='num'>", n_splicing_sig_genes, "</div><div>genes com splicing significativo</div></div>"),
    paste0("<div class='card'><div class='num'>", n_wgcna_genes, "</div><div>genes em modulos WGCNA</div></div>"),
    paste0("<div class='card'><div class='num'>", n_wgcna_hub_genes, "</div><div>genes hub WGCNA</div></div>"),
    paste0("<div class='card'><div class='num'>", n_mfuzz_genes, "</div><div>genes em clusters Mfuzz</div></div>"),
    paste0("<div class='card'><div class='num'>", n_stage_specific_genes, "</div><div>genes com evidencia estagio-especifica</div></div>"),
    paste0("<div class='card'><div class='num'>", n_corrected_figures, "</div><div>figuras corrigidas 055</div></div>"),
    paste0("<div class='card'><div class='num text'>", html_escape(generated_at), "</div><div>gerado em</div></div>"),
    "</div>",
    "<div class='gene-index'>", gene_index, "</div>",
    ifelse(corrected_expression_file != "", paste0("<p class='notice'>Figuras com sufixo <b>corrigido 055</b> usam <code>", html_escape(corrected_expression_file), "</code> como camada visual exploratoria. Tabelas DEG/DTU/WGCNA/Mfuzz permanecem ligadas aos respectivos resultados estatisticos.</p>"), ""),
    "<h2>Visualizacoes por escopo</h2>",
    scope_tabs_html(global_plots, dataset_projects, n_corrected_global_figures),
    figure_details(
      "DEG all_projects raw",
      deg_figure_grid(global_plots),
      count_existing_plots(global_plots[c("deg_heatmap", "deg_tile", "deg_direction")]),
      open = FALSE,
      note = "Resultados estatisticos de DEG permanecem separados da camada visual corrigida."
    ),
    "</section>",
    "<section id='evidence'><h2>Evidencias integradas</h2>",
    "<section class='analysis-panel searchable' data-kind='table' data-search='stage estagio especificidade especifica tau ciclo vida life cycle deg mfuzz'>",
    "<h3>Especificidade por estagio</h3>",
    paste0("<p>Resume genes candidatos com evidencia por expressao concentrada em estagio (tau >= ", stage_tau_threshold, ", expressao media maxima >= ", stage_min_expression, "), DEG significativo em contrastes de <code>stage</code>, ou membership Mfuzz >= ", mfuzz_membership_threshold, ".</p>"),
    table_to_html(stage_specificity_table, 100),
    "</section>",
    "<section class='analysis-panel searchable' data-kind='figure' data-search='batch correction correct corrigido combat pycombat pca antes depois'>",
    "<h3>Batch correction 055: antes e depois</h3>",
    "<p>Inclui figuras nativas da etapa 055 e indica se o relatorio encontrou uma matriz corrigida para gerar graficos exploratorios pareados.</p>",
    paste0("<p><b>Matriz corrigida usada nas figuras:</b> ", ifelse(corrected_expression_file == "", "<em>nao detectada</em>", paste0("<code>", html_escape(corrected_expression_file), "</code>")), "</p>"),
    figure_gallery_html(batch_plot_gallery),
    "</section>",
    "<section class='analysis-panel searchable' data-kind='table' data-search='dtu differential transcript usage gene transcript transcrito uso diferencial de transcritos'>",
    "<h3>DTU: genes e transcritos diferenciais</h3>",
    "<p>Tabela direta para localizar quais pares gene/transcrito tiveram uso diferencial significativo e em qual projeto/variavel apareceram.</p>",
    table_to_html(dtu_overview_table, 100),
    "</section>",
    "<section class='analysis-panel searchable' data-kind='figure' data-search='wgcna weighted gene coexpression network analysis modulo hub soft threshold module trait'>",
    "<h3>WGCNA: modulos, hubs e figuras</h3>",
    "<p>Integra os genes do relatorio aos modulos WGCNA e exibe as figuras geradas pela etapa 085 quando elas existem no diretorio de resultados.</p>",
    table_to_html(wgcna_overview_table, 100),
    figure_gallery_html(wgcna_plot_gallery),
    "</section>",
    "<section class='analysis-panel searchable' data-kind='figure' data-search='mfuzz fuzzy clustering cluster temporal membership centroides'>",
    "<h3>Mfuzz: clusters temporais e figuras</h3>",
    "<p>Mostra a associacao dos genes aos clusters Mfuzz e inclui os centroides dos clusters quando a etapa 086 gerou imagens.</p>",
    table_to_html(mfuzz_overview_table, 100),
    figure_gallery_html(mfuzz_plot_gallery),
    "</section>",
    "</section>",
    "<section id='groups'><h2>Grupos</h2><ul>", group_links, "</ul>", group_sections, "</section>",
    "<section id='genes'><h2>Genes individuais</h2>", gene_sections, "</section>",
    "<section id='tables'><h2>Tabelas</h2>",
    paste0("<p>Arquivos completos: <code>", html_escape(add_table_suffix("tables/gene_catalog.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/gene_expression_summary.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/stage_specificity_summary.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/expression_long.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/expression_summary_by_context.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/deg_hits.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/dtu_hits.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/splicing_hits.tsv")), "</code>, <code>", html_escape(add_table_suffix("tables/wgcna_hits.tsv")), "</code> e <code>", html_escape(add_table_suffix("tables/mfuzz_hits.tsv")), "</code>.</p>"),
    table_to_html(gene_summary, 100),
    "</section>",
    "</body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
}

gene_groups <- parse_gene_groups(genes_file)
tpm <- read_matrix(tpm_file)
samples <- read_samples(samples_file, setdiff(colnames(tpm), "gene_id"))

if (metadata_file != "" && file.exists(metadata_file)) {
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (all(c("dataset", "sample_id") %in% colnames(metadata))) {
    metadata$import_id_combined <- paste(metadata$dataset, metadata$sample_id, sep = "__")
    key <- if (all(samples$import_id %in% metadata$import_id_combined)) "import_id_combined" else "sample_id"
    extra <- metadata[match(samples$import_id, metadata[[key]]), , drop = FALSE]
    add_cols <- setdiff(colnames(extra), colnames(samples))
    samples <- dplyr::bind_cols(samples, extra[, add_cols, drop = FALSE])
  }
}

samples <- complete_sample_fields(samples)
annotations <- load_annotations(gff_file)
gene_catalog <- build_gene_catalog(gene_groups, tpm, annotations)
expr_long <- make_expression_long(tpm, samples, gene_catalog)
expr_summary <- summarise_expression(expr_long)
corrected_expr_long <- NULL
corrected_expr_summary <- NULL
if (corrected_expression_file != "" && file.exists(corrected_expression_file)) {
  corrected_matrix <- read_matrix(corrected_expression_file)
  corrected_sample_cols <- intersect(setdiff(colnames(corrected_matrix), "gene_id"), samples$import_id)
  if (length(corrected_sample_cols) < 2) {
    warning("[WARN] Matriz corrigida encontrada, mas menos de duas amostras batem com a metadata do relatorio: ", corrected_expression_file)
  } else {
    corrected_matrix <- corrected_matrix[, c("gene_id", corrected_sample_cols), drop = FALSE]
    corrected_expr_long <- make_expression_long(corrected_matrix, samples, gene_catalog, expression_source = "batch_corrected_055")
    corrected_expr_summary <- summarise_expression(corrected_expr_long)
  }
}
deg_hits <- load_deg_hits(deg_root, gene_catalog)
deg_hits_annotated <- annotate_deg_hits(deg_hits, gene_catalog)
dtu_hits <- load_dtu_hits(dtu_root, gene_catalog)
dtu_hits_annotated <- annotate_dtu_hits(dtu_hits, gene_catalog)
splicing_hits <- load_splicing_hits(splicing_root, gene_catalog)
splicing_hits_annotated <- annotate_splicing_hits(splicing_hits, gene_catalog)
wgcna_hits <- load_wgcna_hits(wgcna_root, gene_catalog)
wgcna_hits_annotated <- annotate_wgcna_hits(wgcna_hits, gene_catalog)
mfuzz_hits <- load_mfuzz_hits(mfuzz_root, gene_catalog)
mfuzz_hits_annotated <- annotate_mfuzz_hits(mfuzz_hits, gene_catalog)
stage_specificity <- summarise_stage_specificity(expr_long, deg_hits, mfuzz_hits)
gene_summary <- summarise_gene_descriptives(expr_long, expr_summary, deg_hits, gene_catalog)
gene_summary <- add_stage_specificity_summary(gene_summary, stage_specificity)
gene_summary <- add_optional_gene_summaries(gene_summary, dtu_hits, splicing_hits, wgcna_hits, mfuzz_hits)

write_tsv2(gene_catalog, file.path(out_dir, "tables", "gene_catalog.tsv"))
write_tsv2(expr_long, file.path(out_dir, "tables", "expression_long.tsv"))
write_tsv2(expr_summary, file.path(out_dir, "tables", "expression_summary_by_context.tsv"))
if (!is.null(corrected_expr_long) && !is.null(corrected_expr_summary)) {
  write_tsv2(corrected_expr_long, file.path(out_dir, "tables", "expression_long_batch_corrected.tsv"))
  write_tsv2(corrected_expr_summary, file.path(out_dir, "tables", "expression_summary_by_context_batch_corrected.tsv"))
}
write_tsv2(deg_hits_annotated, file.path(out_dir, "tables", "deg_hits.tsv"))
write_tsv2(dtu_hits_annotated, file.path(out_dir, "tables", "dtu_hits.tsv"))
write_tsv2(splicing_hits_annotated, file.path(out_dir, "tables", "splicing_hits.tsv"))
write_tsv2(wgcna_hits_annotated, file.path(out_dir, "tables", "wgcna_hits.tsv"))
write_tsv2(mfuzz_hits_annotated, file.path(out_dir, "tables", "mfuzz_hits.tsv"))
write_tsv2(stage_specificity, file.path(out_dir, "tables", "stage_specificity_summary.tsv"))
write_tsv2(gene_summary, file.path(out_dir, "tables", "gene_expression_summary.tsv"))

global_plots <- list(
  heatmap = file.path("plots", "all_groups_expression_heatmap.png"),
  dotplot = file.path("plots", "all_groups_expression_dotplot.png"),
  annotated_sample_heatmap = file.path("plots", "all_groups_sample_heatmap_annotated.png"),
  gene_correlation = file.path("plots", "all_groups_gene_correlation.png"),
  sample_pca = file.path("plots", "all_groups_sample_pca.png"),
  sample_mds = file.path("plots", "all_groups_sample_mds.png"),
  tissue_sex_heatmap = file.path("plots", "all_groups_tissue_sex_heatmap.png"),
  ovary_testis = file.path("plots", "all_groups_ovary_testis_panel.png"),
  group_aggregate = file.path("plots", "all_groups_aggregate_profile.png"),
  batch_project = file.path("plots", "all_groups_batch_project_boxplot.png"),
  deg_heatmap = file.path("plots", "all_groups_deg_log2fc_heatmap.png"),
  deg_tile = file.path("plots", "all_groups_deg_context_tile.png"),
  deg_direction = file.path("plots", "all_groups_deg_direction_summary.png")
)

corrected_global_plots <- as.list(vapply(global_plots, corrected_plot_path, character(1)))

invisible(plot_or_skip("global expression heatmap", function() plot_expression_heatmap(expr_summary, file.path(out_dir, global_plots$heatmap), "Todos os grupos - expressao media")))
invisible(plot_or_skip("global expression dotplot", function() plot_expression_dotplot(expr_summary, file.path(out_dir, global_plots$dotplot), "Todos os grupos - expressao media e fracao expressa")))
invisible(plot_or_skip("global annotated sample heatmap", function() plot_annotated_sample_heatmap(expr_long, file.path(out_dir, global_plots$annotated_sample_heatmap), "Todos os grupos - amostras anotadas")))
invisible(plot_or_skip("global gene correlation", function() plot_gene_correlation(expr_long, file.path(out_dir, global_plots$gene_correlation), "Todos os grupos - correlacao entre genes")))
invisible(plot_or_skip("global sample PCA", function() plot_sample_ordination(expr_long, file.path(out_dir, global_plots$sample_pca), method = "pca", title = "Todos os grupos - PCA das amostras")))
invisible(plot_or_skip("global sample MDS", function() plot_sample_ordination(expr_long, file.path(out_dir, global_plots$sample_mds), method = "mds", title = "Todos os grupos - MDS das amostras")))
invisible(plot_or_skip("global tissue/sex heatmap", function() plot_tissue_sex_heatmap(expr_summary, file.path(out_dir, global_plots$tissue_sex_heatmap), "Todos os grupos - tecido e sexo")))
invisible(plot_or_skip("global ovary/testis", function() plot_ovary_testis_panel(expr_summary, file.path(out_dir, global_plots$ovary_testis), "Todos os grupos - ovario versus testiculo")))
invisible(plot_or_skip("global group aggregate profile", function() plot_group_aggregate_profile(expr_long, file.path(out_dir, global_plots$group_aggregate), "Todos os grupos - perfil agregado")))
invisible(plot_or_skip("global batch/project", function() plot_batch_project_boxplot(expr_long, file.path(out_dir, global_plots$batch_project), "Todos os grupos - batch/projeto")))

project_names <- sort(unique(as.character(expr_long$dataset)))
project_names <- project_names[project_names != "" & project_names != "unknown" & !is.na(project_names)]
for (project in project_names) {
  project_expr_long <- expr_long %>% dplyr::filter(dataset == project)
  project_expr_summary <- expr_summary %>% dplyr::filter(dataset == project)
  project_deg_hits <- deg_hits %>% dplyr::filter(deg_project == project)
  plot_expression_scope_outputs(project_expr_long, project_expr_summary, project_plot_paths(project), paste("Projeto", project, "raw"))
  plot_group_outputs(
    project_expr_long,
    project_expr_summary,
    project_deg_hits,
    gene_catalog,
    out_dir,
    suffix = project_plot_suffix(project),
    title_suffix = paste("projeto", project),
    include_deg = nrow(project_deg_hits) > 0
  )
  plot_gene_outputs(
    project_expr_long,
    project_deg_hits,
    out_dir,
    suffix = project_plot_suffix(project),
    title_suffix = paste("projeto", project),
    include_deg = nrow(project_deg_hits) > 0
  )
}

invisible(plot_or_skip("global DEG heatmap", function() plot_deg_heatmap(deg_hits, gene_catalog, file.path(out_dir, global_plots$deg_heatmap), "Todos os grupos - log2FC DEG")))
invisible(plot_or_skip("global DEG tile", function() plot_deg_context_tile(deg_hits, gene_catalog, file.path(out_dir, global_plots$deg_tile), "Todos os grupos - DEG por contraste")))
invisible(plot_or_skip("global DEG direction", function() plot_deg_direction_summary(deg_hits, gene_catalog, file.path(out_dir, global_plots$deg_direction), "Todos os grupos - direcao DEG")))
plot_group_outputs(expr_long, expr_summary, deg_hits, gene_catalog, out_dir)
plot_gene_outputs(expr_long, deg_hits, out_dir)

if (!is.null(corrected_expr_long) && !is.null(corrected_expr_summary)) {
  with_expression_unit(corrected_expression_unit, {
    invisible(plot_or_skip("global corrected expression heatmap", function() plot_expression_heatmap(corrected_expr_summary, file.path(out_dir, corrected_global_plots$heatmap), "Todos os grupos - expressao media corrigida 055")))
    invisible(plot_or_skip("global corrected expression dotplot", function() plot_expression_dotplot(corrected_expr_summary, file.path(out_dir, corrected_global_plots$dotplot), "Todos os grupos - expressao corrigida 055 e fracao expressa")))
    invisible(plot_or_skip("global corrected annotated sample heatmap", function() plot_annotated_sample_heatmap(corrected_expr_long, file.path(out_dir, corrected_global_plots$annotated_sample_heatmap), "Todos os grupos - amostras anotadas corrigidas 055")))
    invisible(plot_or_skip("global corrected gene correlation", function() plot_gene_correlation(corrected_expr_long, file.path(out_dir, corrected_global_plots$gene_correlation), "Todos os grupos - correlacao corrigida 055")))
    invisible(plot_or_skip("global corrected sample PCA", function() plot_sample_ordination(corrected_expr_long, file.path(out_dir, corrected_global_plots$sample_pca), method = "pca", title = "Todos os grupos - PCA corrigida 055")))
    invisible(plot_or_skip("global corrected sample MDS", function() plot_sample_ordination(corrected_expr_long, file.path(out_dir, corrected_global_plots$sample_mds), method = "mds", title = "Todos os grupos - MDS corrigida 055")))
    invisible(plot_or_skip("global corrected tissue/sex heatmap", function() plot_tissue_sex_heatmap(corrected_expr_summary, file.path(out_dir, corrected_global_plots$tissue_sex_heatmap), "Todos os grupos - tecido e sexo corrigido 055")))
    invisible(plot_or_skip("global corrected ovary/testis", function() plot_ovary_testis_panel(corrected_expr_summary, file.path(out_dir, corrected_global_plots$ovary_testis), "Todos os grupos - ovario versus testiculo corrigido 055")))
    invisible(plot_or_skip("global corrected group aggregate profile", function() plot_group_aggregate_profile(corrected_expr_long, file.path(out_dir, corrected_global_plots$group_aggregate), "Todos os grupos - perfil agregado corrigido 055")))
    invisible(plot_or_skip("global corrected batch/project", function() plot_batch_project_boxplot(corrected_expr_long, file.path(out_dir, corrected_global_plots$batch_project), "Todos os grupos - batch/projeto corrigido 055")))
    plot_group_outputs(corrected_expr_long, corrected_expr_summary, deg_hits, gene_catalog, out_dir, suffix = "_corrected", title_suffix = "corrigido 055", include_deg = FALSE)
    plot_gene_outputs(corrected_expr_long, deg_hits, out_dir, suffix = "_corrected", title_suffix = "corrigido 055", include_deg = FALSE)
  })
} else {
  stale_corrected <- list.files(out_dir, pattern = "_corrected\\.png$", recursive = TRUE, full.names = TRUE)
  if (length(stale_corrected) > 0) unlink(stale_corrected)
}

write_html_report(file.path(out_dir, "gene_set_report.html"), report_title, gene_catalog, gene_summary, stage_specificity, deg_hits_annotated, dtu_hits_annotated, splicing_hits_annotated, wgcna_hits_annotated, mfuzz_hits_annotated, global_plots, expression_unit, corrected_expression_file, corrected_expression_unit)
log_info(paste("[OK] Relatorio 090 concluido:", file.path(out_dir, "gene_set_report.html")))
