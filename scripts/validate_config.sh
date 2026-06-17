#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE="${1:-config/pipeline_config.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
load_rnaseq_config "${CONFIG_FILE}"

create_rnaseq_output_tree

[[ -d "${PROJECT_DIR}" ]] || die "PROJECT_DIR does not exist: ${PROJECT_DIR}"

if [[ ! -f "${USER_SETTINGS_FILE:-}" ]]; then
  warn "User settings file not found: ${USER_SETTINGS_FILE:-config/user_settings.sh}"
  warn "For simple setup, run: cp config/user_settings_template.sh config/user_settings.sh"
fi

if [[ -z "${PIPELINE_PROJECTS//,/ }" ]]; then
  die "PIPELINE_PROJECTS is empty. Add at least one project accession."
fi

case "${PIPELINE_PROJECTS}" in
  *PRJXXXX*|*PRJYYYY*)
    die "PIPELINE_PROJECTS still contains template values. Edit config/user_settings.sh."
    ;;
esac

if [[ "${ORGANISM_NAME}" == "My organism" || "${ORGANISM_NAME}" == "custom_organism" ]]; then
  die "ORGANISM_NAME still contains a template value. Edit config/user_settings.sh."
fi

for value in "${SCRATCH_ROOT:-}" "${CONDA_BASE:-}" "${GENOME_URL:-}" "${TRANSCRIPTS_URL:-}" "${GFF3_URL:-}" "${REF_GENOME_FA:-}" "${REF_TRANSCRIPTS_FA:-}" "${REF_GFF3:-}" "${REF_GTF:-}"; do
  case "${value}" in
    *example.org*|/path/to/*|*/path/to/*|*/my_user/*)
      die "Configuration still contains a template path or URL: ${value}"
      ;;
  esac
done

if [[ -z "${CONDA_BASE:-}" ]]; then
  warn "CONDA_BASE is empty. Jobs will fail unless conda is available in the Slurm environment."
fi

case "${QUANT_METHOD:-salmon}" in
  salmon|star)
    ;;
  *)
    die "QUANT_METHOD must be 'salmon' or 'star'. Current value: ${QUANT_METHOD}"
    ;;
esac

case "${PIPELINE_COMPRESS_RESULTS:-1}" in
  0|1|true|TRUE|yes|YES|false|FALSE|no|NO|y|Y|n|N)
    ;;
  *)
    die "PIPELINE_COMPRESS_RESULTS must be 0/1 or yes/no. Current value: ${PIPELINE_COMPRESS_RESULTS}"
    ;;
esac

case "${STAR_GENECOUNT_COLUMN:-unstranded}" in
  unstranded|stranded_forward|stranded_reverse|2|3|4)
    ;;
  *)
    die "STAR_GENECOUNT_COLUMN must be unstranded, stranded_forward, stranded_reverse, 2, 3, or 4."
    ;;
esac

if ! [[ "${STAR_LIMIT_BAM_SORT_RAM:-}" =~ ^[0-9]+$ ]] || [[ "${STAR_LIMIT_BAM_SORT_RAM}" -lt 1 ]]; then
  die "STAR_LIMIT_BAM_SORT_RAM must be a positive integer in bytes. Current value: ${STAR_LIMIT_BAM_SORT_RAM:-unset}"
fi

if [[ "${RUN_SALMON_INDEX}" == "1" && -z "${REF_TRANSCRIPTS_FA}" && -z "${TRANSCRIPTS_URL}" ]]; then
  die "Salmon index requested, but REF_TRANSCRIPTS_FA and TRANSCRIPTS_URL are empty."
fi

if [[ "${RUN_STAR_INDEX}" == "1" && -z "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
  die "STAR index requested, but REF_GENOME_FA and GENOME_URL are empty."
fi

if [[ "${RUN_STAR_GTF_INDEX}" == "1" ]]; then
  if [[ -z "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
    die "STAR+GTF index requested, but REF_GENOME_FA and GENOME_URL are empty."
  fi
  if [[ -z "${REF_GTF}" && -z "${GTF_URL}" && -z "${REF_GFF3}" && -z "${GFF3_URL}" ]]; then
    die "STAR+GTF index requested, but no GTF/GFF3 input was configured."
  fi
fi

if [[ "${QUANT_METHOD}" == "star" ]]; then
  if [[ "${RUN_STAR_GTF_INDEX}" != "1" && ! -d "${STAR_QUANT_INDEX_DIR}" ]]; then
    warn "QUANT_METHOD=star but STAR_QUANT_INDEX_DIR does not exist yet: ${STAR_QUANT_INDEX_DIR}"
    warn "Set RUN_STAR_GTF_INDEX=1 to build it, or point STAR_QUANT_INDEX_DIR to an existing STAR index built with a GTF/GFF annotation."
  fi
  if [[ "${RUN_STAR_GTF_INDEX}" == "1" && -z "${REF_GTF}" && -z "${GTF_URL}" && -z "${REF_GFF3}" && -z "${GFF3_URL}" ]]; then
    die "QUANT_METHOD=star needs REF_GTF/GTF_URL or REF_GFF3/GFF3_URL so STAR can produce gene counts."
  fi
fi

case "${PIPELINE_EXECUTOR:-slurm}" in
  slurm|local)
    ;;
  *)
    die "PIPELINE_EXECUTOR must be 'slurm' or 'local'. Current value: ${PIPELINE_EXECUTOR}"
    ;;
esac

case "${PIPELINE_STORAGE_MODE:-full}" in
  full|balanced|minimal)
    ;;
  *)
    die "PIPELINE_STORAGE_MODE must be 'full', 'balanced', or 'minimal'. Current value: ${PIPELINE_STORAGE_MODE}"
    ;;
esac

for flag in \
  RUN_STORAGE_CLEANUP_AFTER_ALIGNMENT \
  CLEANUP_FASTQC_DIRS \
  CLEANUP_TRIMMED_RUNS \
  CLEANUP_TRIMMED_MERGED \
  CLEANUP_FASTQ_FTP \
  CLEANUP_STAR_BAM \
  RUN_DTU_ANALYSIS \
  RUN_SPLICING_ANALYSIS
do
  case "${!flag:-0}" in
    0|1|true|TRUE|yes|YES|false|FALSE|no|NO|y|Y|n|N)
      ;;
    *)
      die "${flag} must be 0/1 or yes/no. Current value: ${!flag}"
      ;;
  esac
done

for numeric_var in \
  DTU_MIN_REPLICATES \
  DTU_MIN_GENE_COUNT \
  DTU_MIN_TRANSCRIPTS_PER_GENE \
  DTU_CONCURRENCY \
  SPLICING_MIN_REPLICATES \
  SPLICING_READ_LENGTH \
  SPLICING_CONCURRENCY
do
  numeric_value="${!numeric_var:-}"
  if ! [[ "${numeric_value}" =~ ^[0-9]+$ ]] || [[ "${numeric_value}" -lt 1 ]]; then
    die "${numeric_var} must be a positive integer. Current value: ${numeric_value:-unset}"
  fi
done
unset numeric_var numeric_value

case "${SPLICING_LIB_TYPE:-fr-unstranded}" in
  fr-unstranded|fr-firststrand|fr-secondstrand)
    ;;
  *)
    die "SPLICING_LIB_TYPE must be fr-unstranded, fr-firststrand, or fr-secondstrand. Current value: ${SPLICING_LIB_TYPE}"
    ;;
esac

if [[ "${RUN_DTU_ANALYSIS:-0}" == "1" ]]; then
  if [[ "${QUANT_METHOD:-salmon}" != "salmon" ]]; then
    warn "RUN_DTU_ANALYSIS=1 but QUANT_METHOD=${QUANT_METHOD}. DTU uses Salmon transcript-level quant.sf files; it can run only if Salmon quantifications already exist under QUANT_DIR."
  fi
  if [[ -z "${REF_TRANSCRIPTS_FA:-}" && -z "${TRANSCRIPTS_URL:-}" && ! -f "${TX2GENE_FILE:-}" ]]; then
    die "RUN_DTU_ANALYSIS=1 needs transcript annotation through REF_TRANSCRIPTS_FA/TRANSCRIPTS_URL or an existing TX2GENE_FILE."
  fi
fi

if [[ "${RUN_SPLICING_ANALYSIS:-0}" == "1" ]]; then
  if [[ "${QUANT_METHOD:-salmon}" != "star" ]]; then
    warn "RUN_SPLICING_ANALYSIS=1 but QUANT_METHOD=${QUANT_METHOD}. Splicing uses STAR sorted BAMs; it can run only if STAR BAMs already exist under STAR_QUANT_DIR."
  fi
  if [[ -z "${REF_GTF:-}" && -z "${GTF_URL:-}" ]]; then
    die "RUN_SPLICING_ANALYSIS=1 needs REF_GTF or GTF_URL because rMATS requires a GTF annotation."
  fi
  if [[ "${PIPELINE_STORAGE_MODE:-full}" == "minimal" ]]; then
    warn "RUN_SPLICING_ANALYSIS=1 with PIPELINE_STORAGE_MODE=minimal: STAR BAM cleanup is disabled by default in config, because rMATS needs BAM files."
  fi
fi

while read -r project; do
  [[ -n "${project}" ]] || continue
  if [[ ! -f "${DATASET_CONFIG_DIR}/${project}/config.yaml" ]]; then
    warn "Download config missing for ${project}: ${DATASET_CONFIG_DIR}/${project}/config.yaml"
  fi
  if [[ ! -f "${METADATA_PARSER_DIR}/${project}/configs/${project}.yaml" ]]; then
    warn "Metadata parser YAML missing for ${project}: ${METADATA_PARSER_DIR}/${project}/configs/${project}.yaml"
  fi
done < <(pipeline_projects)

if [[ "${PIPELINE_EXECUTOR:-slurm}" == "slurm" && "${SKIP_SLURM_CHECK:-false}" != "true" ]]; then
  require_cmd sbatch
fi

log "Configuration validation passed"
