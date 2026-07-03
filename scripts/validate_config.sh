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
  warn "For guided setup, run:"
  warn "  python scripts/bootstrap_project.py --project PRJNA000000 --organism 'Example organism' --scratch-root /scratch/my_user/rnaseq --conda-base /path/to/miniconda3"
  warn "For manual setup, run: cp config/user_settings_template.sh config/user_settings.sh"
fi

if [[ -z "${PIPELINE_PROJECTS//,/ }" ]]; then
  die "PIPELINE_PROJECTS is empty. Add at least one project accession, or run scripts/bootstrap_project.py."
fi

case "${PIPELINE_PROJECTS}" in
  *PRJXXXX*|*PRJYYYY*)
    die "PIPELINE_PROJECTS still contains template values. Edit config/user_settings.sh or run scripts/bootstrap_project.py."
    ;;
esac

if [[ "${ORGANISM_NAME}" == "My organism" || "${ORGANISM_NAME}" == "custom_organism" ]]; then
  die "ORGANISM_NAME still contains a template value. Edit config/user_settings.sh or run scripts/bootstrap_project.py."
fi

for value in "${SCRATCH_ROOT:-}" "${CONDA_BASE:-}" "${GENOME_URL:-}" "${TRANSCRIPTS_URL:-}" "${GFF3_URL:-}" "${REF_GENOME_FA:-}" "${REF_TRANSCRIPTS_FA:-}" "${REF_GFF3:-}" "${REF_GTF:-}"; do
  case "${value}" in
    *example.org*|/path/to/*|*/path/to/*|*/my_user/*)
      die "Configuration still contains a template path or URL: ${value}. Edit config/user_settings.sh."
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

case "${STAR_WRITE_BAM:-0}" in
  0|1|true|TRUE|yes|YES|false|FALSE|no|NO|y|Y|n|N)
    ;;
  *)
    die "STAR_WRITE_BAM must be 0/1 or yes/no. Current value: ${STAR_WRITE_BAM}"
    ;;
esac

if ! [[ "${STAR_LIMIT_BAM_SORT_RAM:-}" =~ ^[0-9]+$ ]] || [[ "${STAR_LIMIT_BAM_SORT_RAM}" -lt 1 ]]; then
  die "STAR_LIMIT_BAM_SORT_RAM must be a positive integer in bytes. Current value: ${STAR_LIMIT_BAM_SORT_RAM:-unset}"
fi

if [[ "${RUN_SALMON_INDEX}" == "1" && -z "${REF_TRANSCRIPTS_FA}" && -z "${TRANSCRIPTS_URL}" ]]; then
  die "Salmon index requested, but REF_TRANSCRIPTS_FA and TRANSCRIPTS_URL are empty."
fi
if [[ "${RUN_SALMON_INDEX}" == "1" && -n "${REF_TRANSCRIPTS_FA}" && ! -f "${REF_TRANSCRIPTS_FA}" && -z "${TRANSCRIPTS_URL}" ]]; then
  die "Salmon index requested, but REF_TRANSCRIPTS_FA does not exist and TRANSCRIPTS_URL is empty: ${REF_TRANSCRIPTS_FA}"
fi

if [[ "${RUN_STAR_INDEX}" == "1" && -z "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
  die "STAR index requested, but REF_GENOME_FA and GENOME_URL are empty."
fi
if [[ "${RUN_STAR_INDEX}" == "1" && -n "${REF_GENOME_FA}" && ! -f "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
  die "STAR index requested, but REF_GENOME_FA does not exist and GENOME_URL is empty: ${REF_GENOME_FA}"
fi

if [[ "${RUN_STAR_GTF_INDEX}" == "1" ]]; then
  if [[ -z "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
    die "STAR+GTF index requested, but REF_GENOME_FA and GENOME_URL are empty."
  fi
  if [[ -n "${REF_GENOME_FA}" && ! -f "${REF_GENOME_FA}" && -z "${GENOME_URL}" ]]; then
    die "STAR+GTF index requested, but REF_GENOME_FA does not exist and GENOME_URL is empty: ${REF_GENOME_FA}"
  fi
  if [[ -z "${REF_GTF}" && -z "${GTF_URL}" && -z "${REF_GFF3}" && -z "${GFF3_URL}" ]]; then
    die "STAR+GTF index requested, but no GTF/GFF3 input was configured."
  fi
  if [[ -n "${REF_GTF}" && ! -f "${REF_GTF}" && -z "${GTF_URL}" ]]; then
    die "STAR+GTF index requested, but REF_GTF does not exist and GTF_URL is empty: ${REF_GTF}"
  fi
  if [[ -n "${REF_GFF3}" && ! -f "${REF_GFF3}" && -z "${GFF3_URL}" ]]; then
    die "STAR+GTF index requested, but REF_GFF3 does not exist and GFF3_URL is empty: ${REF_GFF3}"
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
  RUN_SPLICING_ANALYSIS \
  RUN_WGCNA_ANALYSIS \
  RUN_MFUZZ_ANALYSIS
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
  SPLICING_CONCURRENCY \
  WGCNA_CONCURRENCY \
  WGCNA_MIN_SAMPLES \
  WGCNA_MIN_GENES \
  WGCNA_MIN_MODULE_SIZE \
  WGCNA_TOP_HUBS \
  WGCNA_MAX_BLOCK_SIZE \
  MFUZZ_CONCURRENCY \
  MFUZZ_MIN_SAMPLES \
  MFUZZ_MIN_GENES \
  MFUZZ_CLUSTERS
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

for numeric_var in \
  WGCNA_MIN_EXPRESSION \
  WGCNA_MIN_FRACTION \
  WGCNA_FIT_CUTOFF \
  WGCNA_MERGE_CUT_HEIGHT \
  MFUZZ_M \
  MFUZZ_MIN_EXPRESSION \
  MFUZZ_MIN_FRACTION
do
  numeric_value="${!numeric_var:-}"
  if ! [[ "${numeric_value}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    die "${numeric_var} must be numeric. Current value: ${numeric_value:-unset}"
  fi
done
unset numeric_var numeric_value

if ! [[ "${WGCNA_POWER:-0}" =~ ^[0-9]+$ ]]; then
  die "WGCNA_POWER must be a non-negative integer, or 0 for automatic selection. Current value: ${WGCNA_POWER:-unset}"
fi

case "${WGCNA_NETWORK_TYPE:-signed}" in
  signed|unsigned|signed-hybrid|"signed hybrid")
    ;;
  *)
    die "WGCNA_NETWORK_TYPE must be signed, unsigned, or signed-hybrid. Current value: ${WGCNA_NETWORK_TYPE}"
    ;;
esac

case "${WGCNA_COR_METHOD:-pearson}" in
  pearson|bicor)
    ;;
  *)
    die "WGCNA_COR_METHOD must be pearson or bicor. Current value: ${WGCNA_COR_METHOD}"
    ;;
esac

if [[ "${RUN_MFUZZ_ANALYSIS:-0}" == "1" && -z "${MFUZZ_TIME_VARIABLE:-}" ]]; then
  die "RUN_MFUZZ_ANALYSIS=1 needs MFUZZ_TIME_VARIABLE, for example stage."
fi

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
  case "${STAR_WRITE_BAM:-0}" in
    1|true|TRUE|yes|YES|y|Y)
      ;;
    *)
      warn "RUN_SPLICING_ANALYSIS=1 but STAR_WRITE_BAM=${STAR_WRITE_BAM:-0}. rMATS needs STAR sorted BAMs; enable STAR_WRITE_BAM=1 unless BAMs already exist under STAR_QUANT_DIR."
      ;;
  esac
  if [[ -z "${REF_GTF:-}" && -z "${GTF_URL:-}" ]]; then
    die "RUN_SPLICING_ANALYSIS=1 needs REF_GTF or GTF_URL because rMATS requires a GTF annotation."
  fi
  if [[ "${PIPELINE_STORAGE_MODE:-full}" == "minimal" ]]; then
    warn "RUN_SPLICING_ANALYSIS=1 with PIPELINE_STORAGE_MODE=minimal: STAR BAM cleanup is disabled by default in config, because rMATS needs BAM files."
  fi
fi

missing_project_setup=0
while read -r project; do
  [[ -n "${project}" ]] || continue
  if [[ ! -f "${DATASET_CONFIG_DIR}/${project}/config.yaml" ]]; then
    warn "Download config missing for ${project}: ${DATASET_CONFIG_DIR}/${project}/config.yaml"
    warn "Create starter files with: python scripts/bootstrap_project.py --project ${project} --organism '${ORGANISM_NAME}' --scratch-root '${SCRATCH_ROOT}' --conda-base '${CONDA_BASE:-/path/to/miniconda3}'"
    missing_project_setup=1
  else
    link_file="$(awk -F: '/^[[:space:]]*link_file:/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "${DATASET_CONFIG_DIR}/${project}/config.yaml")"
    if [[ -z "${link_file}" ]]; then
      warn "Download config for ${project} has no download.link_file field: ${DATASET_CONFIG_DIR}/${project}/config.yaml"
    elif [[ ! -f "${DATASET_CONFIG_DIR}/${project}/${link_file}" ]]; then
      warn "FASTQ link file missing for ${project}: ${DATASET_CONFIG_DIR}/${project}/${link_file}"
      warn "Add an ENA-generated wget script or one ftp/https FASTQ URL per line."
    elif ! grep -Eq "(ftp|https)://" "${DATASET_CONFIG_DIR}/${project}/${link_file}"; then
      warn "FASTQ link file for ${project} contains no ftp/https URLs yet: ${DATASET_CONFIG_DIR}/${project}/${link_file}"
    fi
  fi
  parser_yaml="${METADATA_PARSER_DIR}/${project}/configs/${project}.yaml"
  if [[ ! -f "${parser_yaml}" ]]; then
    warn "Metadata parser YAML missing for ${project}: ${parser_yaml}"
    warn "Create starter files with: python scripts/bootstrap_project.py --project ${project} --organism '${ORGANISM_NAME}' --scratch-root '${SCRATCH_ROOT}' --conda-base '${CONDA_BASE:-/path/to/miniconda3}'"
    missing_project_setup=1
  elif grep -Eq "PRJXXXX|SAMPLE|adult: adult|control: control" "${parser_yaml}"; then
    warn "Metadata parser YAML for ${project} still looks like a starter/template. Review regex_map and defaults: ${parser_yaml}"
  fi
done < <(pipeline_projects)
unset parser_yaml link_file

if [[ "${missing_project_setup}" -eq 1 ]]; then
  warn "Project setup is incomplete. Dry-runs may still work, but download/metadata steps will fail until the files above are created."
fi

if [[ -f "$(metadata_default)" ]]; then
  log "Metadata file detected. For a stronger metadata audit, run:"
  log "  python scripts/validate_metadata.py --metadata '$(metadata_default)' --strict"
else
  warn "Final metadata not found yet: $(metadata_default). This is expected before step 025 finishes."
fi

if [[ "${PIPELINE_EXECUTOR:-slurm}" == "slurm" && "${SKIP_SLURM_CHECK:-false}" != "true" ]]; then
  require_cmd sbatch
fi

log "Configuration validation passed"
