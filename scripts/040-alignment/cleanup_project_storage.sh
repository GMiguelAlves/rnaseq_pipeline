#!/bin/bash
#SBATCH --job-name=cleanup_storage
#SBATCH --output=logs/cleanup/cleanup_%j.out
#SBATCH --error=logs/cleanup/cleanup_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00

set -euo pipefail

usage() {
    cat <<'USAGE'
Uso: cleanup_project_storage.sh <PROJECT> [opcoes]

Remove arquivos intermediarios grandes de um projeto depois que Salmon/STAR
terminou com sucesso. O que sera removido depende de PIPELINE_STORAGE_MODE e
das variaveis CLEANUP_* em config/user_settings.sh.

Opcoes:
  --scratch-root PATH    Default: SCRATCH_ROOT do config/pipeline_config.sh
  --storage-mode MODE    Override: full, balanced ou minimal
  --dry-run              Mostra o que seria removido, sem apagar nada
  -h, --help             Mostra esta ajuda
USAGE
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

PROJECT="$1"
shift

case "$PROJECT" in
    ""|"."|".."|*/*|*\\*)
        echo "[ERRO] Nome de projeto invalido para limpeza: '$PROJECT'" >&2
        exit 1
        ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PIPELINE_CONFIG:-}" && -f "$PIPELINE_CONFIG" ]]; then
    source "$PIPELINE_CONFIG"
elif [[ -n "${PROJECT_DIR:-}" && -f "${PROJECT_DIR}/config/pipeline_config.sh" ]]; then
    source "${PROJECT_DIR}/config/pipeline_config.sh"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config/pipeline_config.sh" ]]; then
    PROJECT_DIR="$(cd "$SLURM_SUBMIT_DIR" && pwd)"
    source "${PROJECT_DIR}/config/pipeline_config.sh"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]]; then
    PROJECT_DIR="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    source "${PROJECT_DIR}/config/pipeline_config.sh"
else
    PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
    source "${PROJECT_DIR}/config/pipeline_config.sh"
fi

DRY_RUN=0
SCRATCH_ROOT_VALUE="$SCRATCH_ROOT"
STORAGE_MODE_OVERRIDE=""

while [ $# -gt 0 ]; do
    case "$1" in
        --scratch-root)
            SCRATCH_ROOT_VALUE="$2"
            shift 2
            ;;
        --storage-mode)
            STORAGE_MODE_OVERRIDE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERRO] Opcao desconhecida: $1" >&2
            exit 1
            ;;
    esac
done

apply_storage_mode() {
    local mode="${1,,}"
    case "$mode" in
        full)
            export PIPELINE_STORAGE_MODE="full"
            export RUN_STORAGE_CLEANUP_AFTER_ALIGNMENT=0
            export CLEANUP_FASTQC_DIRS=0
            export CLEANUP_TRIMMED_RUNS=0
            export CLEANUP_TRIMMED_MERGED=0
            export CLEANUP_FASTQ_FTP=0
            export CLEANUP_STAR_BAM=0
            ;;
        balanced)
            export PIPELINE_STORAGE_MODE="balanced"
            export RUN_STORAGE_CLEANUP_AFTER_ALIGNMENT=1
            export CLEANUP_FASTQC_DIRS=1
            export CLEANUP_TRIMMED_RUNS=1
            export CLEANUP_TRIMMED_MERGED=0
            export CLEANUP_FASTQ_FTP=0
            export CLEANUP_STAR_BAM=0
            ;;
        minimal)
            export PIPELINE_STORAGE_MODE="minimal"
            export RUN_STORAGE_CLEANUP_AFTER_ALIGNMENT=1
            export CLEANUP_FASTQC_DIRS=1
            export CLEANUP_TRIMMED_RUNS=1
            export CLEANUP_TRIMMED_MERGED=1
            export CLEANUP_FASTQ_FTP=1
            export CLEANUP_STAR_BAM=1
            ;;
        *)
            echo "[ERRO] Storage mode invalido: $1. Use full, balanced ou minimal." >&2
            exit 1
            ;;
    esac
}

if [[ -n "$STORAGE_MODE_OVERRIDE" ]]; then
    apply_storage_mode "$STORAGE_MODE_OVERRIDE"
fi

truthy() {
    case "${1:-0}" in
        1|true|TRUE|yes|YES|y|Y)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

abs_child_path() {
    local path="$1"
    local parent base
    parent="$(dirname "$path")"
    base="$(basename "$path")"
    if [[ ! -d "$parent" ]]; then
        return 1
    fi
    printf '%s/%s\n' "$(cd "$parent" && pwd -P)" "$base"
}

safe_remove_dir() {
    local path="$1"
    local guard_root="$2"
    local label="$3"
    local abs_path

    if [[ ! -e "$path" ]]; then
        echo "[SKIP] ${label} nao existe: $path"
        return 0
    fi
    if [[ ! -d "$path" ]]; then
        echo "[WARN] ${label} nao e diretorio; ignorando: $path" >&2
        return 0
    fi

    abs_path="$(abs_child_path "$path")" || {
        echo "[ERRO] Nao foi possivel resolver caminho: $path" >&2
        exit 1
    }

    case "$abs_path" in
        "$guard_root"/*)
            ;;
        *)
            echo "[ERRO] Recusando remover fora da area esperada." >&2
            echo "[ERRO] Caminho: $abs_path" >&2
            echo "[ERRO] Area permitida: $guard_root" >&2
            exit 1
            ;;
    esac

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] rm -rf -- $abs_path"
    else
        echo "[CLEANUP] Removendo ${label}: $abs_path"
        rm -rf -- "$abs_path"
    fi
}

safe_remove_star_bams() {
    local root="$1"
    local guard_root="$2"
    local abs_root

    if [[ ! -d "$root" ]]; then
        echo "[SKIP] STAR BAM root nao existe: $root"
        return 0
    fi

    abs_root="$(cd "$root" && pwd -P)"
    case "$abs_root" in
        "$guard_root"|"$guard_root"/*)
            ;;
        *)
            echo "[ERRO] Recusando remover BAMs fora da area esperada." >&2
            echo "[ERRO] Caminho: $abs_root" >&2
            echo "[ERRO] Area permitida: $guard_root" >&2
            exit 1
            ;;
    esac

    if [ "$DRY_RUN" -eq 1 ]; then
        find "$abs_root" -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name '*.bai' \) -print |
            sed 's/^/[DRY-RUN] rm -- /'
    else
        echo "[CLEANUP] Removendo BAMs STAR em: $abs_root"
        find "$abs_root" -type f \( -name '*.bam' -o -name '*.bam.bai' -o -name '*.bai' \) -print -delete
    fi
}

PROJECT_SCRATCH="${SCRATCH_ROOT_VALUE}/${PROJECT}"

if [[ ! -d "$PROJECT_SCRATCH" ]]; then
    echo "[WARN] Scratch do projeto nao existe: $PROJECT_SCRATCH"
    echo "[WARN] Nada a limpar no scratch."
    PROJECT_SCRATCH_ABS=""
else
    PROJECT_SCRATCH_ABS="$(cd "$PROJECT_SCRATCH" && pwd -P)"
fi

echo "[INFO] Projeto: $PROJECT"
echo "[INFO] Storage mode: ${PIPELINE_STORAGE_MODE}"
echo "[INFO] Scratch root: $SCRATCH_ROOT_VALUE"
echo "[INFO] Dry-run: $DRY_RUN"

if [[ -n "$PROJECT_SCRATCH_ABS" ]]; then
    if truthy "$CLEANUP_FASTQC_DIRS"; then
        safe_remove_dir "${PROJECT_SCRATCH}/fastqc_raw" "$PROJECT_SCRATCH_ABS" "FastQC raw"
        safe_remove_dir "${PROJECT_SCRATCH}/fastqc_trimmed_runs" "$PROJECT_SCRATCH_ABS" "FastQC trimmed runs"
        safe_remove_dir "${PROJECT_SCRATCH}/fastqc_merged" "$PROJECT_SCRATCH_ABS" "FastQC merged"
    fi

    if truthy "$CLEANUP_TRIMMED_RUNS"; then
        safe_remove_dir "${PROJECT_SCRATCH}/trimmed_runs" "$PROJECT_SCRATCH_ABS" "trimmed runs"
    fi

    if truthy "$CLEANUP_TRIMMED_MERGED"; then
        safe_remove_dir "${PROJECT_SCRATCH}/trimmed_merged" "$PROJECT_SCRATCH_ABS" "trimmed merged"
    fi

    if truthy "$CLEANUP_FASTQ_FTP"; then
        safe_remove_dir "${PROJECT_SCRATCH}/fastq_ftp" "$PROJECT_SCRATCH_ABS" "FASTQ raw/download"
    fi
fi

if truthy "$CLEANUP_STAR_BAM"; then
    STAR_PROJECT_DIR="${STAR_QUANT_DIR}/${PROJECT}"
    if [[ -d "$STAR_PROJECT_DIR" ]]; then
        STAR_GUARD="$(cd "$STAR_PROJECT_DIR" && pwd -P)"
        safe_remove_star_bams "$STAR_PROJECT_DIR" "$STAR_GUARD"
    else
        echo "[SKIP] STAR quant do projeto nao existe: $STAR_PROJECT_DIR"
    fi
fi

echo "[OK] Limpeza concluida para $PROJECT"
