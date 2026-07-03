#!/bin/bash
#SBATCH --job-name=070_dtu
#SBATCH --output=logs/dtu/dtu_%j.out
#SBATCH --error=logs/dtu/dtu_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

ORIGINAL_ARGS=("$@")

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

SCRIPT_DIR="$DTU_SCRIPTS_DIR"
STEP_DIR="$DTU_DIR"
mkdir -p "$STEP_DIR/logs/dtu"
cd "$STEP_DIR"

DEFAULT_METADATA="$(metadata_default)"
METADATA="$DEFAULT_METADATA"
QUANT_ROOT="$QUANT_DIR"
TX2GENE="$TX2GENE_FILE"
OUTPUT_ROOT="$DTU_DIR"
PROJECTS="auto"
INCLUDE_ALL=0
ALLOW_MISSING=0
SBATCH_DRY_RUN=0
EXECUTOR="$PIPELINE_EXECUTOR"
TEST_VARIABLES="$DTU_TEST_VARIABLES"
MIN_REPLICATES="$DTU_MIN_REPLICATES"
MIN_GENE_COUNT="$DTU_MIN_GENE_COUNT"
MIN_TRANSCRIPTS="$DTU_MIN_TRANSCRIPTS_PER_GENE"
DEPENDENCY=""

usage() {
    echo "Uso: $0 [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH                  Default: METADATA_FINAL_NEW, fallback METADATA_FINAL"
    echo "  --quant-root PATH                Default: QUANT_DIR"
    echo "  --tx2gene PATH                   Default: TX2GENE_FILE"
    echo "  --output-root PATH               Default: DTU_DIR"
    echo "  --projects A,B|auto              Default: auto"
    echo "  --include-all                    Inclui analise combinada all_projects"
    echo "  --test-variables A,B             Default: DTU_TEST_VARIABLES"
    echo "  --min-replicates N               Default: DTU_MIN_REPLICATES"
    echo "  --min-gene-count N               Default: DTU_MIN_GENE_COUNT"
    echo "  --min-transcripts-per-gene N     Default: DTU_MIN_TRANSCRIPTS_PER_GENE"
    echo "  --dependency SPEC                Ex.: afterok:12345 quando submetendo manualmente"
    echo "  --executor slurm|local           Default: PIPELINE_EXECUTOR"
    echo "  --local                          Atalho para --executor local"
    echo "  --allow-missing                  Ignora amostras sem quant.sf"
    echo "  --sbatch-dry-run                 Mostra submissao/execucao sem executar"
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --metadata)
            METADATA=$2
            shift 2
            ;;
        --quant-root)
            QUANT_ROOT=$2
            shift 2
            ;;
        --tx2gene)
            TX2GENE=$2
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT=$2
            shift 2
            ;;
        --projects)
            PROJECTS=$2
            shift 2
            ;;
        --include-all)
            INCLUDE_ALL=1
            shift
            ;;
        --test-variables)
            TEST_VARIABLES=$2
            shift 2
            ;;
        --min-replicates)
            MIN_REPLICATES=$2
            shift 2
            ;;
        --min-gene-count)
            MIN_GENE_COUNT=$2
            shift 2
            ;;
        --min-transcripts-per-gene)
            MIN_TRANSCRIPTS=$2
            shift 2
            ;;
        --dependency)
            DEPENDENCY=$2
            shift 2
            ;;
        --executor)
            EXECUTOR=$2
            shift 2
            ;;
        --local)
            EXECUTOR="local"
            shift
            ;;
        --allow-missing)
            ALLOW_MISSING=1
            shift
            ;;
        --sbatch-dry-run|--dry-run)
            SBATCH_DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage 0
            ;;
        *)
            echo "[ERRO] Opcao desconhecida: $1"
            exit 1
            ;;
    esac
done

case "$EXECUTOR" in
    slurm|local)
        ;;
    *)
        echo "[ERRO] Executor invalido: $EXECUTOR. Use slurm ou local."
        exit 1
        ;;
esac

CMD=(
    Rscript "${SCRIPT_DIR}/dtu_analysis.R"
    --metadata "$METADATA"
    --quant-root "$QUANT_ROOT"
    --tx2gene "$TX2GENE"
    --output-root "$OUTPUT_ROOT"
    --projects "$PROJECTS"
    --test-variables "$TEST_VARIABLES"
    --min-replicates "$MIN_REPLICATES"
    --min-gene-count "$MIN_GENE_COUNT"
    --min-transcripts-per-gene "$MIN_TRANSCRIPTS"
)

if [ "$INCLUDE_ALL" -eq 1 ]; then
    CMD+=(--include-all)
fi

if [ "$ALLOW_MISSING" -eq 1 ]; then
    CMD+=(--allow-missing)
fi

echo "[INFO] Metadata: $METADATA"
echo "[INFO] Quant root: $QUANT_ROOT"
echo "[INFO] tx2gene: $TX2GENE"
echo "[INFO] Output root: $OUTPUT_ROOT"
echo "[INFO] Projetos: $PROJECTS"
echo "[INFO] Include all: $INCLUDE_ALL"
echo "[INFO] Executor: $EXECUTOR"

if [[ "$EXECUTOR" == "slurm" && -z "${SLURM_JOB_ID:-}" && "${DTU_RUNNER_SUBMITTED:-0}" != "1" ]]; then
    SUBMIT_CMD=(
        sbatch --parsable
        --chdir="$STEP_DIR"
        --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh,PIPELINE_EXECUTOR=slurm,DTU_RUNNER_SUBMITTED=1"
    )
    if [ -n "$DEPENDENCY" ]; then
        SUBMIT_CMD+=(--dependency="$DEPENDENCY")
    fi
    SUBMIT_CMD+=("${SCRIPT_DIR}/run_dtu_analysis_slurm.sh" "${ORIGINAL_ARGS[@]}")
    if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] Comando:"
        printf ' %q' "${SUBMIT_CMD[@]}"
        echo
        exit 0
    fi
    JOB_ID=$("${SUBMIT_CMD[@]}" | tail -n 1 | cut -d';' -f1)
    echo "[OK] Job DTU submetido: $JOB_ID"
    exit 0
fi

if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

activate_dtu_analysis
check_command Rscript
echo "+ ${CMD[*]}"
"${CMD[@]}"
echo "[OK] Etapa 070 DTU concluida."
