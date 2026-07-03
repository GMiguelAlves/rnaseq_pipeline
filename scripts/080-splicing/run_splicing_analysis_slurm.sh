#!/bin/bash
#SBATCH --job-name=080_splicing_submit
#SBATCH --output=logs/splicing_submit_%j.out
#SBATCH --error=logs/splicing_submit_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

set -euo pipefail

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

SCRIPT_DIR="$SPLICING_SCRIPTS_DIR"
STEP_DIR="$SPLICING_DIR"
mkdir -p "$STEP_DIR" "$STEP_DIR/logs/splicing"
cd "$STEP_DIR"

DEFAULT_METADATA="$(metadata_default)"
METADATA="$DEFAULT_METADATA"
STAR_ROOT="$STAR_QUANT_DIR"
GTF="$REF_GTF"
OUTPUT_ROOT="$SPLICING_DIR"
PROJECTS="auto"
INCLUDE_ALL=0
ALLOW_MISSING=0
SBATCH_DRY_RUN=0
CONCURRENCY="$SPLICING_CONCURRENCY"
EXECUTOR="$PIPELINE_EXECUTOR"
PLAN=""
TEST_VARIABLES="$SPLICING_TEST_VARIABLES"
MIN_REPLICATES="$SPLICING_MIN_REPLICATES"
READ_LENGTH="$SPLICING_READ_LENGTH"
LIB_TYPE="$SPLICING_LIB_TYPE"
FDR_THRESHOLD="$SPLICING_FDR_THRESHOLD"
RMATS_COMMAND="$SPLICING_RMATS_COMMAND"
DEPENDENCY=""

usage() {
    echo "Uso: $0 [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH          Default: METADATA_FINAL_NEW, fallback METADATA_FINAL"
    echo "  --star-root PATH         Default: STAR_QUANT_DIR"
    echo "  --gtf PATH               Default: REF_GTF"
    echo "  --output-root PATH       Default: SPLICING_DIR"
    echo "  --projects A,B|auto      Default: auto"
    echo "  --include-all            Inclui analise combinada all_projects"
    echo "  --test-variables A,B     Default: SPLICING_TEST_VARIABLES"
    echo "  --min-replicates N       Default: SPLICING_MIN_REPLICATES"
    echo "  --read-length N          Default: SPLICING_READ_LENGTH"
    echo "  --lib-type TYPE          fr-unstranded, fr-firststrand, fr-secondstrand"
    echo "  --fdr N                  Default: SPLICING_FDR_THRESHOLD"
    echo "  --rmats-command CMD      Default: SPLICING_RMATS_COMMAND"
    echo "  --plan PATH              Default: work/rmats_plan.csv"
    echo "  --concurrency N          Default: SPLICING_CONCURRENCY"
    echo "  --dependency SPEC        Ex.: afterok:12345"
    echo "  --executor slurm|local   Default: PIPELINE_EXECUTOR"
    echo "  --local                  Atalho para --executor local"
    echo "  --allow-missing          Usa apenas BAMs existentes"
    echo "  --sbatch-dry-run         Mostra submissao/execucao sem executar"
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --metadata)
            METADATA=$2
            shift 2
            ;;
        --star-root)
            STAR_ROOT=$2
            shift 2
            ;;
        --gtf)
            GTF=$2
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
        --read-length)
            READ_LENGTH=$2
            shift 2
            ;;
        --lib-type)
            LIB_TYPE=$2
            shift 2
            ;;
        --fdr)
            FDR_THRESHOLD=$2
            shift 2
            ;;
        --rmats-command)
            RMATS_COMMAND=$2
            shift 2
            ;;
        --plan)
            PLAN=$2
            shift 2
            ;;
        --concurrency)
            CONCURRENCY=$2
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

if [ -z "$PLAN" ]; then
    PLAN="work/rmats_plan.csv"
fi

if [ -z "$GTF" ]; then
    echo "[ERRO] GTF ausente. Configure REF_GTF/GTF_URL ou passe --gtf."
    exit 1
fi

mkdir -p work logs/splicing

GEN_CMD=(
    python "${SCRIPT_DIR}/generate_rmats_plan.py"
    --metadata "$METADATA"
    --star-root "$STAR_ROOT"
    --gtf "$GTF"
    --output-root "$OUTPUT_ROOT"
    --projects "$PROJECTS"
    --test-variables "$TEST_VARIABLES"
    --min-replicates "$MIN_REPLICATES"
    --output "$PLAN"
)

if [ "$INCLUDE_ALL" -eq 1 ]; then
    GEN_CMD+=(--include-all)
fi

if [ "$ALLOW_MISSING" -eq 1 ]; then
    GEN_CMD+=(--allow-missing)
fi

echo "[INFO] Metadata: $METADATA"
echo "[INFO] STAR root: $STAR_ROOT"
echo "[INFO] GTF: $GTF"
echo "[INFO] Output root: $OUTPUT_ROOT"
echo "[INFO] Plano: $PLAN"
echo "[INFO] Projetos: $PROJECTS"
echo "[INFO] Include all: $INCLUDE_ALL"
echo "[INFO] Executor: $EXECUTOR"

activate_python_env
echo "+ ${GEN_CMD[*]}"
"${GEN_CMD[@]}"

COUNTS=$(python "${SCRIPT_DIR}/splicing_plan_counts.py" "$PLAN")
ANALYSES=$(echo "$COUNTS" | awk -F= '$1=="analyses" {print $2}')

if [ -z "$ANALYSES" ]; then
    echo "[ERRO] Nao foi possivel ler numero de analises em $PLAN"
    echo "$COUNTS"
    exit 1
fi

export SPLICING_READ_LENGTH="$READ_LENGTH"
export SPLICING_LIB_TYPE="$LIB_TYPE"
export SPLICING_FDR_THRESHOLD="$FDR_THRESHOLD"
export SPLICING_RMATS_COMMAND="$RMATS_COMMAND"

if [ "$EXECUTOR" = "local" ]; then
    if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] Exemplo de execucao local:"
        echo "  SLURM_ARRAY_TASK_ID=1 bash ${SCRIPT_DIR}/rmats_plan_job.sh $PLAN"
        exit 0
    fi
    echo "[INFO] Rodando rMATS localmente em ordem sequencial."
    run_local_array "rMATS" "$ANALYSES" "${SCRIPT_DIR}/rmats_plan_job.sh" "$PLAN"
    echo "[OK] Etapa 080 concluida localmente."
    exit 0
fi

CMD=(
    sbatch --parsable
    --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh,STEP_DIR=${SCRIPT_DIR},SPLICING_READ_LENGTH=${READ_LENGTH},SPLICING_LIB_TYPE=${LIB_TYPE},SPLICING_FDR_THRESHOLD=${FDR_THRESHOLD},SPLICING_RMATS_COMMAND=${RMATS_COMMAND}"
    --array="1-${ANALYSES}%${CONCURRENCY}"
)

if [ -n "$DEPENDENCY" ]; then
    CMD+=(--dependency="$DEPENDENCY")
fi

CMD+=("${SCRIPT_DIR}/rmats_plan_job.sh" "$PLAN")

if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

JOB_ID=$("${CMD[@]}" | tail -n 1 | cut -d';' -f1)
echo "[OK] Job rMATS submetido: $JOB_ID"

if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "[INFO] Aguardando jobs rMATS para liberar dependencias downstream."
    wait_for_slurm_jobs "$JOB_ID"
fi
