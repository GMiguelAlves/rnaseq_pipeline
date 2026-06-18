#!/bin/bash
#SBATCH --job-name=086_mfuzz
#SBATCH --output=logs/mfuzz/mfuzz_%j.out
#SBATCH --error=logs/mfuzz/mfuzz_%j.err
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

SCRIPT_DIR="$MFUZZ_SCRIPTS_DIR"
STEP_DIR="$MFUZZ_DIR"
mkdir -p "$STEP_DIR/logs/mfuzz"
cd "$STEP_DIR"

METADATA="$(metadata_default)"
EXPRESSION="$EXPRESSION_MATRIX_FILE"
SAMPLES="$QUANT_SAMPLES_FILE"
QUANTIFICATION_ROOT="$QUANTIFICATION_DIR"
OUTPUT_ROOT="$MFUZZ_DIR"
PROJECTS="auto"
INCLUDE_ALL=0
ALLOW_MISSING=0
SBATCH_DRY_RUN=0
EXECUTOR="$PIPELINE_EXECUTOR"
METHOD="$QUANT_METHOD"
TIME_VARIABLE="$MFUZZ_TIME_VARIABLE"
TIME_LEVELS="$MFUZZ_TIME_LEVELS"
GROUP_COLUMNS="$MFUZZ_GROUP_COLUMNS"
CLUSTERS="$MFUZZ_CLUSTERS"
M_VALUE="$MFUZZ_M"
MIN_SAMPLES="$MFUZZ_MIN_SAMPLES"
MIN_GENES="$MFUZZ_MIN_GENES"
MIN_EXPRESSION="$MFUZZ_MIN_EXPRESSION"
MIN_FRACTION="$MFUZZ_MIN_FRACTION"
DEPENDENCY=""

usage() {
    echo "Uso: $0 [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH              Default: METADATA_FINAL_NEW, fallback METADATA_FINAL"
    echo "  --expression PATH            Default: EXPRESSION_MATRIX_FILE"
    echo "  --samples PATH               Default: QUANT_SAMPLES_FILE"
    echo "  --quantification-dir PATH    Default: QUANTIFICATION_DIR"
    echo "  --output-root PATH           Default: MFUZZ_DIR"
    echo "  --projects A,B|auto          Default: auto"
    echo "  --include-all                Inclui all_projects"
    echo "  --method salmon|star         Default: QUANT_METHOD"
    echo "  --time-variable COL          Default: MFUZZ_TIME_VARIABLE"
    echo "  --time-levels A,B,C          Ordem temporal opcional"
    echo "  --group-columns A,B          Contextos adicionais antes do tempo"
    echo "  --clusters N                 Default: MFUZZ_CLUSTERS"
    echo "  --m N                        0 estima automaticamente"
    echo "  --min-samples N              Default: MFUZZ_MIN_SAMPLES"
    echo "  --min-genes N                Default: MFUZZ_MIN_GENES"
    echo "  --min-expression N           Default: MFUZZ_MIN_EXPRESSION"
    echo "  --min-fraction N             Default: MFUZZ_MIN_FRACTION"
    echo "  --dependency SPEC            Ex.: afterok:12345 quando submetendo manualmente"
    echo "  --executor slurm|local       Default: PIPELINE_EXECUTOR"
    echo "  --local                      Atalho para --executor local"
    echo "  --allow-missing              Ignora escopos sem matriz/tabela"
    echo "  --sbatch-dry-run             Mostra submissao/execucao sem executar"
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --metadata) METADATA=$2; shift 2 ;;
        --expression) EXPRESSION=$2; shift 2 ;;
        --samples) SAMPLES=$2; shift 2 ;;
        --quantification-dir) QUANTIFICATION_ROOT=$2; shift 2 ;;
        --output-root) OUTPUT_ROOT=$2; shift 2 ;;
        --projects) PROJECTS=$2; shift 2 ;;
        --include-all) INCLUDE_ALL=1; shift ;;
        --method) METHOD=$2; shift 2 ;;
        --time-variable) TIME_VARIABLE=$2; shift 2 ;;
        --time-levels) TIME_LEVELS=$2; shift 2 ;;
        --group-columns) GROUP_COLUMNS=$2; shift 2 ;;
        --clusters) CLUSTERS=$2; shift 2 ;;
        --m) M_VALUE=$2; shift 2 ;;
        --min-samples) MIN_SAMPLES=$2; shift 2 ;;
        --min-genes) MIN_GENES=$2; shift 2 ;;
        --min-expression) MIN_EXPRESSION=$2; shift 2 ;;
        --min-fraction) MIN_FRACTION=$2; shift 2 ;;
        --dependency) DEPENDENCY=$2; shift 2 ;;
        --executor) EXECUTOR=$2; shift 2 ;;
        --local) EXECUTOR="local"; shift ;;
        --allow-missing) ALLOW_MISSING=1; shift ;;
        --sbatch-dry-run|--dry-run) SBATCH_DRY_RUN=1; shift ;;
        -h|--help) usage 0 ;;
        *) echo "[ERRO] Opcao desconhecida: $1"; exit 1 ;;
    esac
done

case "$EXECUTOR" in
    slurm|local) ;;
    *) echo "[ERRO] Executor invalido: $EXECUTOR. Use slurm ou local."; exit 1 ;;
esac

CMD=(
    Rscript "${SCRIPT_DIR}/mfuzz_analysis.R"
    --metadata "$METADATA"
    --expression "$EXPRESSION"
    --samples "$SAMPLES"
    --quantification-dir "$QUANTIFICATION_ROOT"
    --output-root "$OUTPUT_ROOT"
    --projects "$PROJECTS"
    --method "$METHOD"
    --time-variable "$TIME_VARIABLE"
    --time-levels "$TIME_LEVELS"
    --group-columns "$GROUP_COLUMNS"
    --clusters "$CLUSTERS"
    --m "$M_VALUE"
    --min-samples "$MIN_SAMPLES"
    --min-genes "$MIN_GENES"
    --min-expression "$MIN_EXPRESSION"
    --min-fraction "$MIN_FRACTION"
)

[ "$INCLUDE_ALL" -eq 1 ] && CMD+=(--include-all)
[ "$ALLOW_MISSING" -eq 1 ] && CMD+=(--allow-missing)

echo "[INFO] Expression: $EXPRESSION"
echo "[INFO] Samples: $SAMPLES"
echo "[INFO] Output root: $OUTPUT_ROOT"
echo "[INFO] Projetos: $PROJECTS"
echo "[INFO] Include all: $INCLUDE_ALL"
echo "[INFO] Time variable: $TIME_VARIABLE"
echo "[INFO] Executor: $EXECUTOR"

if [[ "$EXECUTOR" == "slurm" && -z "${SLURM_JOB_ID:-}" && "${MFUZZ_RUNNER_SUBMITTED:-0}" != "1" ]]; then
    SUBMIT_CMD=(
        sbatch --parsable
        --chdir="$STEP_DIR"
        --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh,PIPELINE_EXECUTOR=slurm,MFUZZ_RUNNER_SUBMITTED=1"
    )
    [ -n "$DEPENDENCY" ] && SUBMIT_CMD+=(--dependency="$DEPENDENCY")
    SUBMIT_CMD+=("${SCRIPT_DIR}/run_mfuzz_analysis_slurm.sh" "${ORIGINAL_ARGS[@]}")
    if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] Comando:"
        printf ' %q' "${SUBMIT_CMD[@]}"
        echo
        exit 0
    fi
    JOB_ID=$("${SUBMIT_CMD[@]}" | tail -n 1 | cut -d';' -f1)
    echo "[OK] Job Mfuzz submetido: $JOB_ID"
    exit 0
fi

if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

activate_mfuzz_analysis
check_command Rscript
echo "+ ${CMD[*]}"
"${CMD[@]}"
echo "[OK] Etapa 086 Mfuzz concluida."
