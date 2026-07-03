#!/bin/bash
#SBATCH --job-name=085_wgcna
#SBATCH --output=logs/wgcna/wgcna_%j.out
#SBATCH --error=logs/wgcna/wgcna_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
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

SCRIPT_DIR="$WGCNA_SCRIPTS_DIR"
STEP_DIR="$WGCNA_DIR"
mkdir -p "$STEP_DIR/logs/wgcna"
cd "$STEP_DIR"

METADATA="$(metadata_default)"
EXPRESSION="$EXPRESSION_MATRIX_FILE"
SAMPLES="$QUANT_SAMPLES_FILE"
QUANTIFICATION_ROOT="$QUANTIFICATION_DIR"
OUTPUT_ROOT="$WGCNA_DIR"
PROJECTS="auto"
INCLUDE_ALL=0
ALLOW_MISSING=0
SBATCH_DRY_RUN=0
EXECUTOR="$PIPELINE_EXECUTOR"
METHOD="$QUANT_METHOD"
TRAIT_COLUMNS="$WGCNA_TRAIT_COLUMNS"
MIN_SAMPLES="$WGCNA_MIN_SAMPLES"
MIN_GENES="$WGCNA_MIN_GENES"
MIN_EXPRESSION="$WGCNA_MIN_EXPRESSION"
MIN_FRACTION="$WGCNA_MIN_FRACTION"
POWER="$WGCNA_POWER"
FIT_CUTOFF="$WGCNA_FIT_CUTOFF"
NETWORK_TYPE="$WGCNA_NETWORK_TYPE"
COR_METHOD="$WGCNA_COR_METHOD"
MIN_MODULE_SIZE="$WGCNA_MIN_MODULE_SIZE"
MERGE_CUT_HEIGHT="$WGCNA_MERGE_CUT_HEIGHT"
TOP_HUBS="$WGCNA_TOP_HUBS"
MAX_BLOCK_SIZE="$WGCNA_MAX_BLOCK_SIZE"
DEPENDENCY=""

usage() {
    echo "Uso: $0 [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH              Default: METADATA_FINAL_NEW, fallback METADATA_FINAL"
    echo "  --expression PATH            Default: EXPRESSION_MATRIX_FILE"
    echo "  --samples PATH               Default: QUANT_SAMPLES_FILE"
    echo "  --quantification-dir PATH    Default: QUANTIFICATION_DIR"
    echo "  --output-root PATH           Default: WGCNA_DIR"
    echo "  --projects A,B|auto          Default: auto"
    echo "  --include-all                Inclui all_projects"
    echo "  --method salmon|star         Default: QUANT_METHOD"
    echo "  --trait-columns A,B          Traits para correlacao modulo-traco"
    echo "  --min-samples N              Default: WGCNA_MIN_SAMPLES"
    echo "  --min-genes N                Default: WGCNA_MIN_GENES"
    echo "  --min-expression N           Default: WGCNA_MIN_EXPRESSION"
    echo "  --min-fraction N             Default: WGCNA_MIN_FRACTION"
    echo "  --power N                    0 seleciona automaticamente"
    echo "  --fit-cutoff N               Default: WGCNA_FIT_CUTOFF"
    echo "  --network-type TYPE          signed, unsigned, signed-hybrid"
    echo "  --cor-method TYPE            pearson ou bicor"
    echo "  --min-module-size N          Default: WGCNA_MIN_MODULE_SIZE"
    echo "  --merge-cut-height N         Default: WGCNA_MERGE_CUT_HEIGHT"
    echo "  --top-hubs N                 Default: WGCNA_TOP_HUBS"
    echo "  --max-block-size N           Default: WGCNA_MAX_BLOCK_SIZE"
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
        --trait-columns) TRAIT_COLUMNS=$2; shift 2 ;;
        --min-samples) MIN_SAMPLES=$2; shift 2 ;;
        --min-genes) MIN_GENES=$2; shift 2 ;;
        --min-expression) MIN_EXPRESSION=$2; shift 2 ;;
        --min-fraction) MIN_FRACTION=$2; shift 2 ;;
        --power) POWER=$2; shift 2 ;;
        --fit-cutoff) FIT_CUTOFF=$2; shift 2 ;;
        --network-type) NETWORK_TYPE=$2; shift 2 ;;
        --cor-method) COR_METHOD=$2; shift 2 ;;
        --min-module-size) MIN_MODULE_SIZE=$2; shift 2 ;;
        --merge-cut-height) MERGE_CUT_HEIGHT=$2; shift 2 ;;
        --top-hubs) TOP_HUBS=$2; shift 2 ;;
        --max-block-size) MAX_BLOCK_SIZE=$2; shift 2 ;;
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
    Rscript "${SCRIPT_DIR}/wgcna_analysis.R"
    --metadata "$METADATA"
    --expression "$EXPRESSION"
    --samples "$SAMPLES"
    --quantification-dir "$QUANTIFICATION_ROOT"
    --output-root "$OUTPUT_ROOT"
    --projects "$PROJECTS"
    --method "$METHOD"
    --trait-columns "$TRAIT_COLUMNS"
    --min-samples "$MIN_SAMPLES"
    --min-genes "$MIN_GENES"
    --min-expression "$MIN_EXPRESSION"
    --min-fraction "$MIN_FRACTION"
    --power "$POWER"
    --fit-cutoff "$FIT_CUTOFF"
    --network-type "$NETWORK_TYPE"
    --cor-method "$COR_METHOD"
    --min-module-size "$MIN_MODULE_SIZE"
    --merge-cut-height "$MERGE_CUT_HEIGHT"
    --top-hubs "$TOP_HUBS"
    --max-block-size "$MAX_BLOCK_SIZE"
)

[ "$INCLUDE_ALL" -eq 1 ] && CMD+=(--include-all)
[ "$ALLOW_MISSING" -eq 1 ] && CMD+=(--allow-missing)

echo "[INFO] Expression: $EXPRESSION"
echo "[INFO] Samples: $SAMPLES"
echo "[INFO] Output root: $OUTPUT_ROOT"
echo "[INFO] Projetos: $PROJECTS"
echo "[INFO] Include all: $INCLUDE_ALL"
echo "[INFO] Executor: $EXECUTOR"

if [[ "$EXECUTOR" == "slurm" && -z "${SLURM_JOB_ID:-}" && "${WGCNA_RUNNER_SUBMITTED:-0}" != "1" ]]; then
    SUBMIT_CMD=(
        sbatch --parsable
        --chdir="$STEP_DIR"
        --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh,PIPELINE_EXECUTOR=slurm,WGCNA_RUNNER_SUBMITTED=1"
    )
    [ -n "$DEPENDENCY" ] && SUBMIT_CMD+=(--dependency="$DEPENDENCY")
    SUBMIT_CMD+=("${SCRIPT_DIR}/run_wgcna_analysis_slurm.sh" "${ORIGINAL_ARGS[@]}")
    if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
        echo "[DRY-RUN] Comando:"
        printf ' %q' "${SUBMIT_CMD[@]}"
        echo
        exit 0
    fi
    JOB_ID=$("${SUBMIT_CMD[@]}" | tail -n 1 | cut -d';' -f1)
    echo "[OK] Job WGCNA submetido: $JOB_ID"
    exit 0
fi

if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

activate_wgcna_analysis
check_command Rscript
echo "+ ${CMD[*]}"
"${CMD[@]}"
echo "[OK] Etapa 085 WGCNA concluida."
