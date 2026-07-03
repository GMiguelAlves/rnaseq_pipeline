#!/bin/bash
#SBATCH --job-name=rmats_plan
#SBATCH --output=logs/splicing/rmats_%A_%a.out
#SBATCH --error=logs/splicing/rmats_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Uso: sbatch --array=1-N $0 <SPLICING_PLAN.csv>"
    exit 1
fi

PLAN=$1

resolve_step_dir() {
    if [[ -n "${STEP_DIR:-}" && -f "${STEP_DIR}/generate_rmats_plan.py" ]]; then
        echo "$STEP_DIR"
        return 0
    fi
    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/generate_rmats_plan.py" ]]; then
        echo "$SLURM_SUBMIT_DIR"
        return 0
    fi
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${script_dir}/generate_rmats_plan.py" ]]; then
        echo "$script_dir"
        return 0
    fi
    return 1
}

resolve_pipeline_config() {
    if [[ -n "${PIPELINE_CONFIG:-}" && -f "$PIPELINE_CONFIG" ]]; then
        echo "$PIPELINE_CONFIG"
        return 0
    fi
    if [[ -n "${PROJECT_DIR:-}" && -f "${PROJECT_DIR}/config/pipeline_config.sh" ]]; then
        echo "${PROJECT_DIR}/config/pipeline_config.sh"
        return 0
    fi
    local step_dir
    step_dir="$(resolve_step_dir)" || return 1
    if [[ -f "${step_dir}/../../config/pipeline_config.sh" ]]; then
        echo "${step_dir}/../../config/pipeline_config.sh"
        return 0
    fi
    return 1
}

STEP_DIR_PATH="$(resolve_step_dir)" || {
    echo "[ERRO] Nao foi possivel localizar a etapa 080."
    exit 1
}

PIPELINE_CONFIG_PATH="$(resolve_pipeline_config)" || {
    echo "[ERRO] Nao foi possivel localizar config/pipeline_config.sh"
    exit 1
}

source "$PIPELINE_CONFIG_PATH"

mkdir -p "${SPLICING_DIR}/logs/splicing"
cd "$SPLICING_DIR"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

read -r ANALYSIS_ID SCOPE PROJECT VARIABLE LEVEL_A LEVEL_B B1_FILE B2_FILE GTF OUTPUT_DIR TMP_DIR < <(
    python -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1], newline='', encoding='utf-8'))); r=rows[int(sys.argv[2])-1]; print(r['analysis_id'], r['scope'], r['project'], r['variable'], r['level_a'], r['level_b'], r['b1_file'], r['b2_file'], r['gtf'], r['output_dir'], r['tmp_dir'], sep='\t')" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

echo "[INFO] Analise: $ANALYSIS_ID"
echo "[INFO] Escopo: $SCOPE/$PROJECT"
echo "[INFO] Variavel: $VARIABLE"
echo "[INFO] Comparacao: $LEVEL_A vs $LEVEL_B"
echo "[INFO] B1: $B1_FILE"
echo "[INFO] B2: $B2_FILE"
echo "[INFO] GTF: $GTF"
echo "[INFO] Output: $OUTPUT_DIR"

activate_splicing
check_command python

RMATS_CMD="${SPLICING_RMATS_COMMAND:-rmats.py}"
check_command "$RMATS_CMD"

THREADS="${SLURM_CPUS_PER_TASK:-${LOCAL_CPUS_PER_TASK:-8}}"
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

"$RMATS_CMD" \
    --b1 "$B1_FILE" \
    --b2 "$B2_FILE" \
    --gtf "$GTF" \
    --od "$OUTPUT_DIR" \
    --tmp "$TMP_DIR" \
    --readLength "$SPLICING_READ_LENGTH" \
    --nthread "$THREADS" \
    --libType "$SPLICING_LIB_TYPE"

SUMMARY_FILE="${OUTPUT_DIR}/splicing_summary.tsv${PIPELINE_TABLE_SUFFIX:-}"
SIGNIFICANT_FILE="${OUTPUT_DIR}/significant_events.tsv${PIPELINE_TABLE_SUFFIX:-}"

python "${STEP_DIR_PATH}/summarize_rmats_results.py" \
    --rmats-dir "$OUTPUT_DIR" \
    --output "$SUMMARY_FILE" \
    --significant-output "$SIGNIFICANT_FILE" \
    --fdr "$SPLICING_FDR_THRESHOLD"

echo "[OK] rMATS concluido: $ANALYSIS_ID"
