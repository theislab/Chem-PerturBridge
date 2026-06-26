#!/bin/bash

set -e

# ============================================================================
# DEG Aggregation Supplementary Pipeline
# ============================================================================
#
# Resumes the DEG pipeline AFTER the per-cell-type array jobs have completed
# but the parent deg.sh process was terminated before STEP 3 (aggregation).
#
# This script:
#   1. Reconstructs the same paths that deg.sh would have used (from flags).
#   2. Searches the DEG .out logs for the "DE analysis completed!" marker
#      and verifies that the number of completed array tasks matches the
#      number of expected pseudobulk input files.
#   3. Submits the aggregation sbatch (intermediate/ -> results/).
#
# It accepts the same path-affecting flags as deg.sh so it can be pointed at
# any (dataset, parameter, qc, filter, subsample, parallel, normalized) run.
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PIPELINE_NAME="deg"
ENV_DIR=./venv
LOGS_DIR=./logs
QOS=cpu_normal
PARTITION=cpu_p

MODE_S=False
MODE_J=False
MODE_Q=False
MODE_N=False
MODE_G=False
MODE_F=False
PAR=""
FILT=""
CONFIG=""
DATASET=""

DEG_PARAMETERS=(
    "group_all_replicates"
    "separate_replicates"
)

while getopts ":sjqp:f:c:d:nhgF" opt; do
    case $opt in
        h)
            echo "Run: $0 [-s] [-j] [-h] [-g] [-F] [-f] [-q] [-n] -p (parameters: ${DEG_PARAMETERS[*]}) -c CONFIG -d DATASET"
            echo ""
            echo "Resume DEG pipeline at the aggregation step (after per-cell-type"
            echo "array jobs completed but deg.sh was terminated before Step 3)."
            echo ""
            echo "Required arguments:"
            echo "  -d  Dataset name (e.g., vcpi_0003)"
            echo "  -p  Design parameter (group_all_replicates, separate_replicates)"
            echo "  -c  Path to config file (auto-set by run_deg_aggregation_supplementary.sh)"
            echo ""
            echo "Optional flags (must match the original deg.sh invocation):"
            echo "  -s  Subsample mode"
            echo "  -j  Parallel-by-cell-type mode (array jobs were run per cell type)"
            echo "  -g  Use GPU QOS/partition (gpu_normal/gpu_p)"
            echo "  -f  VALUE - Filter: minimum cells per sample"
            echo "  -q  Quality control filter"
            echo "  -n  Normalized dataset"
            echo "  -F  Force aggregation even if completion check fails"
            echo "  -h  Show this help message"
            exit 0
            ;;
        s) MODE_S=True ;;
        j) MODE_J=True ;;
        p) PAR=$OPTARG ;;
        f) FILT=$OPTARG ;;
        q) MODE_Q=True ;;
        c) CONFIG=$OPTARG ;;
        d) DATASET=$OPTARG ;;
        n) MODE_N=True ;;
        g) MODE_G=True ;;
        F) MODE_F=True ;;
    esac
done

if [ "$MODE_G" = True ]; then
    QOS=gpu_normal
    PARTITION=gpu_p
fi

# ============================================================================
# Validate arguments
# ============================================================================

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
    echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2
    exit 1
fi

if [ -z "$DATASET" ]; then
    echo "Error: Dataset name must be provided via -d flag" >&2
    exit 1
fi

if ! [[ "$FILT" =~ ^[0-9]+$ ]] && ! [ -z "$FILT" ]; then
    echo "Error: -f filter value is not a number" >&2
    exit 1
fi

if [[ -z "$FILT" ]]; then
    FILTER_FOLDER="filter_min_cells_0"
else
    FILTER_FOLDER="filter_min_cells_${FILT}"
fi

if [ "$MODE_Q" = "True" ]; then
    QC_FOLDER="qc_true"
else
    QC_FOLDER="qc_false"
fi

if [ "$MODE_S" = "True" ]; then
    SUBDIR="subsample"
else
    SUBDIR="full"
fi

if [ "$MODE_J" = "True" ]; then
    echo "> Mode: PARALLEL BY CELL TYPE"
else
    echo "> Mode: SEQUENTIAL BY CELL TYPE"
fi

if [ "$MODE_N" = "True" ]; then
    echo "> Dataset is normalized"
fi

# ============================================================================
# Initialize environment
# ============================================================================

export PATH="$HOME/miniforge3/bin:$PATH"

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

if [ -z "$CONFIG" ]; then
    echo "Error: Config path must be provided via -c flag" >&2
    exit 1
fi

if [ ! -f "$CONFIG" ]; then
    echo "Error: Config file not found: $CONFIG" >&2
    exit 1
fi

echo "> Using config: $CONFIG"

# ============================================================================
# Resolve paths (must match deg.sh)
# ============================================================================

LOGS_BASE_DIR="${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}"
DEG_LOG_DIR="${LOGS_BASE_DIR}/deg"

if ! par_deg=$(jq -e ".$PAR.par_deg" "$CONFIG"); then
    echo "Error: Failed to extract par_deg from $CONFIG" >&2
    exit 1
fi

INPUT_DIR=$(echo "$par_deg" | jq -r '.input_dir')
OUTPUT_DIR=$(echo "$par_deg" | jq -r '.output_dir')

if [ "$MODE_J" = "True" ]; then
    SEARCH_DIR="${INPUT_DIR}/by_celltype"
else
    SEARCH_DIR="${INPUT_DIR}"
fi

BASE_OUTPUT_DIR="${OUTPUT_DIR}/${SUBDIR}/${QC_FOLDER}/${FILTER_FOLDER}"
RESULTS_DIR="${BASE_OUTPUT_DIR}/results"
INTERMEDIATE_DIR="${BASE_OUTPUT_DIR}/intermediate"

echo "> Paths:"
echo "    SEARCH_DIR       = ${SEARCH_DIR}"
echo "    DEG_LOG_DIR      = ${DEG_LOG_DIR}"
echo "    INTERMEDIATE_DIR = ${INTERMEDIATE_DIR}"
echo "    RESULTS_DIR      = ${RESULTS_DIR}"

# ============================================================================
# STEP A: Verify DEG completion from logs
# ============================================================================

echo ""
echo "============================================"
echo "VERIFYING DEG COMPLETION FROM LOGS"
echo "============================================"

if [ ! -d "$SEARCH_DIR" ]; then
    echo "Error: Pseudobulk input directory not found: $SEARCH_DIR" >&2
    exit 1
fi

if [ "$MODE_J" = "True" ]; then
    mapfile -t INPUT_FILES < <(find "$SEARCH_DIR" -mindepth 2 -maxdepth 2 -name "*.h5ad" -type f 2>/dev/null | sort)
else
    mapfile -t INPUT_FILES < <(find "$SEARCH_DIR" -maxdepth 1 -name "*.h5ad" -type f 2>/dev/null | sort)
fi
N_INPUTS=${#INPUT_FILES[@]}

if [ "$N_INPUTS" -eq 0 ]; then
    echo "Error: No .h5ad input files found in ${SEARCH_DIR}" >&2
    exit 1
fi

echo "  Expected DEG tasks (from input files):  ${N_INPUTS}"

if [ ! -d "$DEG_LOG_DIR" ]; then
    echo "Error: DEG log directory not found: ${DEG_LOG_DIR}" >&2
    exit 1
fi

# Find all DEG .out logs (per array task; possibly in cell-type subdirs)
mapfile -t DEG_OUT_FILES < <(find "$DEG_LOG_DIR" -type f -name "deg_analysis.*.out" 2>/dev/null | sort)
N_OUT=${#DEG_OUT_FILES[@]}

if [ "$N_OUT" -eq 0 ]; then
    echo "Error: No DEG .out logs found under ${DEG_LOG_DIR}" >&2
    exit 1
fi

# Success marker emitted by run_deg.R when the run finishes
COMPLETION_MARKER="DE analysis completed!"

# Track which (cell_type, array_task_id) pairs completed successfully.
# We deduplicate so that re-runs of the same task don't inflate the count.
declare -A COMPLETED_TASKS
N_COMPLETED_LOGS=0
N_INCOMPLETE_LOGS=0
INCOMPLETE_LOGS=()

for OUT_FILE in "${DEG_OUT_FILES[@]}"; do
    if grep -Fq "$COMPLETION_MARKER" "$OUT_FILE"; then
        N_COMPLETED_LOGS=$((N_COMPLETED_LOGS + 1))

        # Build a (cell_type, task_id) key.
        # Logs live in either:
        #   ${DEG_LOG_DIR}/deg_analysis.<jobid>_<taskid>.out               (no cell type)
        #   ${DEG_LOG_DIR}/<CELL_TYPE>/deg_analysis.<jobid>_<taskid>.out
        REL_PATH="${OUT_FILE#${DEG_LOG_DIR}/}"
        REL_DIR=$(dirname "$REL_PATH")
        if [ "$REL_DIR" = "." ]; then
            CELL_TYPE_KEY="__nocelltype__"
        else
            CELL_TYPE_KEY="$REL_DIR"
        fi
        BASE=$(basename "$OUT_FILE" .out)         # deg_analysis.<jobid>_<taskid>
        TASK_ID="${BASE##*_}"
        COMPLETED_TASKS["${CELL_TYPE_KEY}/${TASK_ID}"]=1
    else
        N_INCOMPLETE_LOGS=$((N_INCOMPLETE_LOGS + 1))
        INCOMPLETE_LOGS+=("$OUT_FILE")
    fi
done

N_COMPLETED_TASKS=${#COMPLETED_TASKS[@]}

echo "  DEG .out logs found:                    ${N_OUT}"
echo "  Logs with '${COMPLETION_MARKER}':       ${N_COMPLETED_LOGS}"
echo "  Unique completed (cell_type, task_id):  ${N_COMPLETED_TASKS}"
echo "  Logs WITHOUT completion marker:         ${N_INCOMPLETE_LOGS}"

if [ "$N_INCOMPLETE_LOGS" -gt 0 ]; then
    echo "  Incomplete log files:"
    for f in "${INCOMPLETE_LOGS[@]}"; do
        echo "    - $f"
    done
fi

if [ "$N_COMPLETED_TASKS" -lt "$N_INPUTS" ]; then
    echo ""
    echo "WARNING: Only ${N_COMPLETED_TASKS}/${N_INPUTS} unique DEG tasks completed."
    if [ "$MODE_F" != "True" ]; then
        echo "Refusing to aggregate. Re-run the missing tasks, or pass -F to force." >&2
        exit 1
    else
        echo "  -F (force) was set; proceeding with aggregation anyway."
    fi
else
    echo ""
    echo "  All ${N_INPUTS} DEG tasks completed successfully."
fi

# ============================================================================
# STEP B: Aggregate intermediate -> results
# ============================================================================

echo ""
echo "============================================"
echo "AGGREGATING BATCH FILES"
echo "============================================"

SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && export TMPDIR=\${HOME}/tmp && mkdir -p \${TMPDIR} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"

if [ ! -d "$INTERMEDIATE_DIR" ] || [ -z "$(ls -A "$INTERMEDIATE_DIR" 2>/dev/null)" ]; then
    echo "Error: Intermediate directory missing or empty: ${INTERMEDIATE_DIR}" >&2
    exit 1
fi

if [ -d "$RESULTS_DIR" ] && [ -n "$(ls -A "$RESULTS_DIR" 2>/dev/null)" ]; then
    if [ "$MODE_F" != "True" ]; then
        echo "Results directory is not empty: ${RESULTS_DIR}"
        echo "Refusing to overwrite. Pass -F to force re-aggregation." >&2
        exit 1
    else
        echo "  -F (force) was set; re-aggregating despite non-empty results dir."
    fi
fi

mkdir -p "${LOGS_BASE_DIR}/aggregation"
mkdir -p "${RESULTS_DIR}"

echo "  Submitting aggregation job..."
sbatch -W -J deg_aggregation \
    --partition=${PARTITION} \
    --qos=${QOS} \
    --mem=500G \
    --time=2:00:00 \
    --cpus-per-task=2 \
    --exclude=supercpu01,gpusrv36,gpusrv56,gpusrv60,gpusrv45 \
    --output="${LOGS_BASE_DIR}/aggregation/deg_aggregation.%j.out" \
    --error="${LOGS_BASE_DIR}/aggregation/deg_aggregation.%j.err" \
    --wrap="${SBATCH_PREAMBLE} && \
            python3 -m src.deg.aggregating_deg \
            --input_dir ${INTERMEDIATE_DIR} \
            --output_dir ${RESULTS_DIR}"

echo "> Aggregation completed"

echo ""
echo "============================================"
echo "DEG AGGREGATION SUPPLEMENTARY COMPLETED"
echo "============================================"
echo "  Results saved to: ${RESULTS_DIR}"
