#!/bin/bash

set -e

# Wrapper for pipelines/common/deg_aggregation_supplementary.sh
# Resumes a DEG run at the aggregation step when deg.sh was terminated
# after the per-cell-type array jobs completed.
#
# Pass the SAME flags you used with run_deg.sh so paths line up.

LOGS_DIR=./logs
PIPELINE_NAME="deg"
MODE_S=False
MODE_J=False
MODE_Q=False
MODE_N=False
MODE_G=False
MODE_F=False
ARG_S=""
ARG_F=""
ARG_J=""
ARG_N=""
ARG_G=""
ARG_FORCE=""
DATASET=""
PAR=""
FILT=""

VALID_CHOICES=(
    "sciplex"
    "tahoe"
    "dilimap_train"
    "dilimap_train_val"
    "l1000_phase1"
    "l1000_phase2"
    "op3"
    "novartis"
    "vcpi_0001"
    "vcpi_0002"
    "vcpi_0003"
    "gdpx2"
    "cigs_mce"
    "cigs_tcm"
)

DEG_PARAMETERS=(
    "group_all_replicates"
    "separate_replicates"
)

mkdir -p ${LOGS_DIR}

while getopts ":sjd:p:f:qnhgF" opt; do
    case $opt in
        h)
            echo "Run: $0 [-s] [-j] [-h] [-g] [-F] [-f VALUE] [-q] [-n] -d (dataset) -p (parameter)"
            echo ""
            echo "Supplementary script: runs DEG Step 3 (intermediate -> results)"
            echo "after Step 2 array jobs already completed (e.g., parent deg.sh"
            echo "was terminated). Verifies completion via DEG .out logs first."
            echo ""
            echo "  -s Subsample of a dataset for debugging, default=false"
            echo "  -j Parallel mode (array jobs per cell type), default=false"
            echo "  -g Use GPU QOS/partition (gpu_normal/gpu_p), default=CPU"
            echo "  -F Force aggregation even if completion check fails or"
            echo "     results directory is non-empty"
            echo "  -h Help option"
            echo "  -f (value <int>) Min number of cells used in the original run"
            echo "  -q QC filter was used in the original run"
            echo "  -n Normalized data was used in the original run"
            echo "  -d (dataset <str>) Dataset name, required"
            echo "  -p (parameter <str>) DEG design parameter, required"
            echo "     One of: ${DEG_PARAMETERS[*]}"
            exit 0
            ;;
        s) MODE_S=True ;;
        j) MODE_J=True ;;
        p) PAR=$OPTARG ;;
        f) FILT=$OPTARG ;;
        q) MODE_Q=True ;;
        d) DATASET=$OPTARG ;;
        n) MODE_N=True ;;
        g) MODE_G=True ;;
        F) MODE_F=True ;;
    esac
done

if [[ " ${VALID_CHOICES[*]} " != *" $DATASET "* ]]; then
    echo "Error: -d must be set up and one of: ${VALID_CHOICES[*]}" >&2; exit 1
fi

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
    echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2; exit 1
fi

if ! [[ "$FILT" =~ ^[0-9]+$ ]] && ! [ -z "$FILT" ]; then
    echo "Error: -f filter value is not a number" >&2; exit 1
else
    if [[ "$FILT" =~ ^[0-9]+$ ]]; then
        ARG_F="-f $FILT"
    fi
fi

if [[ "$MODE_S" == "True" ]]; then
    SUBDIR="subsample"
    ARG_S="-s"
else
    SUBDIR="full"
fi

if [[ "$MODE_J" == "True" ]]; then ARG_J="-j"; fi
if [[ "$MODE_N" == "True" ]]; then ARG_N="-n"; fi
if [[ "$MODE_G" == "True" ]]; then ARG_G="-g"; fi
if [[ "$MODE_F" == "True" ]]; then ARG_FORCE="-F"; fi

if [[ -z "$FILT" ]]; then
    FILTER_FOLDER="filter_min_cells_0"
else
    FILTER_FOLDER="filter_min_cells_${FILT}"
fi

if [[ "$MODE_Q" == "True" ]]; then
    QC_FOLDER="qc_true"
    ARG_Q="-q"
else
    QC_FOLDER="qc_false"
    ARG_Q=""
fi

LOG_BASE=${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}
mkdir -p ${LOG_BASE}/aggregation

./pipelines/common/deg_aggregation_supplementary.sh $ARG_S -p $PAR $ARG_F $ARG_Q $ARG_J $ARG_N $ARG_G $ARG_FORCE \
    -c ./pipelines/${DATASET}/configs/deg/config.json \
    -d ${DATASET} \
        > ${LOG_BASE}/deg_aggregation_supplementary_${DATASET}.PID$$.out \
        2>${LOG_BASE}/deg_aggregation_supplementary_${DATASET}.PID$$.err &
