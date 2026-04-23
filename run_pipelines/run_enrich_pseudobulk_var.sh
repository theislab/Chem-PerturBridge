#!/bin/bash
# Wrapper: validates -d; launcher logs + background & like run_deg.sh; Slurm job logs use %j in same folder tree.

set -e

LOGS_DIR=./logs
PIPELINE_NAME="enrich_pseudobulk_var"
SUBDIR="full"
DATASET=""
DEG_DIR=""
REFERENCE_H5AD=""
DRY=()
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
    "gdpx2"
)

mkdir -p "$LOGS_DIR"

EXTRA_ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "--dry-run" ]]; then
        DRY=(--dry-run)
    else
        EXTRA_ARGS+=("$arg")
    fi
done
set -- "${EXTRA_ARGS[@]}"

while getopts ":hd:i:r:" opt; do
    case $opt in
        h)
            echo "Run: $0 [-h] [--dry-run] -d DATASET -i DEG_DIR -r REFERENCE_H5AD"
            echo ""
            echo "  -d  Dataset: ${VALID_CHOICES[*]}"
            echo "  -i  Directory tree with DEG *_de.h5ad files (recursive scan)"
            echo "  -r  Processed pseudobulk .h5ad (reference .var)"
            echo ""
            echo "Launcher logs: logs/<dataset>/full/enrich_pseudobulk_var/enrich_<dataset>.PID*.out|.err"
            echo "Slurm job logs:  same folder, enrich_<dataset>.<jobid>.out|.err"
            echo "  --dry-run  Log only, do not write"
            exit 0
            ;;
        d)
            DATASET=$OPTARG
            ;;
        i)
            DEG_DIR=$OPTARG
            ;;
        r)
            REFERENCE_H5AD=$OPTARG
            ;;
        ?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$DATASET" ]]; then
    echo "Error: -d DATASET is required (-h for help)" >&2
    exit 1
fi

if [[ " ${VALID_CHOICES[*]} " != *" $DATASET "* ]]; then
    echo "Error: -d must be one of: ${VALID_CHOICES[*]}" >&2
    exit 1
fi

if [[ -z "$DEG_DIR" || -z "$REFERENCE_H5AD" ]]; then
    echo "Error: -i DEG_DIR and -r REFERENCE_H5AD are required (-h for help)" >&2
    exit 1
fi

LOG_BASE="${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}"
mkdir -p "$LOG_BASE"

./pipelines/common/enrich_pseudobulk_var.sh \
    -d "$DATASET" \
    -i "$DEG_DIR" \
    -r "$REFERENCE_H5AD" \
    "${DRY[@]}" \
    > "${LOG_BASE}/enrich_${DATASET}.PID$$.out" \
    2> "${LOG_BASE}/enrich_${DATASET}.PID$$.err" &