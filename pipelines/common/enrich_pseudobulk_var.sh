#!/bin/bash
# Enrich .var with symbol / is_merged from reference processed pseudobulk (left join).
# Recursively scans DEG_DIR for **/*_de.h5ad.
# Submits work via sbatch -W (like deg.sh).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${PROJECT_ROOT}"

ENV_DIR=./venv
QOS=cpu_normal
PARTITION=cpu_p

VALID_DATASETS=(
    "sciplex"
    "tahoe"
    "l1000_phase1"
    "l1000_phase2"
    "op3"
    "novartis"
    "vcpi_0001"
    "vcpi_0002"
    "gdpx2"
)

DATASET=""
DEG_DIR=""
REFERENCE_H5AD=""
DRY_RUN=()

usage() {
    echo "Run: $0 [-h] -d DATASET -i DEG_DIR -r REFERENCE_H5AD [--dry-run]"
    echo ""
    echo "  -d  Dataset (same as run_deg / run_enrich): ${VALID_DATASETS[*]}"
    echo "  -i  Directory tree containing DEG *_de.h5ad files (searched recursively)"
    echo "  -r  Processed pseudobulk .h5ad (must have .var symbol and is_merged)"
    echo ""
    echo "Slurm logs (same tree as deg): logs/<dataset>/full/enrich_pseudobulk_var/enrich_<dataset>.%j.{out,err}"
    echo "  --dry-run  Only log which files would be updated"
    echo "  -h  Help"
}

EXTRA_ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "--dry-run" ]]; then
        DRY_RUN=(--dry_run)
    else
        EXTRA_ARGS+=("$arg")
    fi
done
set -- "${EXTRA_ARGS[@]}"

while getopts ":hd:i:r:" opt; do
    case $opt in
        h)
            usage
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
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$DATASET" ]]; then
    echo "Error: -d DATASET is required" >&2
    usage >&2
    exit 1
fi

if [[ " ${VALID_DATASETS[*]} " != *" $DATASET "* ]]; then
    echo "Error: -d must be one of: ${VALID_DATASETS[*]}" >&2
    exit 1
fi

if [[ -z "$DEG_DIR" || -z "$REFERENCE_H5AD" ]]; then
    echo "Error: -i DEG_DIR and -r REFERENCE_H5AD are required" >&2
    usage >&2
    exit 1
fi

if [[ ! -d "$DEG_DIR" ]]; then
    echo "Error: DEG directory does not exist: $DEG_DIR" >&2
    exit 1
fi
if [[ ! -f "$REFERENCE_H5AD" ]]; then
    echo "Error: Reference h5ad not found: $REFERENCE_H5AD" >&2
    exit 1
fi

# Same as deg.sh: activate env on submit host before sbatch (compute node uses SBATCH_PREAMBLE too).
export PATH="${HOME}/miniforge3/bin:${PATH}"
eval "$(mamba shell hook --shell bash)"
mamba activate "${ENV_DIR}"

echo "> DEG directory (scan *_de.h5ad): $DEG_DIR"
echo "> Reference pseudobulk: $REFERENCE_H5AD"

DRY_PART=""
if [[ ${#DRY_RUN[@]} -gt 0 ]]; then
    DRY_PART=" --dry_run"
fi

DEG_Q=$(printf '%q' "$DEG_DIR")
REF_Q=$(printf '%q' "$REFERENCE_H5AD")

# Same layout as run_deg: logs/${DATASET}/full/${PIPELINE_NAME}/...
SLURM_LOG_DIR="${PROJECT_ROOT}/logs/${DATASET}/full/enrich_pseudobulk_var"
SLURM_OUT="${SLURM_LOG_DIR}/enrich_${DATASET}.%j.out"
SLURM_ERR="${SLURM_LOG_DIR}/enrich_${DATASET}.%j.err"
mkdir -p "${SLURM_LOG_DIR}"

# Same pattern as pipelines/common/deg.sh (preprocess / aggregation jobs)
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && export TMPDIR=\${HOME}/tmp && mkdir -p \${TMPDIR} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"

echo "> Submitting Slurm job (sbatch -W) enrich_pseudobulk_var"
echo "> Slurm logs: ${SLURM_OUT} / ${SLURM_ERR}"

sbatch -W -J enrich_pseudobulk_var \
    --partition="${PARTITION}" \
    --qos="${QOS}" \
    --mem=250G \
    --time=12:00:00 \
    --cpus-per-task=4 \
    --exclude=supercpu01,gpusrv36,gpusrv56,gpusrv60,gpusrv45 \
    --output="${SLURM_OUT}" \
    --error="${SLURM_ERR}" \
    --wrap="${SBATCH_PREAMBLE} && python3 -m src.deg.enrich_pseudobulk_var_metadata --deg_dir ${DEG_Q} --reference_h5ad ${REF_Q}${DRY_PART}"

echo "> enrich_pseudobulk_var Slurm job finished"
