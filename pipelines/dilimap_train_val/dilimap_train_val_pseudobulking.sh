#!/bin/bash
set -e

DATASET=dilimap_train_val

if [ "$1" = "subsampling" ]; then
	echo "> Warning: subsampling mode is not supported for ${DATASET}, running full version"
fi

DATA_ROOT=./data/${DATASET}/raw/
OUTPUT_DIR=./data/${DATASET}/pseudobulk/full
LOGS_DIR=./logs
ENV_DIR=./venv

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Common sbatch preamble
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"

QOS=cpu_normal
PARTITION=cpu_p

mkdir -p ${DATA_ROOT}
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}/${DATASET}/full

echo "> Running ${DATASET} assembly"
sbatch -W -J pseudobulk_${DATASET} \
       -t 02:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=150G \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/full/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/full/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
               python3 -m src.pseudobulking.datasets.dilimap.run_assembly \
               --mode train_val \
               --output_file ${OUTPUT_DIR}/${DATASET}_assembled.h5ad \
               --data_root ${DATA_ROOT} --annotate_pubchem"

echo "> ${DATASET} assembly is done"
