#!/bin/bash
set -e

if [ "$1" = "subsampling" ]; then
	echo "> Work with a subsample"
	SUFFIX="_subsample"
	SUBDIR="subsample"
	ARG="--subsampling"
else
	echo "> Work with a full version"
	SUFFIX=""
	SUBDIR="full"
	ARG=""
fi

DATASET=dilimap_train
DATA_ROOT=./data/${DATASET}/raw/
OUTPUT_DIR=./data/${DATASET}/pseudobulk/${SUBDIR}
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
mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}

echo "> Running ${DATASET} assembly"
sbatch -W -J pseudobulk_${DATASET} \
       -t 02:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=100G \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
               python3 -m src.pseudobulking.datasets.dilimap_train.run_assembly \
               --mode train \
               --output_file ${OUTPUT_DIR}/${DATASET}_assembled${SUFFIX}.h5ad \
               --data_root ${DATA_ROOT} ${ARG}"

echo "> ${DATASET} assembly is done"
