#!/bin/bash
set -e

if [ "$1" = "subsampling" ]; then
	echo "> Warning: subsampling mode is not supported for ${DATASET}, running full version"
fi

SUFFIX=""
SUBDIR="full"
MEM="500G"

DATASET=novartis

DATA_ROOT=./data/${DATASET}
RAW_DIR=${DATA_ROOT}/raw
OUTPUT_DIR=${DATA_ROOT}/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Common sbatch preamble
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"
QOS=cpu_normal
PARTITION=cpu_p

mkdir -p ${RAW_DIR}
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}

echo "> Running ${DATASET} standardization"

sbatch -W -J pseudobulk_${DATASET} \
       -t 10:00:00 \
       -n 2 \
       --qos=${QOS} \
       --mem=${MEM} \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
               python3 -m src.pseudobulking.datasets.novartis.run_standardization \
               --output_file ${OUTPUT_DIR}/${DATASET}_standardized${SUFFIX}.h5ad \
               --data_root ${DATA_ROOT} \
               --annotate-pubchem \
               --annotate-pubchem-names"

echo "> ${DATASET} standardization is done"
