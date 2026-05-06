#!/bin/bash
set -e

DATASET=op3

if [ "$1" = "subsampling" ]; then
	echo "> Warning: subsampling mode is not supported for ${DATASET}, running full version"
fi

SUFFIX=""
SUBDIR="full"
MEM="100G"

DATA_ROOT=./data/${DATASET}/raw/
OUTPUT_DIR=./data/${DATASET}/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Common sbatch preamble
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"
QOS=${QOS:-normal}
PARTITION=${PARTITION:-compute}


mkdir -p ${DATA_ROOT}
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}

echo "> Running ${DATASET} standardization"

sbatch -W -J pseudobulk_${DATASET} \
       -t 5:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=${MEM} \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
	       python3 -m src.pseudobulking.datasets.op3.run_standardization \
	       --output_file ${OUTPUT_DIR}/${DATASET}_standardized${SUFFIX}.h5ad \
	       --data_root ${DATA_ROOT} \
               --annotate-pubchem"

echo "> ${DATASET} standardization is done"
