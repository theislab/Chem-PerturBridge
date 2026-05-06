#!/bin/bash
set -e

DATASET=vcpi_0002
EXPERIMENT_ID=vcpi-0002

if [ "$1" = "subsampling" ]; then
	echo "> Warning: subsampling mode is not implemented for ${DATASET}; running full pipeline"
fi

SUFFIX=""
SUBDIR="full"
MEM="250G"

DATA_ROOT=./data/${DATASET}
OUTPUT_DIR=${DATA_ROOT}/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"
QOS=${QOS:-normal}
PARTITION=${PARTITION:-compute}

mkdir -p ${DATA_ROOT}/raw
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}

echo "> Running ${DATASET} VCPI standardization (${EXPERIMENT_ID})"

sbatch -W -J pseudobulk_${DATASET} \
       -t 8:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=${MEM} \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
	       python3 -m src.pseudobulking.datasets.vcpi_ginkgo.run_standardization \
	       --experiment ${EXPERIMENT_ID} \
	       --output_file ${OUTPUT_DIR}/${DATASET}_standardized${SUFFIX}.h5ad \
	       --data_root ${DATA_ROOT} \
               --annotate-pubchem --annotate-pubchem-names"

echo "> ${DATASET} standardization job finished"
