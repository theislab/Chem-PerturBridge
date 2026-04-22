#!/bin/bash
set -e

DATASET=cigs_mce
SOURCE=mce

if [ "$1" = "subsampling" ]; then
    echo "> Warning: subsampling mode is not supported for ${DATASET}, running full version"
fi

SUFFIX=""
SUBDIR="full"
MEM="500G"

DATA_ROOT=./data/${DATASET}
OUTPUT_DIR=${DATA_ROOT}/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"
QOS=cpu_normal
PARTITION=cpu_p

mkdir -p ${DATA_ROOT}/raw
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
               python3 -m src.pseudobulking.datasets.cigs.run_standardization \
               --source ${SOURCE} \
               --output_file ${OUTPUT_DIR}/${DATASET}${SUFFIX}.h5ad \
               --data_root ${DATA_ROOT} \
               --annotate-pubchem"

echo "> ${DATASET} standardization is done"
