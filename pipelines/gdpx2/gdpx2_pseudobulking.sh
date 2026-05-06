#!/bin/bash
set -e

DATASET=gdpx2

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

echo "> Running ${DATASET} standardization"

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
           echo '> Downloading GDPx2' && \
           python3 -m src.downloading.run_downloading_datasets \
               --input \
                   'key_adata=ginkgo-datapoints/vcpi/X.h5ad' \
                   'key_obs=ginkgo-datapoints/vcpi/obs.parquet' \
                   'key_var=ginkgo-datapoints/vcpi/var.parquet' \
               --output \
                   'path2adata=${DATA_ROOT}/raw/gdpx2_raw.h5ad' \
                   'path2obs=${DATA_ROOT}/raw/gdpx2_obs.parquet' \
                   'path2var=${DATA_ROOT}/raw/gdpx2_var.parquet' && \
           echo '> Standardizing GDPx2 dataset' && \
           python3 -m src.pseudobulking.datasets.gdpx2.run_standardization \
               --data_root ${DATA_ROOT} \
               --output_file ${OUTPUT_DIR}/${DATASET}_standardized${SUFFIX}.h5ad \
               --annotate-pubchem"

echo "> ${DATASET} standardization job finished"
