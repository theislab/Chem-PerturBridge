#!/bin/bash
set -e

if [ "$1" = "subsampling" ]; then
        echo "> Work with a subsample"
        SUFFIX="_subsample"
        SUBDIR="subsample"
        ARG="--subsampling"
	MEM="250G"
else
	echo "> Work with a full version"
        SUFFIX=""
        SUBDIR="full"
	ARG=""
	MEM="500G"
fi

DATASET=l1000_phase1

DATA_ROOT=./data/${DATASET}/raw/
OUTPUT_DIR=./data/${DATASET}/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Common sbatch preamble
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"
QOS=cpu_preemptible
PARTITION=cpu_p
#QOS=gpu_normal
#PARTITION=gpu_p

mkdir -p ${DATA_ROOT}
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}

echo "> Running ${DATASET} assembly"

sbatch -W -J pseudobulk_${DATASET} \
       -t 10:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=${MEM} \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.err \
       -o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.%j.out \
       --wrap="${SBATCH_PREAMBLE} && \
	       python3 -m src.pseudobulking.datasets.l1000.run_assembling \
	       --output_file ${OUTPUT_DIR}/${DATASET}_level3_deg_ready_landmark${SUFFIX}.h5ad \
	       --data_root ${DATA_ROOT} \
               --annotate-pubchem \
               --dataset ${DATASET} \
               ${ARG}"

echo "> ${DATASET} assembly is done"
