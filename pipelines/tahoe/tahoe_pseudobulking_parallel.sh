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
	MEM="700G"
fi

SC_DIR=./data/tahoe/raw
BULK_DIR=./data/tahoe/pseudobulk_to_merge/${SUBDIR}
OUTPUT_DIR=./data/tahoe/pseudobulk/${SUBDIR}
LOGS_DIR=./logs
ENV_DIR=./venv
QOS=${QOS:-normal}
PARTITION=${PARTITION:-compute}

echo "> Download dataset"
sbatch -W -J pseudobulk_tahoe \
	-t 3:00:00 \
	-n 1 \
	--array=1-14 \
	--qos=${QOS} \
	--mem=250G \
	--partition=${PARTITION} \
	--cpus-per-task=2 \
	-e ${LOGS_DIR}/tahoe/${SUBDIR}/plate%a/pseudobulk_tahoe_download.%A_%a.err \
	-o ${LOGS_DIR}/tahoe/${SUBDIR}/plate%a/pseudobulk_tahoe_download.%A_%a.out \
	--wrap="eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
		python3 -m src.downloading.run_downloading_datasets \
		--input key_adata=2025-02-25/h5ad/plate\${SLURM_ARRAY_TASK_ID}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad \
			key_obs=tahoe100/obs/plate\${SLURM_ARRAY_TASK_ID}.parquet \
			key_var=tahoe100/var.parquet \
		--output path2adata=${SC_DIR}/plate\${SLURM_ARRAY_TASK_ID}.h5ad \
                	path2obs=${SC_DIR}/plate\${SLURM_ARRAY_TASK_ID}_obs.parquet \
			path2var=${SC_DIR}/var.parquet ${ARG}"


echo "> Running pseudobulking"

sbatch -W -J pseudobulk_tahoe \
       -t 10:00:00 \
       -n 1 \
       --array=1-14 \
       --qos=${QOS} \
       --mem=${MEM} \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/tahoe/${SUBDIR}/plate%a/pseudobulk_tahoe_process.%A_%a.err \
       -o ${LOGS_DIR}/tahoe/${SUBDIR}/plate%a/pseudobulk_tahoe_process.%A_%a.out \
       --wrap="eval \"\$(mamba shell hook --shell bash)\" && \
               mamba activate ${ENV_DIR} && \
	       python3 -m src.pseudobulking.common.run_pseudobulking \
	       --dataset_name tahoe \
	       --input path2adata=${SC_DIR}/plate\${SLURM_ARRAY_TASK_ID}${SUFFIX}.h5ad \
               	       path2obs=${SC_DIR}/plate\${SLURM_ARRAY_TASK_ID}_obs.parquet \
                       path2var=${SC_DIR}/var.parquet \
		--output ${BULK_DIR}/plate\${SLURM_ARRAY_TASK_ID}${SUFFIX}.h5ad \
		--filter_malat1 \
        	--filter_low_counts \
        	--filter_nans"


echo "> Pseudobulking is done"

echo "> Running combining datasets"
sbatch -W -J pseudobulk_tahoe \
       -t 3:00:00 \
       -n 1 \
       --qos=${QOS} \
       --mem=250G \
       --partition=${PARTITION} \
       --cpus-per-task=2 \
       -e ${LOGS_DIR}/tahoe/${SUBDIR}/pseudobulk_tahoe_unite.%j.err \
       -o ${LOGS_DIR}/tahoe/${SUBDIR}/pseudobulk_tahoe_unite.%j.out \
       --wrap="eval \"\$(mamba shell hook --shell bash)\" && \
               mamba activate ${ENV_DIR} && \
	       python3 -m src.pseudobulking.common.run_combining_datasets \
	       --input ${BULK_DIR}/ \
	       --output ${OUTPUT_DIR}/tahoe${SUFFIX}.h5ad"

echo "> Combining datasets is done"
