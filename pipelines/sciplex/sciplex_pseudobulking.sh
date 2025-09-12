#!/bin/bash
#SBATCH -J pseudobulk_sciplex
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=100G
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=2

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

NAME="srivatsan20_sciplex3"
SC_DIR=./data/sciplex/raw
BULK_DATA=./data/sciplex/pseudobulk/${SUBDIR}/${NAME}
ENV_DIR=./venv

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

echo "> Download dataset"
python3 -m src.downloading.run_downloading_datasets \
	--input "key_adata=scperturb/srivatsan20_sciplex3.h5ad" \
                "key_obs=scperturb/obs/srivatsan20_sciplex3.parquet" \
                "key_var=scperturb/var/srivatsan20_sciplex3.parquet" \
	--output "path2adata=${SC_DIR}/${NAME}.h5ad" \
		 "path2obs=${SC_DIR}/${NAME}_obs.parquet" \
		 "path2var=${SC_DIR}/${NAME}_var.parquet" ${ARG}

echo "> Running pseudobulking"
python3 -m src.pseudobulking.common.run_pseudobulking \
	--dataset_name sciplex \
	--input "path2adata=${SC_DIR}/${NAME}${SUFFIX}.h5ad" \
		"path2obs=${SC_DIR}/${NAME}_obs.parquet" \
		"path2var=${SC_DIR}/${NAME}_var.parquet" \
	--output "${BULK_DATA}${SUFFIX}.h5ad" \
	--filter_malat1 \
	--filter_low_counts \
	--filter_nans

echo "> Pseudobulking is done"
