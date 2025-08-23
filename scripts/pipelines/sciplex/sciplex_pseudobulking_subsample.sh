#!/bin/bash
#SBATCH -J pseudobulk_sciplex
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=20G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=5

set -e

source ~/.venv/bin/activate
NAME="srivatsan20_sciplex3"
SC_DIR=~/data/Sciplex/raw
BULK_DATA=~/data/Sciplex/pseudobulk/${NAME}

echo "> Download dataset"
python3 -m src.downloading.run_downloading_datasets \
	--output "path2adata=${SC_DIR}/${NAME}.h5ad" \
		 "path2obs=${SC_DIR}/${NAME}_obs.parquet" \
		 "path2var=${SC_DIR}/${NAME}_var.parquet" \
	--config ./pipelines/sciplex/configs/sciplex_downloading.json \
	--subsampling

echo "> Running pseudobulking"
python3 -m src.pseudobulking.run_pseudobulking \
	--dataset sciplex \
	--input "path2adata=${SC_DIR}/${NAME}_subsample.h5ad" \
		"path2obs=${SC_DIR}/${NAME}_obs.parquet" \
		"path2var=${SC_DIR}/${NAME}_var.parquet" \
	--output "${BULK_DATA}_subsample.h5ad" \
	--filter_malat1 \
	--filter_low_counts \
	--filter_nans
echo "> Pseudobulking is done"
