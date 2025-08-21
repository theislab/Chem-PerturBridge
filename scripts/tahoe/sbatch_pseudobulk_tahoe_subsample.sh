#!/bin/bash
#SBATCH -J pseudobulk_tahoe
#SBATCH -t 5:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=250G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=5

set -e

source ~/.venv/bin/activate 
SC_DIR=~/data/Tahoe/raw
BULK_DIR=~/data/Tahoe/pseudobulk_to_merge
OUTPUT_DIR=~/data/Tahoe/pseudobulk

echo "> Download dataset"
for i in $(seq 1 14)
do
	echo "> Download plate $i"

	python3 ./utils/download_datasets.py \
		--input "key_adata=2025-02-25/h5ad/plate${i}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad" \
			"key_obs=tahoe100/obs/plate${i}.parquet" \
			"key_var=tahoe100/var.parquet" \
		--output "path2adata=${SC_DIR}/plate${i}.h5ad" \
                	 "path2obs=${SC_DIR}/plate${i}_obs.parquet" \
                 	"path2var=${SC_DIR}/var.parquet" \
        	--subsampling
done

echo "> Running pseudobulking"

for i in $(seq 1 14)
do
	echo "> Processing plate $i"

	python3 ./utils/pseudobulking.py \
		--dataset tahoe \
		--input "path2adata=${SC_DIR}/plate${i}_subsample.h5ad" \
                        "path2obs=${SC_DIR}/plate${i}_obs.parquet" \
                        "path2var=${SC_DIR}/var.parquet" \
		--output "$BULK_DIR/subsample/plate${i}_subsample.h5ad" \
		--filter_malat1 \
        	--filter_low_counts \
        	--filter_nans
done
echo "> Pseudobulking is done"

echo "> Running combining datasets"
python3 ./utils/combining_datasets.py --input "$BULK_DIR/subsample/" --output "$OUTPUT_DIR/subsample/tahoe.h5ad"
echo "> Combining datasets is done"
