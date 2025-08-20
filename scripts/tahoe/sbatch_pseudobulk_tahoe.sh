#!/bin/bash
#SBATCH -J pseudobulk_tahoe
#SBATCH -o ./logs/pseudobulk_tahoe.%j.out
#SBATCH -e ./logs/pseudobulk_tahoe.%j.err
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_short
#SBATCH --mem=20G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=1

set -e

#source /home/icb/olga.novitskaia/.env/bin/activate 
SC_DIR="/home/ubuntu/data/Tahoe/raw"
BULK_DIR="/home/ubuntu/data/Tahoe/pseudobulk_to_merge"
OUTPUT_DIR="/home/ubuntu/data/Tahoe/pseudobulk"

echo "> Download dataset"
for i in $(seq 1 14)
do
	python3 ./utils/download_datasets.py \
		--input "key_adata=2025-02-25/h5ad/plate${i}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad" \
			"key_obs=tahoe100/obs/plate${i}.parquet" \
			"key_var=tahoe100/var.parquet" \
		--output "path2adata=${SC_DIR}/plate${i}.h5ad" \
                	 "path2obs=${SC_DIR}/plate${i}_obs.parquet" \
                 	"path2var=${SC_DIR}/var.parquet"
done

echo "> Running pseudobulking"

for i in $(seq 1 14)
do
	echo "> Processing plate $i"

	python3 ./utils/pseudobulking.py \
		--dataset tahoe \
		--input "path2adata=${SC_DIR}/plate${i}.h5ad" \
                        "path2obs=${SC_DIR}/plate${i}_obs.parquet" \
                        "path2var=${SC_DIR}/var.parquet" \
		--output "$BULK_DIR/full/plate${i}.h5ad" \
		--filter_malat1 \
        	--filter_low_counts \
        	--filter_nans
done
echo "> Pseudobulking is done"

echo "> Running combining datasets"
python3 ./utils/combining_datasets.py --input "$BULK_DIR/full/" --output "$OUTPUT_DIR/full/tahoe.h5ad"
echo "> Combining datasets is done"
