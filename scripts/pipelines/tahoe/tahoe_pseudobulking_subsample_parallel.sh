#!/bin/bash
set -e

source ~/.venv/bin/activate 
SC_DIR=~/data/Tahoe/raw
BULK_DIR=~/data/Tahoe/pseudobulk_to_merge
OUTPUT_DIR=~/data/Tahoe/pseudobulk
N=14

echo "> Download dataset"

for i in $(seq 1 14)
do
	(
		echo "> Download plate $i"

		sbatch -W -J pseudobulk_tahoe \
			-t 10:00:00 \
			-n 1 \
			--qos=cpu_normal \
			--mem=250G \
			--partition=cpu_p \
			--cpus-per-task=2 \
			-e ~/logs/pseudobulk_tahoe_download.%j.err \
			-o ~/logs/pseudobulk_tahoe_download.%j.out \
			--wrap="source ~/.venv/bin/activate &&\
				python3 -m src.downloading.run_downloading_datasets \
				--input key_adata=2025-02-25/h5ad/plate${i}_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad \
					key_obs=tahoe100/obs/plate${i}.parquet \
					key_var=tahoe100/var.parquet \
				--output path2adata=${SC_DIR}/plate${i}.h5ad \
                			 path2obs=${SC_DIR}/plate${i}_obs.parquet \
			 		 path2var=${SC_DIR}/var.parquet \
				--subsampling"
	) &
    	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        	wait -n
    	fi
done

wait

echo "> Running pseudobulking"

for i in $(seq 1 14)
do
	(

		echo "> Processing plate $i"
		sbatch -W -J pseudobulk_tahoe \
                        -t 10:00:00 \
                        -n 1 \
                        --qos=cpu_normal \
                        --mem=250G \
                        --partition=cpu_p \
                        --cpus-per-task=2 \
			-e ~/logs/pseudobulk_tahoe_process.%j.err \
                        -o ~/logs/pseudobulk_tahoe_process.%j.out \
			--wrap="source ~/.venv/bin/activate &&\
				python3 -m src.pseudobulking.run_pseudobulking \
				--dataset tahoe \
				--input path2adata=${SC_DIR}/plate${i}_subsample.h5ad \
                        		path2obs=${SC_DIR}/plate${i}_obs.parquet \
                        		path2var=${SC_DIR}/var.parquet \
				--output $BULK_DIR/subsample/plate${i}_subsample.h5ad \
				--filter_malat1 \
        			--filter_low_counts \
        			--filter_nans"
	) &
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
                wait -n
        fi
done

wait

echo "> Pseudobulking is done"

echo "> Running combining datasets"
sbatch -W -J pseudobulk_tahoe \
                        -t 2:00:00 \
                        -n 1 \
                        --qos=cpu_normal \
                        --mem=250G \
                        --partition=cpu_p \
                        --cpus-per-task=2 \
			-e ~/logs/pseudobulk_tahoe_unite.%j.err \
                        -o ~/logs/pseudobulk_tahoe_unite.%j.out \
			--wrap="source ~/.venv/bin/activate &&\
				python3 -m src.pseudobulking.run_combining_datasets \
					--input $BULK_DIR/subsample/ \
					--output $OUTPUT_DIR/subsample/tahoe_subsample.h5ad"
echo "> Combining datasets is done"
