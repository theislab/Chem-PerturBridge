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

source /home/icb/olga.novitskaia/.env/bin/activate 
SC_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/filtered/"
BULK_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/pseudobulk_to_merge/"
OUTPUT="/home/icb/olga.novitskaia/data/Tahoe/subsamples/pseudobulk/tahoe.h5ad"

echo "> Running pseudobulking"
for file in $SC_DIR*; do
	filename=$(basename -- "$file")
	echo "> Processing $filename"
	python3 ./utils/pseudobulking.py --input $file --output $BULK_DIR$filename --filter_cells_params min_counts=1 max_genes=5 --config ./configs/tahoe_pseudobulk.json
done
echo "> Pseudobulking is done"

echo "> Running combining datasets"
python3 ./utils/combining_datasets.py --input $BULK_DIR --output $OUTPUT --config ./configs/tahoe_combining_datasets.json
echo "> Combining datasets is done"
