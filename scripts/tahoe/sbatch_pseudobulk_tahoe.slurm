#!/bin/bash
#SBATCH -J pseudobulk_tahoe
#SBATCH -o ./logs/pseudobulk_tahoe.%j.out
#SBATCH -e ./logs/pseudobulk_tahoe.%j.err
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=350G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=5

set -e

source /home/icb/olga.novitskaia/.env/bin/activate 
SC_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/filtered/"
BULK_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/bulk/"
CL_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/cell_lines/"

echo "> Running pseudobulking"
for file in $SC_DIR*; do
	filename=$(basename -- "$file")
	echo "> Processing $filename"
	python3 ./utils/pseudobulking.py --input $file --output $BULK_DIR$file --min_cells 50 --config ./configs/subsample/tahoe_pseudobulk_subsample.json
done
echo "> Pseudobulking is done"

echo "> Running combining replicates"
python3 ./utils/combining_replicates.py --input $BULK_DIR --output $CL_DIR --config ./configs/subsample/tahoe_combining_replicates_subsample.json
echo "> Combining replicates is done"
