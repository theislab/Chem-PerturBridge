#!/bin/bash
#SBATCH -J pseudobulk_tahoe
#SBATCH -o pseudobulk_tahoe.%j.out
#SBATCH -e pseudobulk_tahoe.%j.err
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_short
#SBATCH --mem=20G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=1

set -e

source /home/icb/olga.novitskaia/.env/bin/activate 
SC_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/filtered/"
BULK_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/bulk/"
CL_DIR="/home/icb/olga.novitskaia/data/Tahoe/subsamples/cell_lines/"

echo "> Running pseudobulking"
for file in $SC_DIR*; do
	filename=$(basename -- "$file")
	echo "> Processing $filename"
	python3 ./utils/pseudobulking.py --input $file --output $BULK_DIR$file --config ./configs/tahoe_pseudobulk.json
done
echo "> Pseudobulking is done"

echo "> Running combining replicates"
python3 ./utils/combining_replicates.py --input $BULK_DIR --output $CL_DIR --config ./configs/tahoe_combining_replicates.json
echo "> Combining replicates is done"
