#!/bin/bash
#SBATCH -J pseudobulk_sciplex
#SBATCH -o pseudobulk_sciplex.%j.out
#SBATCH -e pseudobulk_sciplex.%j.err
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_short
#SBATCH --mem=20G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=1

set -e

source /home/icb/olga.novitskaia/.env/bin/activate
INPUT="/home/icb/olga.novitskaia/data/Sciplex/subsamples/raw/subsample.h5ad"
OUTPUT="/home/icb/olga.novitskaia/data/Sciplex/subsamples/bulk/subsample.h5ad" 

echo "> Running pseudobulking"
python3 ./utils/pseudobulking.py --input $INPUT --output $OUTPUT --config ./configs/sciplex_pseudobulk.json
echo "> Pseudobulking is done"
