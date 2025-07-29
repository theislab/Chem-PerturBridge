#!/bin/bash
#SBATCH -J pseudobulk_sciplex
#SBATCH -o ./logs/pseudobulk_sciplex.%j.out
#SBATCH -e ./logs/pseudobulk_sciplex.%j.err
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=100G
#SBATCH --partition=cpu_p 
#SBATCH --cpus-per-task=5

set -e

source /home/icb/olga.novitskaia/.env/bin/activate 
INPUT="/lustre/groups/ml01/workspace/manuel.gander/sciplex/SrivatsanTrapnell2020_sciplex3.h5ad"
OUTPUT="/home/icb/olga.novitskaia/data/Sciplex/full/pseudobulk/sciplex3.h5ad"

echo "> Running pseudobulking"
python3 ./utils/pseudobulking.py --input $INPUT --output $OUTPUT --config ./configs/sciplex_pseudobulk.json
echo "> Pseudobulking is done"
