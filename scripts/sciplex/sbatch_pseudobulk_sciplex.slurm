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

source /home/icb/olga.novitskaia/.env/bin/activate 

python3 ./utils/pseudobulking.py parse_args --input /home/icb/olga.novitskaia/data/Sciplex/subsamples/raw/subsample.h5ad --output /home/icb/olga.novitskaia/data/Sciplex/subsamples/bulk/subsample.h5ad --dropna cell_line --drop_cols cancer organism tissue_type perturbation_type --min_cells 0 --drug_col perturbation --bulk_fields well plate cell_line replicate dose_value perturbation time
