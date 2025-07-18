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

source /home/icb/olga.novitskaia/.env/bin/activate 

python3 ./utils/pseudobulking.py parse_args --input /home/icb/olga.novitskaia/data/Tahoe/subsamples/filtered/plate1.h5ad --output /home/icb/olga.novitskaia/data/Tahoe/subsamples/bulk/plate1.h5ad --select_cols sample cell_name bc1_well drugname_drugconc drug --min_cells 0 --drug_col drug --bulk_fields bc1_well cell_name

python3 ./utils/combining_replicates.py parse_args --input /home/icb/olga.novitskaia/data/Tahoe/subsamples/bulk --output /home/icb/olga.novitskaia/data/Tahoe/subsamples/cell_lines --min_cells 100 
