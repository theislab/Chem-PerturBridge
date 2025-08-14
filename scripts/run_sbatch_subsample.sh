#rm -rf logs
mkdir -p logs
sbatch tahoe/sbatch_pseudobulk_tahoe_subsample.sh
sbatch sciplex/sbatch_pseudobulk_sciplex_subsample.sh
