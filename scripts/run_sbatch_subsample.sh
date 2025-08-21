#rm -rf ~/logs
mkdir -p ~/logs

if [ ! -d ~/.venv ]; then
    echo "create .venv"
    chmod u+x ./set_env.sh
    ./set_env.sh
fi

#sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./sciplex/sbatch_pseudobulk_sciplex_subsample.sh
#sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./sciplex/sbatch_pseudobulk_sciplex.sh
#sbatch  -e ~/logs/pseudobulk_tahoe.%j.err -o ~/logs/pseudobulk_tahoe.%j.out ./tahoe/sbatch_pseudobulk_tahoe_subsample.sh
sbatch  -e ~/logs/pseudobulk_tahoe.%j.err -o ~/logs/pseudobulk_tahoe.%j.out ./tahoe/sbatch_pseudobulk_tahoe.sh
