#rm -rf ~/logs
mkdir -p ~/logs

if [ ! -d ~/.venv ]; then
    echo "create .venv"
    chmod u+x $SCRIPT_DIR/set_env.sh
    $SCRIPT_DIR/set_env.sh
fi

sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking_subsample.sh
#sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking.sh
#sbatch  -e ~/logs/pseudobulk_tahoe.%j.err -o ~/logs/pseudobulk_tahoe.%j.out ./pipelines/tahoe/tahoe_pseudobulking_subsample.sh
#sbatch  -e ~/logs/pseudobulk_tahoe.%j.err -o ~/logs/pseudobulk_tahoe.%j.out ./pipelines/tahoe/tahoe_pseudobulking.sh
