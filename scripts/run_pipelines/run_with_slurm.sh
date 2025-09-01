set -e

#rm -rf ~/logs
mkdir -p ~/logs

SCRIPT_DIR="$(dirname "$0")"

#Create environment
echo "create .venv"
chmod u+x $SCRIPT_DIR/set_env.sh
$SCRIPT_DIR/set_env.sh

#-------------------SCIPLEX-------------------
#sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking.sh

#-------------------TAHOE---------------------
./pipelines/tahoe/tahoe_pseudobulking_parallel.sh > ~/logs/pseudobulk_tahoe.out 2>~/logs/pseudobulk_tahoe.err &


##################--DEBUGGING--OPTIONS--##################


#-------------------SCIPLEX-------------------
#sbatch -e ~/logs/pseudobulk_sciplex.%j.err -o ~/logs/pseudobulk_sciplex.%j.out ./pipelines/sciplex/sciplex_pseudobulking_subsample.sh

#-------------------TAHOE---------------------
#./pipelines/tahoe/tahoe_pseudobulking_subsample_parallel.sh > ~/logs/pseudobulk_tahoe.out 2>~/logs/pseudobulk_tahoe.err &
