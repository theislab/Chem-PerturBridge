#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=100G
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=2

set -e

if [ "$1" = "subsampling" ]; then
	echo "> Work with a subsample"
	SUFFIX="_subsample"
	SUBDIR="subsample"
	ARG="--subsampling"

else	
	echo "> Work with a full version"
	SUFFIX=""
	SUBDIR="full"
	ARG=""
fi

ENV_DIR=./venv

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

echo "> Run DEG"
Rscript ./src/deg/run_deg.r --config $2 $ARG

echo "> DEG calculations are completed"

