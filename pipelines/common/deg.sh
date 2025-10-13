#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -n 1
#SBATCH --qos=cpu_normal
#SBATCH --mem=100G
#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=2

set -e

MODE=""
PAR=""
FILT=""
CONFIG=""
ARG_S=""
ARG_F=""
ENV_DIR=./venv

DEG_PARAMETERS=(
        "group_all_replicates"
        "separate_replicates"
)

for item in "$@"; do
	shift
	case "$item" in
		subsampling)   set -- "$@" "-s" ;;
		*)        set -- "$@" "$item"
	esac
done

while getopts ":sp:f:c:h" opt; do
        case $opt in
                h)
                        echo "Run: $0 [-s] [-h] [-f] -p (parameters: ${DEG_PARAMETERS[*]})"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -h Help option, default=false"
                        echo "  -f Min number of cells in pseudobulk to filter samples with the lower number"
                        echo "  -p Parameter for DEG pipeline, required"
			echo "  -c Path to the config file"
                        exit 0
                        ;;
                s)
                        MODE="subsampling"
                        ;;
                p)
                        PAR=$OPTARG
                        ;;
                f)
                        FILT=$OPTARG
                        ;;
		c)	CONFIG=$OPTARG
        esac
done

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
        echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2; exit 1
fi

if ! [[ "$FILT" =~ ^[0-9]+$ ]] && ! [ -z "$FILT" ]; then
   echo "Error: Not a number" >&2; exit 1
else
   if [[ "$FILT" =~ ^[0-9]+$ ]]; then
   	ARG_F="--min_cells $FILT"
   fi
fi


if [ "$MODE" = "subsampling" ]; then
	echo "> Work with a subsample"
	SUFFIX="_subsample"
	SUBDIR="subsample"
	ARG_S="--subsampling"

else
	echo "> Work with a full version"
	SUFFIX=""
	SUBDIR="full"
	ARG_S=""
fi

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}


echo "> Preprocess pseudobulk with a $PAR parameter"
if ! par_process=$(jq -e ".$PAR.par_process" $CONFIG); then
  echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

#python3 -m src.deg.run_processing_pseudobulk --config <(echo "$par_process" | jq ".")


echo "> Run DGE with a $PAR parameter"
if ! par_deg=$(jq -e ".$PAR.par_deg" $CONFIG); then
  echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

# Create temporary file instead of using process substitution
temp_config=$(mktemp)
echo "$par_deg" | jq "." > "$temp_config"

Rscript ./src/deg/run_deg.R $ARG_F $ARG_S --config "$temp_config"

# Clean up temporary file
rm "$temp_config"
echo "> DGE calculations are completed"
