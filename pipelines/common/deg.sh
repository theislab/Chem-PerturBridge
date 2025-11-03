#!/bin/bash

set -e

PIPELINE_NAME="deg"
MODE_S=False
MODE_J=False
MODE_Q=False
PAR=""
FILT=""
CONFIG=""
DATASET=""
ARG_S=""
ARG_F=""
ENV_DIR=./venv
LOGS_DIR=./logs

MAX_CONCURRENT=10
QOS=cpu_normal
PARTITION=cpu_p

DEG_PARAMETERS=(
        "group_all_replicates"
        "separate_replicates"
)

while getopts ":sjp:f:q:c:d:h" opt; do
        case $opt in
                h)
                        echo "Run: $0 [-s] [-j] [-h] [-f] [-q] -p (parameters: ${DEG_PARAMETERS[*]}) -c CONFIG -d DATASET"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -j Parallel mode (array jobs per cell type), default=false"
                        echo "  -h Help option, default=false"
                        echo "  -f Min number of cells in pseudobulk to filter samples with the lower number"
                        echo "  -q Filter samples that did not pass quality control"
                        echo "  -p Parameter for DEG pipeline, required"
			echo "  -c Path to the config file, required"
			echo "  -d Dataset name, required"
                        exit 0
                        ;;
                s)
                        MODE_S=True
                        ;;
                j)
                        MODE_J=True
                        ;;
                p)
                        PAR=$OPTARG
                        ;;
                f)
                        FILT=$OPTARG
                        ;;
		q)
			MODE_Q=True
			;;
		c)	CONFIG=$OPTARG
			;;
		d)	DATASET=$OPTARG
			;;
        esac
done

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
        echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2; exit 1
fi

# Validate DATASET
if [ -z "$DATASET" ]; then
  echo "Error: Dataset name must be provided via -d flag" >&2
  exit 1
fi

if ! [[ "$FILT" =~ ^[0-9]+$ ]] && ! [ -z "$FILT" ]; then
        echo "Error: Not a number" >&2; exit 1
else
        if [[ "$FILT" =~ ^[0-9]+$ ]]; then
   	        ARG_F="--min_cells $FILT"
        fi
fi

# Determine filter folder name
if [[ -z "$FILT" ]]; then
	FILTER_FOLDER="filter_min_cells_0"
else
	FILTER_FOLDER="filter_min_cells_${FILT}"
fi

# Set QC argument if flag is provided
if [[ "$MODE_Q" == "True" ]]; then
        echo "> Filter samples that did not pass quality control"
        QC_FOLDER="qc_true"
        ARG_Q="--qc"
else
        echo "> Do not filter samples that did not pass quality control"
        QC_FOLDER="qc_false"
        ARG_Q=""
fi


if [ "$MODE_S" = "True" ]; then
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

# Validate CONFIG
if [ -z "$CONFIG" ]; then
  echo "Error: Config path must be provided via -c flag" >&2
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "Error: Config file not found: $CONFIG" >&2
  exit 1
fi

echo "> Using config: $CONFIG"

echo "> Preprocess pseudobulk with a $PAR parameter"
if ! par_process=$(jq -e ".$PAR.par_process" $CONFIG); then
  echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

# Add split_by_celltype flag if parallel mode
if [ "$MODE_J" = "True" ]; then
        echo "> Mode: PARALLEL (array jobs per cell type)"
        ARG_J="--split_by_celltype"
else
        echo "> Mode: SEQUENTIAL (single job, all cell types)"
        ARG_J=""
fi


eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

# Create tmp directory for configs
TMP_DIR="./tmp"
TMP_DIR_DEG="${TMP_DIR}/deg_$$"
# Clean up if it exists from a previous failed run
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"
mkdir -p "${TMP_DIR_DEG}"
# Set up automatic cleanup on exit (even if script fails)
trap "rm -rf ${TMP_DIR}" EXIT

# Create config file in tmp directory
preprocess_config="${TMP_DIR_DEG}/preprocess_config.json"
echo "$par_process" | jq "." > "$preprocess_config"

sbatch -W -J deg_processing_pseudobulk \
       --partition=${PARTITION} \
       --qos=${QOS} \
       --mem=250G \
       --time=2:00:00 \
       --cpus-per-task=2 \
       -o "${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_processing_pseudobulk.%j.out" \
       -e "${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_processing_pseudobulk.%j.err" \
       --wrap="eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                python3 -m src.deg.run_processing_pseudobulk $ARG_J \
                --config ${preprocess_config}"

# Keep tmp directory for DEG step (cleanup at the end)


echo "> Run DGE with a $PAR parameter"
if ! par_deg=$(jq -e ".$PAR.par_deg" $CONFIG); then
        echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

# Get input directory from config
INPUT_DIR=$(echo "$par_deg" | jq -r '.input_dir')

# Set mode-specific parameters now that INPUT_DIR is defined
if [ "$MODE_J" = "True" ]; then
	LOG_PREFIX="deg_celltype"
	SEARCH_DIR="${INPUT_DIR}/by_celltype"
	FILE_PATTERN="*"
        MEM=150G
else
	LOG_PREFIX="deg_sequential"
	SEARCH_DIR="${INPUT_DIR}"
	FILE_PATTERN="*.h5ad"
        MEM=250G
fi

# Check if search directory exists
if [ ! -d "$SEARCH_DIR" ]; then
	echo "Error: Directory not found: $SEARCH_DIR" >&2
	exit 1
fi

# List files to process
mapfile -t INPUT_FILES < <(ls "$SEARCH_DIR"/$FILE_PATTERN 2>/dev/null | sort)
N_FILES=${#INPUT_FILES[@]}

if [ "$N_FILES" -eq 0 ]; then
	echo "Error: No files found in $SEARCH_DIR" >&2
	exit 1
fi

echo "  Found $N_FILES files to process"

# Create file list in tmp directory
FILE_LIST="${TMP_DIR_DEG}/file_list.txt"
printf '%s\n' "${INPUT_FILES[@]}" > "$FILE_LIST"

# Create DEG config in tmp directory
deg_config="${TMP_DIR_DEG}/deg_config.json"
echo "$par_deg" | jq "." > "$deg_config"

# Submit array job with inline command
echo "  Submitting array job [${N_FILES} jobs, max ${MAX_CONCURRENT} concurrent]..."
sbatch -W \
	-J deg_analysis \
	--array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
	--partition=${PARTITION} \
	--qos=${QOS} \
	--mem=${MEM} \
	--time=6:00:00 \
	--cpus-per-task=2 \
	--output="${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_analysis_%A_%a.out" \
	--error="${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_analysis_%A_%a.err" \
	--wrap="eval \"\$(mamba shell hook --shell bash)\" && \
		mamba activate ${ENV_DIR} && \
		INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${FILE_LIST}) && \
		Rscript ./src/deg/run_deg.R \
			--input \"\$INPUT_FILE\" \
			--config ${deg_config} \
			$ARG_F $ARG_S $ARG_Q"
rm -rf "${TMP_DIR}"
echo "> DGE calculations are completed"
