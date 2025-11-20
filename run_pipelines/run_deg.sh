#!/bin/bash

set -e

LOGS_DIR=./logs
PIPELINE_NAME="deg"
MODE_S=False
MODE_J=False
MODE_Q=False
ARG_S=""
ARG_F=""
ARG_J=""
DATASET=""
PAR=""
FILT=""
VALID_CHOICES=(
        "sciplex"
        "tahoe"
        "l1000"
)

DEG_PARAMETERS=(
	"group_all_replicates"
	"separate_replicates"
)

mkdir -p $LOGS_DIR

while getopts ":sjd:p:f:qh" opt; do
  	case $opt in
		h)
			echo "Run: $0 [-s] [-j] [-h] [-f] [-q] -d (dataset: ${VALID_CHOICES[*]}) -p (parameters: ${DEG_PARAMETERS[*]})"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -j Parallel mode (array jobs per cell type), default=false"
                        echo "  -h Help option"
			echo "  -f (value <int>) Min number of cells in pseudobulk to filter samples with the lower number"
			echo "  -q Filter samples that did not pass quality control"
			echo "  -d (dataset <str>) Name of dataset to process, required"
			echo "  -p (parameter <str>) Parameter for DEG pipeline, required"
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
    		d)
        		DATASET=$OPTARG
      			;;
  	esac
done


if [[ " ${VALID_CHOICES[*]} " != *" $DATASET "* ]]; then
        echo "Error: -d must be set up and one of: ${VALID_CHOICES[*]}" >&2; exit 1
fi

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
        echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2; exit 1
fi

if ! [[ "$FILT" =~ ^[0-9]+$ ]] && ! [ -z "$FILT" ]; then
  	echo "Error: Not a number" >&2; exit 1
else
	if [[ "$FILT" =~ ^[0-9]+$ ]]; then
		ARG_F="-f $FILT"
	fi
fi


if [[ "$MODE_S" == "True" ]]; then
	SUBDIR="subsample"
	ARG_S="-s"
else
	SUBDIR="full"
fi

if [[ "$MODE_J" == "True" ]]; then
	ARG_J="-j"
fi

# Determine filter folder names
if [[ -z "$FILT" ]]; then
	FILTER_FOLDER="filter_min_cells_0"
else
	FILTER_FOLDER="filter_min_cells_${FILT}"
fi

# Determine QC folder name
if [[ "$MODE_Q" == "True" ]]; then
	QC_FOLDER="qc_true"
	ARG_Q="-q"
else
	QC_FOLDER="qc_false"
	ARG_Q=""
fi

mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}
./pipelines/common/deg.sh $ARG_S -p $PAR $ARG_F $ARG_Q $ARG_J \
	-c ./pipelines/${DATASET}/configs/deg/config.json \
	-d ${DATASET} \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_${DATASET}.PID$$.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}/deg_${DATASET}.PID$$.err &
