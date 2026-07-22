#!/bin/bash

set -e

LOGS_DIR=./logs
PIPELINE_NAME="deg"
MODE_S=False
MODE_J=False
MODE_Q=False
MODE_N=False
MODE_G=False
MODE_P=False
MODE_C=False
ARG_S=""
ARG_F=""
ARG_J=""
ARG_N=""
ARG_G=""
ARG_P=""
ARG_C=""
DATASET=""
PAR=""
FILT=""
VALID_CHOICES=(
        "sciplex"
        "tahoe"
        "dilimap_train"
        "dilimap_train_val"
        "l1000_phase1"
        "l1000_phase2"
        "op3"
        "novartis"
        "vcpi_0001"
        "vcpi_0002"
        "vcpi_0003"
        "gdpx2"
        "cigs_mce"
        "cigs_tcm"
)

DEG_PARAMETERS=(
	"group_all_replicates"
	"separate_replicates"
)

mkdir -p $LOGS_DIR

while getopts ":sjd:p:f:qnhgPC" opt; do
  	case $opt in
		h)
			echo "Run: $0 [-s] [-j] [-h] [-g] [-P] [-C] [-f] [-q] [-n] -d (dataset: ${VALID_CHOICES[*]}) -p (parameters: ${DEG_PARAMETERS[*]})"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -j Parallel mode (array jobs per cell type), default=false"
                        echo "  -g Use GPU QOS/partition (gpu_normal/gpu_p), default=CPU"
                        echo "  -P Use cpu_preemptible QOS for DEG step (100G, 48h, 10 CPUs), default=false"
                        echo "  -C Per-plate controls: pass --per_plate_control to preprocessing and R DEG"
                        echo "  -h Help option"
			echo "  -f (value <int>) Min number of cells in pseudobulk to filter samples with the lower number"
			echo "  -q Filter samples that did not pass quality control"
			echo "  -n Normalized data: preprocessing + R skip expression normalization; Ensembl drops ambiguous targets"
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
		n)
			MODE_N=True
      			;;
		g)
			MODE_G=True
			;;
		P)
			MODE_P=True
			;;
		C)
			MODE_C=True
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
	ARG_S=""
fi

if [[ "$MODE_J" == "True" ]]; then
	ARG_J="-j"
else
	ARG_J=""
fi

if [[ "$MODE_N" == "True" ]]; then
	ARG_N="-n"
else
	ARG_N=""
fi

if [[ "$MODE_G" == "True" ]]; then
	ARG_G="-g"
else
	ARG_G=""
fi

if [[ "$MODE_P" == "True" ]]; then
	ARG_P="-P"
else
	ARG_P=""
fi

if [[ "$MODE_C" == "True" ]]; then
	ARG_C="-C"
	PAR_FOLDER="${PAR}_per_plate_control"
else
	ARG_C=""
	PAR_FOLDER="${PAR}"
fi

if [[ -z "$FILT" ]]; then
	FILTER_FOLDER="filter_min_cells_0"
else
	FILTER_FOLDER="filter_min_cells_${FILT}"
fi

if [[ "$MODE_Q" == "True" ]]; then
	QC_FOLDER="qc_true"
	ARG_Q="-q"
else
	QC_FOLDER="qc_false"
	ARG_Q=""
fi

mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR_FOLDER}/${QC_FOLDER}/${FILTER_FOLDER}
./pipelines/common/deg.sh $ARG_S -p $PAR $ARG_F $ARG_Q $ARG_J $ARG_N $ARG_G $ARG_P $ARG_C \
	-c ./pipelines/${DATASET}/configs/deg/config.json \
	-d ${DATASET} \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR_FOLDER}/${QC_FOLDER}/${FILTER_FOLDER}/deg_${DATASET}.PID$$.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR_FOLDER}/${QC_FOLDER}/${FILTER_FOLDER}/deg_${DATASET}.PID$$.err &
