#!/bin/bash

set -e

LOGS_DIR=./logs
ARG=""
DATASET=""
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
        "gdpx2"
        "cigs_mce"
        "cigs_tcm"
)

mkdir -p $LOGS_DIR

while getopts ":sd:h" opt; do
  	case $opt in
		h)
			echo "Run: $0 [-s] [-h] -d (dataset: ${VALID_CHOICES[*]})"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -h Help option, default=false"
			echo "  -d (dataset)  Name of dataset to process, required"
                        exit 0
                        ;;
    		s)
      			ARG="subsampling"
      			;;
    		d)
        		DATASET=$OPTARG
      			;;
  	esac
done


if [[ " ${VALID_CHOICES[*]} " != *" $DATASET "* ]]; then
        echo "Error: -d must be set up and one of: ${VALID_CHOICES[*]}" >&2
        exit 1
fi

if [[ "$ARG" == "subsampling" ]]; then
	SUBDIR="subsample"
else
	SUBDIR="full"

fi

if [[ "$DATASET" == "sciplex" ]]; then
	mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
  	sbatch -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_sciplex.%j.err \
		-o ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_sciplex.%j.out \
		./pipelines/sciplex/sciplex_pseudobulking.sh $ARG

elif [[ "$DATASET" == "tahoe" ]]; then
	mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
    	./pipelines/tahoe/tahoe_pseudobulking_parallel.sh $ARG \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_tahoe.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_tahoe.err &

elif [[ "$DATASET" == "dilimap_train" ]]; then
        mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
        ./pipelines/dilimap_train/dilimap_train_pseudobulking.sh $ARG \
                >  ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_dilimap_train.out \
                2> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_dilimap_train.err &

elif [[ "$DATASET" == "dilimap_train_val" ]]; then
        mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
        ./pipelines/dilimap_train_val/dilimap_train_val_pseudobulking.sh $ARG \
                >  ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_dilimap_train_val.out \
                2> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_dilimap_train_val.err &
		
elif [[ "$DATASET" == "l1000_phase1" ]] || [[ "$DATASET" == "l1000_phase2" ]] || [[ "$DATASET" == "op3" ]] \
	|| [[ "$DATASET" == "vcpi_0001" ]] || [[ "$DATASET" == "vcpi_0002" ]] \
	|| [[ "$DATASET" == "cigs_mce" ]] || [[ "$DATASET" == "cigs_tcm" ]]; then
	mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
    	./pipelines/${DATASET}/${DATASET}_pseudobulking.sh $ARG \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.err &

elif [[ "$DATASET" == "novartis" ]]; then
	mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
	./pipelines/${DATASET}/${DATASET}_pseudobulking.sh $ARG \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.err &

elif [[ "$DATASET" == "gdpx2" ]]; then
	mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
	./pipelines/${DATASET}/${DATASET}_pseudobulking.sh $ARG \
		> ${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.out \
		2>${LOGS_DIR}/${DATASET}/${SUBDIR}/pseudobulk_${DATASET}.PID$$.err &
else
    	echo "Invalid choice"
fi
