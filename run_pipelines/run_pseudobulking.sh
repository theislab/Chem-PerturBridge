#!/bin/bash

set -e

LOGS_DIR=./logs
ARG=""
DATASET=""
VALID_CHOICES=(
        "sciplex"
        "tahoe"
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
        echo "Error: -d must be set up and one of: ${VALID_CHOICES[*]}"
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
else
    	echo "Invalid choice"
fi
