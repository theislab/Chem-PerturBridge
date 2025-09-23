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

mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
sbatch -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/deg_${DATASET}.%j.err \
	-o ${LOGS_DIR}/${DATASET}/${SUBDIR}/deg_${DATASET}.%j.out \
	./pipelines/common/deg.sh $ARG ./pipelines/${DATASET}/configs/deg.json

