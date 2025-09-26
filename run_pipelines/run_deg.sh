#!/bin/bash

set -e

LOGS_DIR=./logs
ARG_S=""
ARG_F=""
DATASET=""
PAR=""
FILT=""
VALID_CHOICES=(
        "sciplex"
        "tahoe"
)

DEG_PARAMETERS=(
	"group_all_replicates"
	"separate_replicates"
)

mkdir -p $LOGS_DIR

while getopts ":sd:p:f:h" opt; do
  	case $opt in
		h)
			echo "Run: $0 [-s] [-h] [-f] -d (dataset: ${VALID_CHOICES[*]}) -p (parameters: ${DEG_PARAMETERS[*]})"
                        echo "  -s Subsample of a dataset for debugging, default=false"
                        echo "  -h Help option"
			echo "  -f (value <int>) Min number of cells in pseudobulk to filter samples with the lower number"
			echo "  -d (dataset <str>) Name of dataset to process, required"
			echo "  -p (parameter <str>) Parameter for DEG pipeline, required"
                        exit 0
                        ;;
    		s)
      			ARG_S="subsampling"
      			;;
		p)
			PAR=$OPTARG
			;;
		f)
			FILT=$OPTARG
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


if [[ "$ARG_S" == "subsampling" ]]; then
	SUBDIR="subsample"
else
	SUBDIR="full"

fi

mkdir -p ${LOGS_DIR}/${DATASET}/${SUBDIR}
sbatch -e ${LOGS_DIR}/${DATASET}/${SUBDIR}/deg_${DATASET}.%j.err \
	-o ${LOGS_DIR}/${DATASET}/${SUBDIR}/deg_${DATASET}.%j.out \
	./pipelines/common/deg.sh $ARG_S -p $PAR $ARG_F -c ./pipelines/${DATASET}/configs/${DATASET}_deg.json
