#!/bin/bash

set -e

# Get project root directory (parent of pipelines directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

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
N_OBS_JSON_FILE_NAME="n_obs_table.json"

# Threshold for choosing parallel vs sequential DEG
N_OBS_THRESHOLD=5000

MAX_CONCURRENT=10
QOS=cpu_preemptible
PARTITION=cpu_p

# Number of batches for step 2 (parallel processing) - hardcoded
N_BATCHES=20

DEG_PARAMETERS=(
        "group_all_replicates"
        "separate_replicates"
)

while getopts ":sjqp:f:c:d:h" opt; do
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
                c)
                        CONFIG=$OPTARG
                        ;;
                d)
                        DATASET=$OPTARG
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

# Add split_by_celltype flag if parallel mode (needed early for preprocessing)
if [ "$MODE_J" = "True" ]; then
        echo "> Mode: PARALLEL BY CELL TYPE (array jobs per cell type)"
        ARG_J="--split_by_celltype"
else
        echo "> Mode: SEQUENTIAL BY CELL TYPE (single job, all cell types)"
        ARG_J=""
fi

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

# Create tmp directory for configs
TMP_DIR="${PROJECT_ROOT}/tmp"
TMP_DIR_DEG="${TMP_DIR}/deg_$$"
# Clean up if it exists from a previous failed run
mkdir -p "${TMP_DIR}"
mkdir -p "${TMP_DIR_DEG}"
# Set up automatic cleanup on exit (even if script fails)
trap "rm -rf ${TMP_DIR_DEG}" EXIT

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

# Create log base directory (needed for preprocessing logs)
LOGS_BASE_DIR="${LOGS_DIR}/${DATASET}/${SUBDIR}/${PIPELINE_NAME}/${PAR}/${QC_FOLDER}/${FILTER_FOLDER}"
mkdir -p "${LOGS_BASE_DIR}/preprocessing"

# ============================================================================
# STEP 1: PREPROCESS PSEUDOBULK
# ============================================================================

echo "> Preprocess pseudobulk with a $PAR parameter"
if ! par_process=$(jq -e ".$PAR.par_process" $CONFIG); then
  echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

# Create config file in tmp directory
preprocess_config="${TMP_DIR_DEG}/preprocess_config.json"
echo "$par_process" | jq "." > "$preprocess_config"

echo "> Running preprocessing pseudobulk"
sbatch -W -J deg_processing_pseudobulk \
       --partition=${PARTITION} \
       --qos=${QOS} \
       --mem=50G \
       --time=2:00:00 \
       --cpus-per-task=2 \
       --output="${LOGS_BASE_DIR}/preprocessing/deg_processing_pseudobulk.%j.out" \
       --error="${LOGS_BASE_DIR}/preprocessing/deg_processing_pseudobulk.%j.err" \
       --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                python3 -m src.deg.run_processing_pseudobulk ${ARG_J} \
                --config ${preprocess_config}"

echo "> Preprocessing pseudobulk is done"

# ============================================================================
# STEP 2: DGE ANALYSIS (Parallel or Sequential based on n_obs)
# ============================================================================

echo "> Starting DGE analysis with a $PAR parameter"
if ! par_deg=$(jq -e ".$PAR.par_deg" $CONFIG); then
        echo "Error: Failed to extract parameters from $CONFIG" >&2; exit 1
fi

# Get input directory from config
INPUT_DIR=$(echo "$par_deg" | jq -r '.input_dir')

# Set mode-specific parameters now that INPUT_DIR is defined
if [ "$MODE_J" = "True" ]; then
        
	LOG_PREFIX="deg_celltype"
	SEARCH_DIR="${INPUT_DIR}/by_celltype"
	FILE_PATTERN="*.h5ad"
        MEM=50G
else
	LOG_PREFIX="deg_sequential"
	SEARCH_DIR="${INPUT_DIR}"
	FILE_PATTERN="*.h5ad"
        MEM=100G
fi

# Check if search directory exists
if [ ! -d "$SEARCH_DIR" ]; then
	echo "Error: Directory not found: $SEARCH_DIR" >&2
	exit 1
fi

# List files to process (only .h5ad files, exclude .json and other files)
mapfile -t INPUT_FILES < <(ls "$SEARCH_DIR"/$FILE_PATTERN 2>/dev/null | grep -E '\.h5ad$' | sort)
N_FILES=${#INPUT_FILES[@]}

if [ "$N_FILES" -eq 0 ]; then
	echo "Error: No files found in $SEARCH_DIR" >&2
	exit 1
fi

echo "  Found $N_FILES files to process"

# Create DEG config in tmp directory
deg_config="${TMP_DIR_DEG}/deg_config.json"
echo "$par_deg" | jq "." > "$deg_config"

# Get output directory from config
OUTPUT_DIR=$(echo "$par_deg" | jq -r '.output_dir')

# Determine full output directory path
FULL_OUTPUT_DIR="${OUTPUT_DIR}/${SUBDIR}/${QC_FOLDER}/${FILTER_FOLDER}/results"

# Intermediate directory is now a subdirectory of output_dir
FULL_INTERMEDIATE_DIR="${OUTPUT_DIR}/${SUBDIR}/${QC_FOLDER}/${FILTER_FOLDER}/intermediate"

# Create log directories for each step (additional directories for DEG steps)
mkdir -p "${LOGS_BASE_DIR}/step1_filter"
mkdir -p "${LOGS_BASE_DIR}/step2_voom"
mkdir -p "${LOGS_BASE_DIR}/step3_lmfit"
mkdir -p "${LOGS_BASE_DIR}/step4_batch"
mkdir -p "${LOGS_BASE_DIR}/step5_aggregate"
mkdir -p "${LOGS_BASE_DIR}/deg_sequential"

# ============================================================================
# Function: Run Parallel DEG (Step 1, 2, 3, 4, 5)
# Steps 1-3: Array jobs across all cell lines (parallel)
# Step 4: Sequential per cell line (with internal batch parallelization)
# Step 5: Sequential per cell line
# ============================================================================

run_parallel_deg() {
    local FILES_ARRAY_NAME=$1
    local N_FILES=$2
    
    # Get the array by name
    local -n FILES_ARRAY=$FILES_ARRAY_NAME
    
    echo ""
    echo "============================================"
    echo "PARALLEL BY CONTRASTS (Step 1, 2, 3, 4, 5)"
    echo "============================================"
    echo "  Processing ${N_FILES} cell line(s)"
    
    # Create directories
    STEP1_OUTPUT_DIR="${FULL_INTERMEDIATE_DIR}/step1_filter"
    mkdir -p "${STEP1_OUTPUT_DIR}"
    STEP2_OUTPUT_DIR="${FULL_INTERMEDIATE_DIR}/step2_voom"
    mkdir -p "${STEP2_OUTPUT_DIR}"
    STEP3_OUTPUT_DIR="${FULL_INTERMEDIATE_DIR}/step3_lmfit"
    mkdir -p "${STEP3_OUTPUT_DIR}"
    STEP5_OUTPUT_DIR="${FULL_OUTPUT_DIR}"
    
    # Create file lists for array jobs
    STEP1_FILE_LIST="${TMP_DIR_DEG}/step1_file_list.txt"
    printf '%s\n' "${FILES_ARRAY[@]}" > "$STEP1_FILE_LIST"
    
    # ============================================================================
    # STEP 1: FILTER - Array job across all cell lines (parallel)
    # ============================================================================
    echo ""
    echo "  STEP 1: FILTER (filterByExpr + calcNormFactors) - Array job [${N_FILES} cell lines]"
    STEP1_ARRAY_JOB_ID=$(sbatch -W --parsable \
        -J deg_step1_all \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=50G \
        --time=24:00:00 \
        --cpus-per-task=4 \
        --output="${LOGS_BASE_DIR}/step1_filter/step1.%A_%a.out" \
        --error="${LOGS_BASE_DIR}/step1_filter/step1.%A_%a.err" \
        --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${STEP1_FILE_LIST}) && \
                INPUT_BASENAME=\$(basename \"\$INPUT_FILE\" .h5ad | sed 's/_processed//') && \
                Rscript ./src/deg/run_deg_parallel_contrasts_step1_filter.R \
                    --input \"\$INPUT_FILE\" \
                    --output_dir \"${STEP1_OUTPUT_DIR}\" \
                    --config \"${deg_config}\" \
                    $ARG_F $ARG_S $ARG_Q")
    
    # ============================================================================
    # STEP 2: VOOM - Array job across all cell lines (parallel, depends on step1)
    # ============================================================================
    echo ""
    echo "  STEP 2: VOOM (voom transformation) - Array job [${N_FILES} cell lines]"
    STEP2_FILE_LIST="${TMP_DIR_DEG}/step2_file_list.txt"
    # Create list of filtered file paths
    for INPUT_FILE in "${FILES_ARRAY[@]}"; do
        INPUT_BASENAME=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        echo "${STEP1_OUTPUT_DIR}/${INPUT_BASENAME}_filtered.rds" >> "$STEP2_FILE_LIST"
    done
    
    STEP2_ARRAY_JOB_ID=$(sbatch -W --parsable \
        -J deg_step2_all \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=50G \
        --time=24:00:00 \
        --cpus-per-task=4 \
        --output="${LOGS_BASE_DIR}/step2_voom/step2.%A_%a.out" \
        --error="${LOGS_BASE_DIR}/step2_voom/step2.%A_%a.err" \
        --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${STEP2_FILE_LIST}) && \
                INPUT_BASENAME=\$(basename \"\$INPUT_FILE\" _filtered.rds) && \
                Rscript ./src/deg/run_deg_parallel_contrasts_step2_voom.R \
                    --input_file \"\$INPUT_FILE\" \
                    --output_dir \"${STEP2_OUTPUT_DIR}\"")
    
    # ============================================================================
    # STEP 3: LMFIT - Array job across all cell lines (parallel, depends on step2)
    # ============================================================================
    echo ""
    echo "  STEP 3: LMFIT (linear model fitting) - Array job [${N_FILES} cell lines]"
    STEP3_FILE_LIST="${TMP_DIR_DEG}/step3_file_list.txt"
    # Create list of voom file paths
    for INPUT_FILE in "${FILES_ARRAY[@]}"; do
        INPUT_BASENAME=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        echo "${STEP2_OUTPUT_DIR}/${INPUT_BASENAME}_voom.rds" >> "$STEP3_FILE_LIST"
    done
    
    STEP3_ARRAY_JOB_ID=$(sbatch -W --parsable \
        -J deg_step3_all \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=50G \
        --time=24:00:00 \
        --cpus-per-task=4 \
        --output="${LOGS_BASE_DIR}/step3_lmfit/step3.%A_%a.out" \
        --error="${LOGS_BASE_DIR}/step3_lmfit/step3.%A_%a.err" \
        --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${STEP3_FILE_LIST}) && \
                INPUT_BASENAME=\$(basename \"\$INPUT_FILE\" _voom.rds) && \
                Rscript ./src/deg/run_deg_parallel_contrasts_step3_lmfit.R \
                    --input_file \"\$INPUT_FILE\" \
                    --output_dir \"${STEP3_OUTPUT_DIR}\"")
    
    # ============================================================================
    # STEP 4: BATCH PROCESSING - Sequential per cell line (with internal parallelization)
    # ============================================================================
    echo ""
    echo "  STEP 4: BATCH PROCESSING - Sequential per cell line [${N_BATCHES} batches per cell line]"
    echo "  Step 3 completed, starting Step 4 sequentially..."
    
    for i in "${!FILES_ARRAY[@]}"; do
        INPUT_FILE="${FILES_ARRAY[$i]}"
        CELL_TYPE=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        INPUT_BASENAME=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        STEP4_OUTPUT_DIR="${FULL_INTERMEDIATE_DIR}/step4_batch/${CELL_TYPE}"
        mkdir -p "${STEP4_OUTPUT_DIR}"
        mkdir -p "${LOGS_BASE_DIR}/step4_batch/${CELL_TYPE}"
        
        echo "    [$((i+1))/${N_FILES}] Processing ${CELL_TYPE}..."
        
        STEP4_JOB_ID=$(sbatch -W --parsable \
            -J deg_step4_${CELL_TYPE} \
            --array=0-$((N_BATCHES-1))%${MAX_CONCURRENT} \
            --partition=${PARTITION} \
            --qos=${QOS} \
            --mem=50G \
            --time=24:00:00 \
            --cpus-per-task=4 \
            --output="${LOGS_BASE_DIR}/step4_batch/${CELL_TYPE}/step4_${CELL_TYPE}.%A_%a.out" \
            --error="${LOGS_BASE_DIR}/step4_batch/${CELL_TYPE}/step4_${CELL_TYPE}.%A_%a.err" \
            --wrap="export TMPDIR=\${HOME}/tmp && \
                    mkdir -p \${TMPDIR} && \
                    cd ${PROJECT_ROOT} && \
                    eval \"\$(mamba shell hook --shell bash)\" && \
                    mamba activate ${ENV_DIR} && \
                    Rscript ./src/deg/run_deg_parallel_contrasts_step4_batch.R \
                        --input_file \"${STEP3_OUTPUT_DIR}/${INPUT_BASENAME}.rds\" \
                        --output_dir \"${STEP4_OUTPUT_DIR}\" \
                        --batch_id \${SLURM_ARRAY_TASK_ID} \
                        --n_batches ${N_BATCHES}")
    done
    
    # ============================================================================
    # STEP 5: AGGREGATE - Array job across all cell lines (parallel)
    # ============================================================================
    echo ""
    echo "  STEP 5: AGGREGATE - Array job [${N_FILES} cell lines]"
    echo "  All Step 4 jobs completed, starting Step 5..."
    
    # Create file list for step 5 array job
    # Format: CELL_TYPE|INPUT_BASENAME
    STEP5_FILE_LIST="${TMP_DIR_DEG}/step5_file_list.txt"
    for i in "${!FILES_ARRAY[@]}"; do
        INPUT_FILE="${FILES_ARRAY[$i]}"
        CELL_TYPE=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        INPUT_BASENAME=$(basename "$INPUT_FILE" .h5ad | sed 's/_processed//')
        echo "${CELL_TYPE}|${INPUT_BASENAME}" >> "$STEP5_FILE_LIST"
    done
    
    STEP5_ARRAY_JOB_ID=$(sbatch -W --parsable \
        -J deg_step5_all \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=50G \
        --time=10:00:00 \
        --cpus-per-task=2 \
        --output="${LOGS_BASE_DIR}/step5_aggregate/step5.%A_%a.out" \
        --error="${LOGS_BASE_DIR}/step5_aggregate/step5.%A_%a.err" \
        --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                LINE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${STEP5_FILE_LIST}) && \
                CELL_TYPE=\$(echo \"\$LINE\" | cut -d'|' -f1) && \
                INPUT_BASENAME=\$(echo \"\$LINE\" | cut -d'|' -f2) && \
                STEP4_OUTPUT_DIR=\"${FULL_INTERMEDIATE_DIR}/step4_batch/\${CELL_TYPE}\" && \
                Rscript ./src/deg/run_deg_parallel_contrasts_step5_aggregate.R \
                    --input_file \"${STEP3_OUTPUT_DIR}/\${INPUT_BASENAME}.rds\" \
                    --input_dir \"\${STEP4_OUTPUT_DIR}\" \
                    --output_file \"${STEP5_OUTPUT_DIR}/\${CELL_TYPE}_de.h5ad\"")
    
    echo "  Step 5 Array Job ID: ${STEP5_ARRAY_JOB_ID}"
    
    echo ""
    echo "  ✓ Parallel DEG pipeline completed for all ${N_FILES} cell line(s)"
}

# ============================================================================
# Function: Run Sequential DEG
# ============================================================================

run_sequential_deg() {
    local FILES_ARRAY_NAME=$1
    local N_FILES=$2
    
    # Get the array by name
    local -n FILES_ARRAY=$FILES_ARRAY_NAME
    
    echo ""
    echo "============================================"
    echo "SEQUENTIAL BY CONTRASTS : Array Job Mode"
    echo "============================================"
    echo "  Processing ${N_FILES} file(s) with array job..."
    
    # Create file list in tmp directory
    SEQ_FILE_LIST="${TMP_DIR_DEG}/sequential_file_list.txt"
    printf '%s\n' "${FILES_ARRAY[@]}" > "$SEQ_FILE_LIST"
    
    echo "  Submitting array job [${N_FILES} jobs, max ${MAX_CONCURRENT} concurrent]..."
    sbatch -W \
        -J deg_sequential_analysis \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=${MEM} \
        --time=10:00:00 \
        --cpus-per-task=2 \
        --output="${LOGS_BASE_DIR}/deg_sequential/deg_sequential.%A_%a.out" \
        --error="${LOGS_BASE_DIR}/deg_sequential/deg_sequential.%A_%a.err" \
        --wrap="export TMPDIR=\${HOME}/tmp && \
                mkdir -p \${TMPDIR} && \
                cd ${PROJECT_ROOT} && \
                eval \"\$(mamba shell hook --shell bash)\" && \
                mamba activate ${ENV_DIR} && \
                INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${SEQ_FILE_LIST}) && \
                Rscript ./src/deg/run_deg_sequential_contrasts.R \
                    --input_file \"\$INPUT_FILE\" \
                    --output_dir \"${FULL_OUTPUT_DIR}\" \
                    --config \"${deg_config}\" \
                    ${ARG_F} ${ARG_S} ${ARG_Q}"
    
    echo "  ✓ Sequential DEG array job completed"
}

# ============================================================================
# Function: Count observations in h5ad file using JSON table
# ============================================================================

count_observations() {
    local H5AD_FILE=$1
    local JSON_FILE=$2
    
    # Get basename of the file for lookup
    local FILE_BASENAME=$(basename "$H5AD_FILE")
    
    # Check if JSON file exists
    if [ ! -f "$JSON_FILE" ]; then
        echo "Error: n_obs JSON file not found: $JSON_FILE" >&2
        return 1
    fi
    
    # Extract n_obs using basename as key
    local N_OBS=$(jq -r ".[\"${FILE_BASENAME}\"] // empty" "$JSON_FILE" 2>/dev/null)
    
    # Check if we got a valid number
    if [ -z "$N_OBS" ] || [ "$N_OBS" = "null" ]; then
        echo "Error: Could not find n_obs for file ${FILE_BASENAME} in JSON table" >&2
        return 1
    fi
    
    echo "$N_OBS"
}

# ============================================================================
# Main: Process each cell type file
# ============================================================================

echo ""
echo "============================================"
echo "DGE ANALYSIS"
echo "============================================"
echo "  Threshold for parallel processing: ${N_OBS_THRESHOLD} observations"
echo "  Processing ${N_FILES} file(s)..."
echo ""

# Determine JSON file path - look for .json file in the same directory as h5ad files
JSON_FILE=$(ls "$INPUT_DIR"/*.json 2>/dev/null | head -1)
if [ -n "$JSON_FILE" ]; then
    N_OBS_JSON_FILE="$JSON_FILE"
else
    # Fallback to n_obs_table.json if no .json file found
    N_OBS_JSON_FILE="${INPUT_DIR}/${N_OBS_JSON_FILE_NAME}"
fi

if [ "$MODE_J" = "True" ]; then
    # Arrays to store files for parallel and sequential processing
    PARALLEL_FILES=()
    SEQUENTIAL_FILES=()

    # First pass: classify all files
    for i in "${!INPUT_FILES[@]}"; do
        INPUT_FILE="${INPUT_FILES[$i]}"
        CELL_TYPE=$(basename "$INPUT_FILE" | sed 's/_processed\.h5ad$//' | sed 's/\.h5ad$//')
        
        echo "[$((i+1))/${N_FILES}] Classifying: $CELL_TYPE"
        
        # Count observations from JSON table
        echo "  Reading n_obs from JSON table..."
        N_OBS=$(count_observations "${INPUT_FILE}" "${N_OBS_JSON_FILE}")
        
        if [ -z "$N_OBS" ] || ! [[ "$N_OBS" =~ ^[0-9]+$ ]]; then
            echo "  Error: Failed to get n_obs for ${CELL_TYPE} from JSON table" >&2
            continue
        fi
        
        echo "  Found ${N_OBS} observations"
        
        # Classify based on threshold
        if [ "$N_OBS" -ge "$N_OBS_THRESHOLD" ]; then
            echo "  → PARALLEL BY CONTRASTS (n_obs=${N_OBS} >= ${N_OBS_THRESHOLD})"
            PARALLEL_FILES+=("${INPUT_FILE}")
        else
            echo "  → SEQUENTIAL BY CONTRASTS (n_obs=${N_OBS} < ${N_OBS_THRESHOLD})"
            SEQUENTIAL_FILES+=("${INPUT_FILE}")
        fi
    done

    # Process parallel DEG files with array jobs (steps 1-3 parallel, step 4 sequential per cell line)
    if [ ${#PARALLEL_FILES[@]} -gt 0 ]; then
        run_parallel_deg PARALLEL_FILES "${#PARALLEL_FILES[@]}"
    fi

    # Process sequential DEG files with array job
    if [ ${#SEQUENTIAL_FILES[@]} -gt 0 ]; then
        run_sequential_deg SEQUENTIAL_FILES "${#SEQUENTIAL_FILES[@]}"
    fi
else
    # MODE_J=False: Process all files with array job (no classification)
    run_sequential_deg INPUT_FILES "${N_FILES}"
fi

echo ""
echo "============================================"
echo "DGE PIPELINE COMPLETED"
echo "============================================"
echo "  Results saved to: ${FULL_OUTPUT_DIR}"
echo ""
echo "> DGE pipeline completed successfully!"
echo "> Cleaning up temporary files..."
rm -rf "${TMP_DIR_DEG}"
echo "> Cleanup complete"

