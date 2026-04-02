#!/bin/bash

set -e

# ============================================================================
# DEG (Differential Expression) Pipeline
# ============================================================================
#
# This pipeline runs differential expression analysis on pseudobulk data:
#   STEP 1: Preprocess pseudobulk (add labels, split by cell type if -j)
#   STEP 2: Run DEG analysis (R-based, per cell type in parallel or sequential)
#   STEP 3: Aggregate batch results (if parallel mode was used)
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Pipeline settings
PIPELINE_NAME="deg"
ENV_DIR=./venv
LOGS_DIR=./logs
MAX_CONCURRENT=50
# QOS/PARTITION: set by -g flag (gpu_normal/gpu_p) or default (cpu_normal/cpu_p)
QOS=cpu_normal
PARTITION=cpu_p

# Default values
MODE_S=False
MODE_J=False
MODE_Q=False
MODE_N=False
MODE_G=False
PAR=""
FILT=""
CONFIG=""
DATASET=""
ARG_S=""
ARG_F=""
ARG_N=""


# Valid parameters
DEG_PARAMETERS=(
    "group_all_replicates"
    "separate_replicates"
)

while getopts ":sjqp:f:c:d:nhg" opt; do
        case $opt in
                h)
                        echo "Run: $0 [-s] [-j] [-h] [-g] [-f] [-q] [-n] -p (parameters: ${DEG_PARAMETERS[*]}) -c CONFIG -d DATASET"
                        echo ""
                        echo "Required arguments:"
                        echo "  -d  Dataset name (sciplex, tahoe, l1000_phase1, l1000_phase2)"
                        echo "  -p  Design parameter (group_all_replicates, separate_replicates)"
                        echo "  -c  Path to config file (auto-set by run_deg.sh)"
                        echo ""
                        echo "Optional flags:"
                        echo "  -s  Subsample mode - limits DEG analysis for testing/debugging"
                        echo "      NOTE: Does NOT change input files (those are defined in config)"
                        echo "  -j  Parallel mode - run array jobs per cell type (faster)"
                        echo "  -g  Use GPU QOS/partition (gpu_normal/gpu_p); default is CPU (cpu_normal/cpu_p)"
                        echo "  -f  VALUE - Filter: minimum cells per sample (e.g., -f 50)"
                        echo "  -q  Quality control filter - remove samples that failed QC"
                        echo "  -n  Normalized dataset: pass --normalized to preprocessing (Ensembl: drop ambiguous"
                        echo "      targets, no count aggregation) and to R DEG (skip voom/normalization)"
                        echo "  -h  Show this help message"
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
                n)
                        MODE_N=True
                        ;;
                g)
                        MODE_G=True
                        ;;
        esac
done

if [ "$MODE_G" = True ]; then
    QOS=gpu_normal
    PARTITION=gpu_p
fi

# ============================================================================
# Validate arguments
# ============================================================================

if [[ " ${DEG_PARAMETERS[*]} " != *" $PAR "* ]]; then
    echo "Error: -p must be set up and one of: ${DEG_PARAMETERS[*]}" >&2
    exit 1
fi

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
if [ "$MODE_Q" = "True" ]; then
    echo "> Filter samples that did not pass quality control"
    QC_FOLDER="qc_true"
    ARG_Q="--qc"
else
    echo "> Do not filter samples that did not pass quality control"
    QC_FOLDER="qc_false"
    ARG_Q=""
fi

# Set subsample flag if flag is provided
if [ "$MODE_S" = "True" ]; then
    echo "> Work with a subsample"
    SUBDIR="subsample"
    ARG_S="--subsampling"
else
	echo "> Work with a full version"
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

# Set normalized flag argument
if [ "$MODE_N" = "True" ]; then
        echo "> Dataset is normalized, skipping normalization steps"
        ARG_N="--normalized"
else
        echo "> Dataset is not normalized, running normalization steps"
        ARG_N=""
fi

# ============================================================================
# Initialize environment
# ============================================================================

# Add miniforge3/bin to PATH so mamba can be found
export PATH="$HOME/miniforge3/bin:$PATH"

eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

# Setup temporary directories
TMP_DIR="${PROJECT_ROOT}/tmp"
TMP_DIR_DEG="${TMP_DIR}/deg_$$"
mkdir -p "${TMP_DIR}" "${TMP_DIR_DEG}"
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
# Helper Functions
# ============================================================================

extract_cell_type() {
    local INPUT_FILE=$1
    # If file is in a /by_celltype/ directory, extract cell type from directory name
    # Otherwise, return empty string (file is not split by cell type)
    if [[ "$INPUT_FILE" == *"/by_celltype/"* ]]; then
        basename "$(dirname "$INPUT_FILE")"
    else
        # For files not in /by_celltype/, return empty string
        # This indicates the file contains multiple cell types or is not split by cell type
        echo ""
    fi
}

# Common sbatch preamble
SBATCH_PREAMBLE="export PATH=\${HOME}/miniforge3/bin:\${PATH} && export TMPDIR=\${HOME}/tmp && mkdir -p \${TMPDIR} && cd ${PROJECT_ROOT} && eval \"\$(mamba shell hook --shell bash)\" && mamba activate ${ENV_DIR}"

# ============================================================================
# STEP 1: PREPROCESS PSEUDOBULK
# ============================================================================

echo "> Preprocess pseudobulk with parameter: $PAR"

if ! par_process=$(jq -e ".$PAR.par_process" "$CONFIG"); then
    echo "Error: Failed to extract par_process from $CONFIG" >&2
    exit 1
fi

preprocess_config="${TMP_DIR_DEG}/preprocess_config.json"
echo "$par_process" | jq "." > "$preprocess_config"



echo "> Running preprocessing pseudobulk"
sbatch -W -J deg_processing_pseudobulk \
    --partition=${PARTITION} \
    --qos=${QOS} \
    --mem=500G \
    --time=2:00:00 \
    --cpus-per-task=2 \
    --exclude=supercpu01,gpusrv36,gpusrv56,gpusrv60,gpusrv45 \
    --output="${LOGS_BASE_DIR}/preprocessing/deg_processing_pseudobulk.%j.out" \
    --error="${LOGS_BASE_DIR}/preprocessing/deg_processing_pseudobulk.%j.err" \
    --wrap="${SBATCH_PREAMBLE} && \
            python3 -m src.deg.run_processing_pseudobulk ${ARG_J} ${ARG_N} \
            --config ${preprocess_config}"

echo "> Preprocessing completed"

# ============================================================================
# STEP 2: DGE ANALYSIS
# ============================================================================

echo "> Starting DGE analysis with parameter: $PAR"

if ! par_deg=$(jq -e ".$PAR.par_deg" "$CONFIG"); then
    echo "Error: Failed to extract par_deg from $CONFIG" >&2
    exit 1
fi

# Extract paths from config
INPUT_DIR=$(echo "$par_deg" | jq -r '.input_dir')
OUTPUT_DIR=$(echo "$par_deg" | jq -r '.output_dir')

# Determine search directory and memory based on mode
if [ "$MODE_J" = "True" ]; then
    SEARCH_DIR="${INPUT_DIR}/by_celltype"
    MEM=50G
else
    SEARCH_DIR="${INPUT_DIR}"
    MEM=250G
fi

# Validate search directory
if [ ! -d "$SEARCH_DIR" ]; then
    echo "Error: Directory not found: $SEARCH_DIR" >&2
    exit 1
fi

# Find input files
if [ "$MODE_J" = "True" ]; then
    mapfile -t INPUT_FILES < <(find "$SEARCH_DIR" -mindepth 2 -maxdepth 2 -name "*.h5ad" -type f 2>/dev/null | sort)
else
    mapfile -t INPUT_FILES < <(find "$SEARCH_DIR" -maxdepth 1 -name "*.h5ad" -type f 2>/dev/null | sort)
fi

N_FILES=${#INPUT_FILES[@]}

if [ "$N_FILES" -eq 0 ]; then
    echo "Error: No .h5ad files found in $SEARCH_DIR" >&2
    exit 1
fi

echo "  Found $N_FILES file(s) to process"

# Create DEG config
deg_config="${TMP_DIR_DEG}/deg_config.json"
echo "$par_deg" | jq "." > "$deg_config"

# Setup output directories
BASE_OUTPUT_DIR="${OUTPUT_DIR}/${SUBDIR}/${QC_FOLDER}/${FILTER_FOLDER}"
RESULTS_DIR="${BASE_OUTPUT_DIR}/results"
INTERMEDIATE_DIR="${BASE_OUTPUT_DIR}/intermediate"
mkdir -p "${LOGS_BASE_DIR}/deg"

# ============================================================================
# Function: Run DEG
# ============================================================================

run_deg() {
    local FILES_ARRAY_NAME=$1
    local N_FILES=$2
    
    local -n FILES_ARRAY=$FILES_ARRAY_NAME
    
    echo ""
    echo "============================================"
    echo "DEG ANALYSIS: Array Job Mode"
    echo "============================================"
    echo "  Processing ${N_FILES} file(s) with array job..."
    
    # Extract unique cell types and create directories
    declare -A CELL_TYPES
    for INPUT_FILE in "${FILES_ARRAY[@]}"; do
        CELL_TYPE=$(extract_cell_type "$INPUT_FILE")
        # Only add non-empty cell types to the array (empty means file is not in by_celltype/)
        if [[ -n "$CELL_TYPE" ]]; then
            CELL_TYPES["$CELL_TYPE"]=1
        fi
    done
    
    # Create results directory for combined files (no cell type)
    mkdir -p "${RESULTS_DIR}"
    
    # Create base log directory (needed for files without cell types)
    mkdir -p "${LOGS_BASE_DIR}/deg"
    
    # Create intermediate directories for separate files (with cell types)
    # Logs go directly to cell_type folders, no cleanup needed
    for CELL_TYPE in "${!CELL_TYPES[@]}"; do
        mkdir -p "${INTERMEDIATE_DIR}/${CELL_TYPE}" "${LOGS_BASE_DIR}/deg/${CELL_TYPE}"
    done
    
    # Create file list and cell types list
    DEG_FILE_LIST="${TMP_DIR_DEG}/deg_file_list.txt"
    DEG_CELL_TYPES_LIST="${TMP_DIR_DEG}/deg_cell_types_list.txt"
    printf '%s\n' "${FILES_ARRAY[@]}" > "$DEG_FILE_LIST"
    for INPUT_FILE in "${FILES_ARRAY[@]}"; do
        extract_cell_type "$INPUT_FILE"
    done > "$DEG_CELL_TYPES_LIST"
    
    echo "  Submitting array job [${N_FILES} jobs, max ${MAX_CONCURRENT} concurrent]..."
    sbatch -W \
        -J deg_analysis \
        --array=0-$((N_FILES-1))%${MAX_CONCURRENT} \
        --partition=${PARTITION} \
        --qos=${QOS} \
        --mem=${MEM} \
        --time=10:00:00 \
        --cpus-per-task=2 \
        --exclude=supercpu01,gpusrv36,gpusrv56,gpusrv60,gpusrv45 \
        --output=/dev/null \
        --error=/dev/null \
        --wrap="${SBATCH_PREAMBLE} && \
                INPUT_FILE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${DEG_FILE_LIST}) && \
                CELL_TYPE=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID + 1))p\" ${DEG_CELL_TYPES_LIST}) && \
                if [[ -z \"\${CELL_TYPE}\" ]]; then \
                    LOG_DIR=\"${LOGS_BASE_DIR}/deg\"; \
                    OUTPUT_DIR=\"${RESULTS_DIR}\"; \
                else \
                    LOG_DIR=\"${LOGS_BASE_DIR}/deg/\${CELL_TYPE}\"; \
                    OUTPUT_DIR=\"${INTERMEDIATE_DIR}/\${CELL_TYPE}\"; \
                fi && \
                mkdir -p \"\${LOG_DIR}\" \"\${OUTPUT_DIR}\" && \
                exec > \"\${LOG_DIR}/deg_analysis.\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out\" 2> \"\${LOG_DIR}/deg_analysis.\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.err\" && \
                Rscript ./src/deg/run_deg.R \
                    --input_file \"\$INPUT_FILE\" \
                    --output_dir \"\${OUTPUT_DIR}\" \
                    --config \"${deg_config}\" \
                    ${ARG_F} ${ARG_S} ${ARG_Q} ${ARG_N}"
    
    echo "> DEG array job completed"
}

# ============================================================================
# Main execution
# ============================================================================

echo ""
echo "============================================"
echo "DGE ANALYSIS"
echo "============================================"
echo "  Processing ${N_FILES} file(s)..."
echo ""

run_deg INPUT_FILES "${N_FILES}"

# ============================================================================
# STEP 3: AGGREGATE BATCH FILES
# ============================================================================

echo ""
echo "============================================"
echo "AGGREGATING BATCH FILES"
echo "============================================"

# Check if intermediate directory has files and results directory is empty
if [ -d "$INTERMEDIATE_DIR" ] && [ "$(ls -A $INTERMEDIATE_DIR 2>/dev/null)" ]; then
    # Check if results directory is empty
    if [ ! -d "$RESULTS_DIR" ] || [ -z "$(ls -A $RESULTS_DIR 2>/dev/null)" ]; then
        echo "  Found files in ${INTERMEDIATE_DIR} and results directory is empty"
        echo "  Aggregating batch files from intermediate/ to results/..."
        
        mkdir -p "${LOGS_BASE_DIR}/aggregation"
        
        echo "  Running aggregation..."
        sbatch -W -J deg_aggregation \
            --partition=${PARTITION} \
            --qos=${QOS} \
            --mem=500G \
            --time=2:00:00 \
            --cpus-per-task=2 \
            --exclude=supercpu01,gpusrv36,gpusrv56,gpusrv60,gpusrv45 \
            --output="${LOGS_BASE_DIR}/aggregation/deg_aggregation.%j.out" \
            --error="${LOGS_BASE_DIR}/aggregation/deg_aggregation.%j.err" \
            --wrap="${SBATCH_PREAMBLE} && \
                    python3 -m src.deg.aggregating_deg \
                    --input_dir ${INTERMEDIATE_DIR} \
                    --output_dir ${RESULTS_DIR}"
        
        echo "> Aggregation completed"
    else
        echo "  Results directory is not empty, skipping aggregation"
    fi
else
    echo "  No intermediate directory found or it is empty, skipping aggregation"
fi

echo ""
echo "============================================"
echo "DGE PIPELINE COMPLETED"
echo "============================================"
echo "  Results saved to: ${BASE_OUTPUT_DIR}"
echo "    - intermediate/ for separate files"
echo "    - results/ for combined"
echo ""
echo "> DGE pipeline is finished"
echo "> Cleaning up temporary files..."
rm -rf "${TMP_DIR_DEG}"
echo "> Cleanup complete"
