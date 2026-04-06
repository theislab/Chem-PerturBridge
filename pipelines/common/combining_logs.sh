#!/bin/bash

# Universal log combiner for OP3_v2 pipelines
# This script provides a convenient wrapper around the Python log aggregator

set -e

# Get project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${PROJECT_ROOT}"

# Activate mamba environment
ENV_DIR=./venv
export PATH="${HOME}/miniforge3/bin:${PATH}"
eval "$(mamba shell hook --shell bash)"
mamba activate ${ENV_DIR}

# Default values
LOG_DIR=""
PIPELINE_TYPE=""
OUTPUT_NAME=""
RECURSIVE=""

# Show help
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Aggregate pipeline log files into combined logs."
    echo ""
    echo "Options:"
    echo "  -r              Recursive mode: aggregate all logs"
    echo "  -l LOG_DIR      Log directory path"
    echo "  -t TYPE         Pipeline type: 'pseudobulk', 'deg', or 'enrich'"
    echo "  -o OUTPUT       Custom output filename"
    echo "  -h              Show this help"
    echo ""
    echo "Examples:"
    echo "  $0 -r"
    echo "  $0 -l logs/tahoe/full -t pseudobulk"
    echo "  $0 -l logs/sciplex/full/deg/.../filter_min_cells_10 -t deg -o my_log.log"
}

# Parse command line arguments
while getopts ":rl:t:o:h" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        r)
            RECURSIVE="true"
            ;;
        l)
            LOG_DIR="$OPTARG"
            ;;
        t)
            PIPELINE_TYPE="$OPTARG"
            ;;
        o)
            OUTPUT_NAME="$OPTARG"
            ;;
        \?)
            echo "Error: Invalid option: -$OPTARG" >&2
            echo "Use -h for help"
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument" >&2
            echo "Use -h for help"
            exit 1
            ;;
    esac
done

# Validate arguments
if [ "$RECURSIVE" = "true" ]; then
    # Recursive mode
    if [ -n "$PIPELINE_TYPE" ] || [ -n "$OUTPUT_NAME" ]; then
        echo "Error: Recursive mode (-r) cannot be combined with -t or -o" >&2
        echo "Recursive mode auto-detects pipeline types and creates multiple files"
        echo "Use -h for help"
        exit 1
    fi
    
    echo "Running recursive log aggregation..."
    
    # Build command
    CMD=(python3 -m src.utils.aggregate_logs --recursive)
    
    if [ -n "$LOG_DIR" ]; then
        CMD+=(--log_dir "$LOG_DIR")
    fi
    
    "${CMD[@]}"
    
else
    # Targeted mode
    if [ -z "$LOG_DIR" ]; then
        echo "Error: -l (log directory) is required for targeted mode" >&2
        echo "Use -h for help"
        exit 1
    fi
    
    if [ -z "$PIPELINE_TYPE" ]; then
        echo "Error: -t (pipeline type) is required for targeted mode" >&2
        echo "Use -h for help"
        exit 1
    fi
    
    # Validate pipeline type
    if [ "$PIPELINE_TYPE" != "pseudobulk" ] && [ "$PIPELINE_TYPE" != "deg" ] && [ "$PIPELINE_TYPE" != "enrich" ]; then
        echo "Error: Pipeline type must be 'pseudobulk', 'deg', or 'enrich'" >&2
        exit 1
    fi
    
    # Build command
    CMD=(python3 -m src.utils.aggregate_logs --log_dir "$LOG_DIR" --pipeline_type "$PIPELINE_TYPE")
    
    if [ -n "$OUTPUT_NAME" ]; then
        CMD+=(--output_name "$OUTPUT_NAME")
    fi
    
    echo "Running targeted log aggregation..."
    echo "  Log directory: $LOG_DIR"
    echo "  Pipeline type: $PIPELINE_TYPE"
    if [ -n "$OUTPUT_NAME" ]; then
        echo "  Output name: $OUTPUT_NAME"
    fi
    echo ""
    
    "${CMD[@]}"
fi

echo ""
echo "Log aggregation complete!"
