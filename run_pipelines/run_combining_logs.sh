#!/bin/bash

# Wrapper script for running log aggregation
# Validates arguments and calls the common combining_logs.sh script

set -e

# Get the script directory and change to project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# Valid pipeline type choices
VALID_TYPES=(
    "pseudobulk"
    "deg"
)

# Show help
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Aggregate pipeline log files into combined logs."
    echo ""
    echo "Options:"
    echo "  -r              Recursive mode: aggregate all logs under ./logs"
    echo "  -l LOG_DIR      Log directory path (relative to project root)"
    echo "  -t TYPE         Pipeline type: 'pseudobulk' or 'deg'"
    echo "  -o OUTPUT       Custom output filename (default: combined.log)"
    echo "  -h              Show this help"
    echo ""
    echo "Examples:"
    echo "  $0 -r                                                # Aggregate all logs"
    echo "  $0 -l logs/tahoe/full -t pseudobulk                 # Specific dataset"
    echo "  $0 -l logs/sciplex/full/deg/.../filter_min_cells_10 -t deg  # Specific DEG config"
    echo ""
}

# Initialize variables
RECURSIVE=""
LOG_DIR=""
PIPELINE_TYPE=""
OUTPUT_NAME=""

# Parse command line arguments
while getopts ":rl:t:o:h" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        r)
            RECURSIVE="-r"
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

# Validate pipeline type if provided
if [ -n "$PIPELINE_TYPE" ]; then
    if [[ " ${VALID_TYPES[*]} " != *" $PIPELINE_TYPE "* ]]; then
        echo "Error: Pipeline type must be one of: ${VALID_TYPES[*]}" >&2
        exit 1
    fi
fi

# Build command for common script
CMD=(./pipelines/common/combining_logs.sh)

if [ -n "$RECURSIVE" ]; then
    CMD+=("$RECURSIVE")
fi

if [ -n "$LOG_DIR" ]; then
    CMD+=(-l "$LOG_DIR")
fi

if [ -n "$PIPELINE_TYPE" ]; then
    CMD+=(-t "$PIPELINE_TYPE")
fi

if [ -n "$OUTPUT_NAME" ]; then
    CMD+=(-o "$OUTPUT_NAME")
fi

# Execute the command
"${CMD[@]}"
