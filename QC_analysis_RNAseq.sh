#!/bin/bash

# RNA-seq Quality Control Pipeline

set -e # Exit on error
set -u # Exit on undefined variable
set -o pipefail # Exit on pipe failure

# Parameters
THREADS=8
ADAPTER="TruSeq"
MIN_LENGTH=50
MIN_QUALITY=20

# Function to print usage
usage() {
    echo "Usage: $0 <input_dir> <output_dir> <genome_index>"
    echo "Options:"
    echo "  -t|--threads INT        Number of threads (default: 8)"
    echo "  -a|--adapter STR       Adapter type (default: TruSeq)"
    echo "  -l|--min-length INT    Minimum read length (default: 50)"
    echo "  -q|--min-quality INT   Minimum base quality (default: 20)"
    echo "  -h|--help             Print this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -a|--adapter)
            ADAPTER="$2"
            shift 2
            ;;
        -l|--min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        -q|--min-quality)
            MIN_QUALITY="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            INPUT_DIR="$1"
            OUTPUT_DIR="$2"
            GENOME_INDEX="$3"
            shift 3
            ;;
    esac
done

if [[ -z "${INPUT_DIR:-}" ]] || [[ -z "${OUTPUT_DIR:-}" ]] || [[ -z "${GENOME_INDEX:-}" ]]; then
    echo "Error: Missing required arguments"
    usage
fi
