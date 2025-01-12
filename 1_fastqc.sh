#!/bin/bash
set -e
set -u
set -o pipefail

# Default parameters
THREADS=8
INPUT_DIR=""
OUTPUT_DIR=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 -i <input_dir> -o <output_dir> [-p threads]"
    exit 1
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}/fastqc"
mkdir -p "${OUTPUT_DIR}/logs"

# Run FastQC
echo "=== Running FastQC on raw reads ==="
for fastq in "$INPUT_DIR"/*_L001_R[12]_001.fastq.gz; do
    echo "Processing $(basename "$fastq")..."
    fastqc \
        --outdir "${OUTPUT_DIR}/fastqc" \
        --threads "$THREADS" \
        "$fastq"
done