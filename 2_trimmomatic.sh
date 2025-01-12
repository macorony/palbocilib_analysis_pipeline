#!/bin/bash
set -e
set -u
set -o pipefail

# Default parameters
THREADS=8
ADAPTER="TruSeq"
MIN_LENGTH=50
MIN_QUALITY=20
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
        --adapter)
            ADAPTER="$2"
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
    echo "Usage: $0 -i <input_dir> -o <output_dir> [-p threads] [--adapter adapter_type]"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/trimmed"

# Run Trimmomatic
for r1 in "$INPUT_DIR"/*_R1.fastq.gz; do
    r2="${r1/_R1/_R2}"
    basename=$(basename "$r1" _R1.fastq.gz)
    
    echo "Trimming $basename..."
    trimmomatic PE \
        -threads "$THREADS" \
        "$r1" "$r2" \
        "${OUTPUT_DIR}/trimmed/${basename}_R1_paired.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/${basename}_R1_unpaired.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/${basename}_R2_paired.fastq.gz" \
        "${OUTPUT_DIR}/trimmed/${basename}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"${ADAPTER}:2:30:10" \
        LEADING:"$MIN_QUALITY" \
        TRAILING:"$MIN_QUALITY" \
        SLIDINGWINDOW:4:"$MIN_QUALITY" \
        MINLEN:"$MIN_LENGTH"
done