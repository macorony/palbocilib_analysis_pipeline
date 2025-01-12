#!/bin/bash
set -e
set -u
set -o pipefail

# Default parameters
THREADS=8
LIBTYPE="A"
INPUT_DIR=""
OUTPUT_DIR=""
SALMON_INDEX=""

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
        -x|--index)
            SALMON_INDEX="$2"
            shift 2
            ;;
        -p|--threads)
            THREADS="$2"
            shift 2
            ;;
        --libtype)
            LIBTYPE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$SALMON_INDEX" ]]; then
    echo "Usage: $0 -i <input_dir> -o <output_dir> -x <salmon_index> [-p threads] [--libtype type]"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/salmon"

# Run Salmon quantification
for r1 in "${INPUT_DIR}"/trimmed/*_R1_paired.fastq.gz; do
    r2="${r1/_R1/_R2}"
    basename=$(basename "$r1" _R1_paired.fastq.gz)
    
    echo "Quantifying $basename..."
    salmon quant \
        -i "$SALMON_INDEX" \
        -l "$LIBTYPE" \
        -1 "$r1" \
        -2 "$r2" \
        -p "$THREADS" \
        --validateMappings \
        -o "${OUTPUT_DIR}/salmon/${basename}"
done