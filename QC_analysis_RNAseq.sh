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