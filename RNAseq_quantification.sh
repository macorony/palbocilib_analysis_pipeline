#!/bin/bash

# RNA-seq Analysis Pipeline using Salmon
# This script performs RNA-seq analysis including QC, trimming, and Salmon quantification

# Set up error handling
set -e
set -u
set -o pipefail

# Default parameters
THREADS=8
ADAPTER="TruSeq"
MIN_LENGTH=50
MIN_QUALITY=20
LIBTYPE="A"  # Automatically detect library type

# Configuration directories
INPUT_DIR=""
OUTPUT_DIR=""
SALMON_INDEX=""
GTF_FILE=""
TRANSCRIPTOME=""

# Create help message
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -x <salmon_index> -t <transcriptome> [options]"
    echo
    echo "Required arguments:"
    echo "  -i, --input        Input directory containing raw fastq files"
    echo "  -o, --output       Output directory for results"
    echo "  -x, --index        Salmon index directory"
    echo "  -t, --transcriptome Reference transcriptome (FASTA)"
    echo
    echo "Optional arguments:"
    echo "  -p, --threads      Number of threads (default: 8)"
    echo "  --adapter          Adapter type for trimming (default: TruSeq)"
    echo "  --min-length       Minimum read length after trimming (default: 50)"
    echo "  --min-quality      Minimum base quality for trimming (default: 20)"
    echo "  --libtype          Library type (default: A - auto-detect)"
    echo "  -h, --help         Show this help message"
    exit 1
}

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
        -t|--transcriptome)
            TRANSCRIPTOME="$2"
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
        --libtype)
            LIBTYPE="$2"
            shift 2
            ;;
        --min-length)
            MIN_LENGTH="$2"
            shift 2
            ;;
        --min-quality)
            MIN_QUALITY="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$SALMON_INDEX" ]] || [[ -z "$TRANSCRIPTOME" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Create output directory structure
mkdir -p "${OUTPUT_DIR}"/{fastqc,trimmed,salmon,logs,multiqc}
LOG_DIR="${OUTPUT_DIR}/logs"

# Start logging
exec 1> >(tee "${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log")
exec 2>&1

# Function to check dependencies
check_dependencies() {
    local required_tools=("fastqc" "trimmomatic" "salmon" "multiqc")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        echo "ERROR: Missing required tools: ${missing_tools[*]}"
        exit 1
    fi
}

# Function to build Salmon index if it doesn't exist
build_salmon_index() {
    if [ ! -d "$SALMON_INDEX" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building Salmon index..."
        salmon index \
            -t "$TRANSCRIPTOME" \
            -i "$SALMON_INDEX" \
            -p "$THREADS"
    fi
}

# Function to run FastQC
run_fastqc() {
    local input_file="$1"
    local output_dir="$2"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FastQC on $(basename "$input_file")..."
    fastqc \
        --outdir "$output_dir" \
        --threads "$THREADS" \
        "$input_file"
}

# Function to run Trimmomatic
run_trimmomatic() {
    local r1="$1"
    local r2="${r1/_R1/_R2}"
    local basename=$(basename "$r1" _R1.fastq.gz)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimming $basename..."
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
}

# Function to run Salmon quantification
run_salmon() {
    local r1="$1"
    local r2="${r1/_R1/_R2}"
    local basename=$(basename "$r1" _R1_paired.fastq.gz)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Salmon quantification on $basename..."
    salmon quant \
        -i "$SALMON_INDEX" \
        -l "$LIBTYPE" \
        -1 "$r1" \
        -2 "$r2" \
        -p "$THREADS" \
        --validateMappings \
        -o "${OUTPUT_DIR}/salmon/${basename}"
}

# Main pipeline execution
main() {
    echo "=== RNA-seq Pipeline Started at $(date) ==="
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    echo "Salmon index: $SALMON_INDEX"
    echo "Threads: $THREADS"
    
    # Check dependencies
    check_dependencies
    
    # Build Salmon index if needed
    build_salmon_index
    
    # Run FastQC on raw reads
    echo "=== Running FastQC on raw reads ==="
    find "$INPUT_DIR" -name "*_R[12].fastq.gz" | \
        parallel -j "$THREADS" run_fastqc {} "${OUTPUT_DIR}/fastqc"
    
    # Run Trimmomatic
    echo "=== Running Trimmomatic ==="
    find "$INPUT_DIR" -name "*_R1.fastq.gz" | \
        parallel -j 1 run_trimmomatic {}
    
    # Run FastQC on trimmed reads
    echo "=== Running FastQC on trimmed reads ==="
    find "${OUTPUT_DIR}/trimmed" -name "*_paired.fastq.gz" | \
        parallel -j "$THREADS" run_fastqc {} "${OUTPUT_DIR}/fastqc"
    
    # Run Salmon quantification
    echo "=== Running Salmon Quantification ==="
    find "${OUTPUT_DIR}/trimmed" -name "*_R1_paired.fastq.gz" | \
        parallel -j 1 run_salmon {}
    
    # Run MultiQC
    echo "=== Running MultiQC ==="
    multiqc \
        --force \
        --outdir "${OUTPUT_DIR}/multiqc" \
        "$OUTPUT_DIR"
    
    # Generate summary report
    echo "=== Generating Summary Report ==="
    {
        echo "Sample,Raw Reads,Trimmed Reads,Mapped Reads,Mapping Rate"
        for quant_file in "${OUTPUT_DIR}"/salmon/*/quant.sf; do
            sample=$(basename $(dirname "$quant_file"))
            raw_reads=$(zcat "${INPUT_DIR}/${sample}_R1.fastq.gz" | echo $((`wc -l`/4)))
            trimmed_reads=$(zcat "${OUTPUT_DIR}/trimmed/${sample}_R1_paired.fastq.gz" | echo $((`wc -l`/4)))
            mapped_reads=$(grep "NumMappedFragments" "$(dirname "$quant_file")/aux_info/meta_info.json" | cut -d'"' -f4)
            mapping_rate=$(grep "MappingRate" "$(dirname "$quant_file")/aux_info/meta_info.json" | cut -d'"' -f4)
            echo "$sample,$raw_reads,$trimmed_reads,$mapped_reads,$mapping_rate"
        done
    } > "${OUTPUT_DIR}/pipeline_summary.csv"
    
    echo "=== Pipeline completed at $(date) ==="
    echo "Summary report saved to ${OUTPUT_DIR}/pipeline_summary.csv"
}

# Execute main pipeline
main