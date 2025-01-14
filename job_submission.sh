#!/bin/bash
#SBATCH --account=def-lebrun    
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=gang.yan@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=RNA_seq_quantification
#SBATCH --output=RNA_seq_%j.out
#SBATCH --error=RNA_seq_%j.err

# Load required modules
module load StdEnv/2020
module load fastqc/0.11.9
module load trimmomatic/0.39
module load salmon/1.4.0
module load python/3.8.10
module load multiqc/1.9

# Create project directory structure in scratch
PROJECT_DIR=/home/gyan/projects/def-lebrun/gyan/RNAseq_CDK4_6

# Set working directory
RNA_SEQ_SCRIPT= $PROJECT_DIR/RNAseq_quantification.sh

INPUT_DIR=$PROJECT_DIR/raw_data
OUTPUT_DIR=$PROJECT_DIR/output
SALMON_INDEX=$PROJECT_DIR/salmon_index
TRANSCRIPTOME=$PROJECT_DIR/transcriptome"



CMD="$RNA_SEQ_SCRIPT \
  -i '$INPUT__DIR' \
  -o '$OUTPUT_DIR' \
  -x '$SALMON_INDEX' \
  -t '$TRANSCRIPTOME' \
  -p '$SLURM_CPUS_PER_TASK' \
  --libtype 'IU'" 

eval "$CMD"
