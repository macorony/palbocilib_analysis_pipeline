#!/bin/bash
#SBATCH --account=def-lebrun    
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=124G
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

# Set working directory
cd $SLURM_SUBMIT_DIR

# Create project directory structure in scratch
PROJECT_DIR=$SCRATCH/RNA_seq_project
mkdir -p $PROJECT_DIR

# Run the pipeline
./RNA_seq_pipeline.sh \
    -i $PROJECT_DIR/raw_data \
    -o $PROJECT_DIR/results \
    -x $PROJECT_DIR/salmon_index \
    -t $PROJECT_DIR/reference/transcriptome.fa \
    -p $SLURM_CPUS_PER_TASK 