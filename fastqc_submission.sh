#!/bin/bash
#SBATCH --account=def-username     # Replace with your account/allocation
#SBATCH --time=4:00:00            # Adjust time based on your data size
#SBATCH --mem=16G                 # Memory per node
#SBATCH --cpus-per-task=8         # Match with THREADS in your script
#SBATCH --job-name=fastqc         # Job name
#SBATCH --output=%x-%j.out        # Standard output log
#SBATCH --error=%x-%j.err         # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL    # Email notification
#SBATCH --mail-user=your.email@institution.ca  # Replace with your email

# Load required modules
module load StdEnv/2020
module load fastqc/0.11.9

# Navigate to your project directory
cd $SLURM_SUBMIT_DIR

# Run the FastQC script
./1_fastqc.sh -i /path/to/input/directory \
              -o /path/to/output/directory \
              -p 8