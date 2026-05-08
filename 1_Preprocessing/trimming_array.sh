#!/bin/bash
#SBATCH --job-name=trimming_array
#SBATCH --output=trimming_%A_%a.out
#SBATCH --error=trimming_%A_%a.err
#SBATCH --array=1-12             
#SBATCH --time=12:00:00          
#SBATCH --mem=8G                 
#SBATCH --cpus-per-task=4        
#SBATCH --partition=main 


# Get input and output directories from arguments
INPUT_LIST=$1
OUTPUT_DIR=$2

TRIMMED_DIR="$OUTPUT_DIR/adapter_trimmed"  # Intermediate output
FINAL_DIR="$OUTPUT_DIR/hard_trimmed"  # Final output

mkdir -p $OUTPUT_DIR $TRIMMED_DIR $FINAL_DIR

# Get the correct line (1-based index for SLURM)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")

# Parse out the two FASTQ file paths
R1_FILE=$(echo $LINE | cut -d ' ' -f 1)
R2_FILE=$(echo $LINE | cut -d ' ' -f 2)

# Extract a sample name (remove directory and _1.fastq/_2.fastq)
BASENAME=$(basename "$R1_FILE")
SAMPLE=${BASENAME%_1.fastq}


# Step 1: Adapter removal and quality trimming with Trim Galore
trim_galore \
    --paired \
    --fastqc \
    --cores 4 \
    --illumina \
    -o "$TRIMMED_DIR" \
    "$R1_FILE" "$R2_FILE"

# Define the output files from Trim Galore
TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE}_1_val_1.fq"
TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE}_2_val_2.fq"

# Step 2: Hard trim the reads to 75 bp using FASTX Toolkit
fastx_trimmer -i "$TRIMMED_R1" -o "${FINAL_DIR}/${SAMPLE}_1_trimmed.fastq" -l 75 -f 7
fastx_trimmer -i "$TRIMMED_R2" -o "${FINAL_DIR}/${SAMPLE}_2_trimmed.fastq" -l 75 -f 7
