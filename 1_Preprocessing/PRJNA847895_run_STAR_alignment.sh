#!/bin/bash
#SBATCH --job-name=star_alignment
#SBATCH --output=star_alignment_%A_%a.out
#SBATCH --error=star_alignment_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-12  

# Set the number of threads
THREADS=8

# Read input FASTQ directory from argument
INPUT_DIR=$1

if [ -z "$INPUT_DIR" ]; then
  echo "Usage: sbatch --array=1-N star_alignment.sh /path/to/fastqs"
  exit 1
fi

GENOME_PATH='/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/reference_genomes/mm39'
GENOME_DIR="${GENOME_PATH}/STAR_genome_index"

# Get sample list
SAMPLES=($(ls $INPUT_DIR/*_1.fastq | xargs -n 1 basename | sed 's/_1.fastq//'))

# Select sample for this SLURM array task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define file paths
R1_FASTQ="${INPUT_DIR}/${SAMPLE}_1.fastq"
R2_FASTQ="${INPUT_DIR}/${SAMPLE}_2.fastq"
OUTPUT_DIR="/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/new_data/PRJNA847895_BMDM_multilane/STAR_alignment_UCSC_genome/${SAMPLE}"

mkdir -p $OUTPUT_DIR

# Run STAR alignment with specified defaults
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_DIR \
     --readFilesIn $R1_FASTQ $R2_FASTQ \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 0.06 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFileNamePrefix $OUTPUT_DIR/${SAMPLE}_ \
     --outSAMtype BAM SortedByCoordinate

