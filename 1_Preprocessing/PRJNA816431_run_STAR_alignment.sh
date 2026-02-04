#!/bin/bash
#SBATCH --job-name=star_alignment
#SBATCH --output=star_alignment_%A_%a.out
#SBATCH --error=star_alignment_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-21  # Adjust this to the number of samples

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


# Get full sample names (without .fastq)
SAMPLES=($(ls $INPUT_DIR/*.fastq | xargs -n 1 basename | sed 's/.fastq//'))

# Select sample based on SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Build full path to FASTQ
FASTQ="${INPUT_DIR}/${SAMPLE}.fastq"

# Final check: file must exist
if [ ! -f "$FASTQ" ]; then
  echo "FASTQ file not found: $FASTQ"
  exit 1
fi

# Set output directory
OUTPUT_DIR="/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/new_data/PRJNA816431_callus/STAR_alignment_UCSC_genome/${SAMPLE}"

mkdir -p $OUTPUT_DIR

# Run STAR alignment
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_DIR \
     --readFilesIn $FASTQ \
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
