#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --output=featurecounts_%A_%a.out
#SBATCH --error=featurecounts_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-9

# Paths 
GTF="/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/reference_genomes/mm39/mm39_refGene.gtf"
OUTDIR="/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/new_data/PRJNA1173774_BMDM/UCSC_featurecounts_out"
THREADS=8

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Get list of BAM files
BAM_FILES=(/scratch1/eschwab/Benayoun_Lee_Lab/ATOM/new_data/PRJNA1173774_BMDM/STAR_alignment_UCSC_genome/*/*Aligned.sortedByCoord.out.bam)

# Get BAM file for this array task
BAM=${BAM_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Extract sample name
SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing sample: $SAMPLE"
echo "BAM file: $BAM"

featureCounts \
    -T $THREADS \
    -t exon \
    -g gene_name \
    -p \
    --primary \
    -Q 20 \
    -a "$GTF" \
    -o "$OUTDIR/${SAMPLE}_featureCounts.txt" \
    "$BAM"

echo "Done processing $SAMPLE"

