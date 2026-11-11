#!/bin/bash
#SBATCH --array=0-25
#SBATCH --partition=Standard.2.0
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G



N=${SLURM_ARRAY_TASK_ID} 

FASTA="$1"
KMERS_PREFIX="$2"
OUTDIR="$3"

./fastaRegex.to.bed.sh "$FASTA" "$KMERS_PREFIX"$(printf "%03d" ${N}) "$OUTDIR"
