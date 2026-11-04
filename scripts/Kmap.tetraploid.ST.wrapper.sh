#!/bin/bash
#SBATCH --array=0-190
#SBATCH --partition=Standard.2.0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rdalipo1@binghamton.edu
#SBATCH --cpus-per-task=1

# Run from the split files folder

N=${SLURM_ARRAY_TASK_ID}

~/scripts/Kmap.tetraploid.ST.R ~/v10.XlaXtr.homology/Xla/Xenla10.chr.size.noscaff.txt split.stcount.$(printf "%03d" ${N}) 20000000 1908673 ~/v10.XlaXtr.homology/Xla/Xenla10.Subgenome.txt Xla.st13.tukey.split.$(printf "%03d" ${N}).tbl Xla.st13.norm.$(printf "%03d" ${N}).tbl