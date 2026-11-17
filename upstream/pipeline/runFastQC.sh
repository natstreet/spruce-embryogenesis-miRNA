#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 1
#SBATCH -t 72:00:00

module purge
module load FastQC/0.12.1-Java-11

f=$1
out=$2

fastqc --outdir $out -t 1 --noextract $f

