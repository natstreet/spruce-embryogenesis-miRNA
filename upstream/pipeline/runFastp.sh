#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 4
#SBATCH -t 12:00:00

module purge
module load GCC/12.3.0
module load GCC/13.2.0
module load fastp

f=$1
out=$2
adap=$3

## taking options from nf-core/smrnaseq
fastp \
  --in1 $f \
  --out1  $out/$(basename ${f/.fastq.gz//}).fastp.fastq.gz \
  --thread 4 \
  --json $out/$(basename ${f/.fastq.gz//}).fastp.json \
  --html $out/$(basename ${f/.fastq.gz//}).fastp.html \
  --adapter_fasta $adap \
  -l 17 --max_len1 40 \
  2> $out/$(basename ${f/.fastq.gz//}).fastp.log

