#!/bin/bash -l

in=${PROJECT_DIR}/data/rawUgertsdotter
out=${PROJECT_DIR}/data/analysis/trimmed
adap=${PROJECT_DIR}/doc/nfcore_known_adapters.fa
# adapters taken from nf-core smrnaseq pipeline v2.4.0

mkdir -p $out

# Use only forward reads
for f in $(find $in -name "*_R1_001.fastq.gz")
do
echo $f
#echo $(basename ${f/.fastq.gz//})
sbatch ./runFastp.sh $f $out $adap
done
