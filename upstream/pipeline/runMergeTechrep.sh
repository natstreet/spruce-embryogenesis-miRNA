#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 1
#SBATCH -t 12:00:00

# folders
in=${PROJECT_DIR}/data/analysis/sortmerna
out=${PROJECT_DIR}/data/analysis/mergeTechrep

mkdir -p $out

# SE
sampsheet=${PROJECT_DIR}/doc/U.Egertsdotter_17_02_sample_info.txt
for f in $(sed 1d $sampsheet | cut -f1)
do
echo $f
ls $in/*/${f}_*.fastp_sortmerna.fq.gz
cat $in/*/${f}_*.fastp_sortmerna.fq.gz > $out/${f}.fq.gz
done

# ZE
sampsheet=${PROJECT_DIR}/doc/U.Egertsdotter_18_02_sample_info.txt
for f in $(sed 1d $sampsheet | cut -f1)
do
echo $f
ls $in/*/${f}_*.fastp_sortmerna.fq.gz
cat $in/*/${f}_*.fastp_sortmerna.fq.gz > $out/${f}.fq.gz
done
