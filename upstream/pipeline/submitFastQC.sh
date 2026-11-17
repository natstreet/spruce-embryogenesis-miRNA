#!/bin/bash -l


out=${PROJECT_DIR}/data/analysis/fastqc
mkdir -p $out

#raw
in=${PROJECT_DIR}/data/rawUgertsdotter
mkdir -p $out/raw
for f in $(find $in -name "*_R1_001.fastq.gz")
do
echo $f
sbatch ./runFastQC.sh $f $out/raw
done

#trim
in=${PROJECT_DIR}/data/analysis/trimmed
mkdir -p $out/trimmed
for f in $(find $in -name "*.fastq.gz")
do
echo $f
sbatch ./runFastQC.sh $f $out/trimmed
done

#remove contaminents
in=${PROJECT_DIR}/data/analysis/sortmerna
mkdir -p $out/sortme
for f in $(find $in -name "*_sortmerna.fq.gz")
do
echo $f
sbatch ./runFastQC.sh $f $out/sortme
done

#merge techrep
in=${PROJECT_DIR}/data/analysis/mergeTechrep
mkdir -p $out/mergedtechrep
for f in $(find $in -name "*.fq.gz")
do
echo $f
sbatch ./runFastQC.sh $f $out/mergedtechrep
done
