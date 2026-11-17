#!/bin/bash -l

# folders
in=${PROJECT_DIR}/data/analysis/trimmed
out=${PROJECT_DIR}/data/analysis/sortmerna

ref=${REFERENCE_DIR}/smr_with_Picab_rRNAsnoRNAtRNA.fa

for f in $(ls $in/*.fastq.gz)
do
echo $f
#echo $(basename ${f/.fastq.gz//})
sbatch ./runSortMeRNA.sh $f $out $ref
done
