#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 8
#SBATCH --mem=100GB
#SBATCH -t 168:00:00
#SBATCH -J mrdp

# Bowtie version
# ref=${REFERENCE_DIR}/Picab02_chromosomes_and_unplaced.fa.gz
ref=${REFERENCE_DIR}/Picab02_chromosomes_and_unplaced.fa
refidx=${REFERENCE_DIR}/bowtieIdx/Picab02_chromosomes_and_unplaced_index
mrdp2=${TOOLS_DIR}/miRDP2-v1.1.4
scripts=${PROJECT_DIR}
out=${PROJECT_DIR}/data/analysis/results

module purge
module load GCC/12.3.0
module load Bowtie/1.3.1
module load Perl/5.38.0

mkdir -p $out
cd $out
$mrdp2/miRDP2-v1.1.4_pipeline.bash \
  -g $ref -x $refidx -f -b $scripts/doc/fastalist.txt \
  --large-index -o $out \
  -p 8

# bowtie2 version
# ref=${REFERENCE_DIR}/Picab02_chromosomes_and_unplaced.fa
# refidx=${REFERENCE_DIR}/bowtie2Idx/Picab02_chromosomes_and_unplaced_index
# mrdp2=${TOOLS_DIR}/miRDP2-v1.1.4
# out=${PROJECT_DIR}/data/results2
# 
# module load GCC/13.2.0
# module load Bowtie2/2.5.4
# module load Perl/5.38.0
# 
# cd $out
# $mrdp2/miRDP2-v1.1.4_pipeline.bash \
#   -g $ref -x $refidx -f -b ${PROJECT_DIR}/doc/fastalist.txt \
#   --large-index -o $out \
#   --bowtie2 -p 8
