#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J idx
#SBATCH --mem=100GB

module load GCC/12.3.0
module load Bowtie/1.3.1

ref=${REFERENCE_DIR}/Picab02_chromosomes_and_unplaced.fa.gz
refidx=${REFERENCE_DIR}/bowtieIdx/Picab02_chromosomes_and_unplaced_index
bowtie-build --threads 8 $ref $refidx

mirdpidx=${TOOLS_DIR}/miRDP2-v1.1.4/scripts/index
bowtie-build --threads 8 $mirdpidx/ncRNA_rfam.fa $mirdpidx/rfam_index

# mirdpidx=${TOOLS_DIR}/miRDP2-v1.1.4/scripts/index
# bowtie-inspect $mirdpidx/mature_index > $mirdpidx/mature.fa

# module purge
# module load GCC/13.2.0
# module load Bowtie2/2.5.4

# ref=${REFERENCE_DIR}/Picab02_chromosomes_and_unplaced.fa.gz
# refidx=${REFERENCE_DIR}/bowtie2Idx/Picab02_chromosomes_and_unplaced_index
# bowtie2-build --threads 8 $ref $refidx

# bowtie2-build --threads 8 $mirdpidx/mature.fa $mirdpidx/mature_index
# bowtie2-build --threads 8 $mirdpidx/ncRNA_rfam.fa $mirdpidx/rfam_index
