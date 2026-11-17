#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 10
#SBATCH --mem=50GB
#SBATCH -t 12:00:00

f=$1
out=$2
ref=$3

module purge

mkdir $out/$(basename ${f/.fastp.fastq.gz//})

apptainer exec -B ${BIND_DIR}:${BIND_DIR} ${CONTAINER_DIR}/nf-core-sortmerna-4.3.7--6502243397c065ba.img \
  sortmerna --ref $ref --reads $f --workdir $out/$(basename ${f/.fastp.fastq.gz//}) --other $out/$(basename ${f/.fastp.fastq.gz//})/$(basename ${f/.fastq.gz//})_sortmerna --fastx --threads 10

