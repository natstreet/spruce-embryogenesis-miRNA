#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 1
#SBATCH --mem=50GB
#SBATCH -t 24:00:00
#SBATCH -J dedup

# folders
in=${PROJECT_DIR}/data/analysis/mergeTechrep
out=${PROJECT_DIR}/data/analysis/fasta
scripts=${PROJECT_DIR}

module purge
module load GCC/13.2.0  OpenMPI/4.1.6
module load Biopython/1.84

# deduplication using the modified version of Ahsan's script
# now it read fastq directly without the need of seqtk
for f in $(ls $in/*.fq.gz)
do
echo $(basename ${f/.fq.gz//})
python $scripts/pipeline/preprocessing_for_mirDP2_fastqinput.py $f
done

# moving output files away
mkdir -p $out
mv $in/*.dedup.fa $out/

# make filelist for mrdp2
realpath $out/*.dedup.fa > $scripts/doc/fastalist.txt
