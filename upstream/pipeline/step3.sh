#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 4
#SBATCH --mem=100GB
#SBATCH -t 72:00:00
#SBATCH -J align

in=${PROJECT_DIR}/data/analysis/mergeTechrep
res=${PROJECT_DIR}/data/analysis/results
mrdp2pip=${TOOLS_DIR}/miRDeep-P2_pipeline/script

module purge
module load GCC/12.3.0 Perl/5.38.0
module load Bowtie/1.3.1 SAMtools/1.18 BioPerl/1.7.8

cd $res
bowtie-build Picea_mature.fa miRNA_ref

mkdir -p $res/bowtie
for f in $(ls $in/*.fq.gz)
do
echo $(basename ${f/.fq.gz//})
bowtie -p 4 -v 0 --norc -S -x miRNA_ref $f | samtools view -Sb - > bowtie/$(basename ${f/.fq.gz//}).bam
perl $mrdp2pip/bam2ref_counts.pl -bam bowtie/$(basename ${f/.fq.gz//}).bam -f Picea_mature.fa > bowtie/$(basename ${f/.fq.gz//})_count.txt
done

cd $res/bowtie
realpath *_count.txt > count_list.txt
perl $mrdp2pip/combine_htseq_counts.pl count_list.txt count_table.txt
# This last perl step took suspiciously long when submitted as slurm job, I ran it separately in terminal it was quickly done

# Remove the last four columns (mean, median, variance and coefficient of variation) before DE analysis
awk 'BEGIN {OFS="\t"} {NF -= 4; print}' count_table.txt > count_DE.txt
