#!/bin/bash -l

set -eux
export SINGULARITY_BINDPATH="/mnt:/mnt"

#path
proj=kaw
salmon=${CONTAINER_DIR}/salmon.sif

#spruce
salmon_index=${REFERENCE_DIR}/Picea-abies/v2.0/indices/salmon/Picab02_230926_at01_longest_no_TE_sorted_mRNA_salmon-version-1-dot-10-dot-3
in=${PROJECT_DIR}/data/mRNA/ZE_RNA_fastq
outdir=${PROJECT_DIR}/data/mRNA/TestLongestTx_ZEsalmon

[[ ! -d $outdir ]] && mkdir -p $outdir

for f in $in/*_1.fq.gz ; do
 sbatch -A <account> -p main ./runSalmon.sh \
   $salmon $salmon_index $f ${f/_1.fq.gz/_2.fq.gz} $outdir
done
