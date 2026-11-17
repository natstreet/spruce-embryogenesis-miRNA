#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 8
#SBATCH --mem=50GB
#SBATCH -t 168:00:00
#SBATCH -J tapir

tapir="apptainer exec -B ${BIND_DIR}:${BIND_DIR} ${CONTAINER_DIR}/tapir.sif"

out=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/TAPIRfasta
#in=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picea_mature.fa
#ref=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa

# ml SeqKit/2.8.2
# seqkit split Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa -p 40 -O PicabTx40

mkdir -p $out
cd $out

for f in $(ls ../PicabTx40/*.fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
echo $n
  for g in $(ls ../miR10/*.fa)
  do
  m=$(basename ${g/Picea_mature_shortnames.part_//})
  echo $m
  $tapir tapir_fasta --mir_file $g --target_file $f > miR${m}_tapir_fasta_${n}.txt
  echo "done"
  done
echo "done"
done
