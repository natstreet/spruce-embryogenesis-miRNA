#!/bin/bash -l

out=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/TAPIRhybrid
mkdir -p $out
cd $out

mir=../Picea_mature.fa
in=../PicabTx40

# seqkit split Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa -p 40 -O PicabTx40

for f in $(ls $in/*.fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
echo $n
sbatch -J $n ${PROJECT_DIR}/pipeline/runTapirHybrid.sh $mir $f All_tapir_hybrid_${n}.txt
done
