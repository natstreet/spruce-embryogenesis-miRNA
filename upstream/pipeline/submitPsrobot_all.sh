#!/bin/bash -l

out=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/PSROBOT
#in=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picea_mature.fa
#ref=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa

# ml SeqKit/2.8.2
# seqkit split Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa -p 40 -O PicabTx40

# mkdir -p $out
cd $out

# it seems to have problem when target file is 2 folders away with a lot of . ?

for f in $(ls ../PicabTx40/*_01[89].fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
cp $f ./target_${n}
sbatch -J $n ${PROJECT_DIR}/pipeline/runPsrobot.sh ../Picea_mature_shortnames.fa target_${n} All_psrobot_${n}.txt
done

for f in $(ls ../PicabTx40/*_02[0123456789].fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
cp $f ./target_${n}
sbatch -J $n ${PROJECT_DIR}/pipeline/runPsrobot.sh ../Picea_mature_shortnames.fa target_${n} All_psrobot_${n}.txt
done

for f in $(ls ../PicabTx40/*_03[0123456789].fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
cp $f ./target_${n}
sbatch -J $n ${PROJECT_DIR}/pipeline/runPsrobot.sh ../Picea_mature_shortnames.fa target_${n} All_psrobot_${n}.txt
done

for f in $(ls ../PicabTx40/*_040.fa)
do
n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
cp $f ./target_${n}
sbatch -J $n ${PROJECT_DIR}/pipeline/runPsrobot.sh ../Picea_mature_shortnames.fa target_${n} All_psrobot_${n}.txt
done
