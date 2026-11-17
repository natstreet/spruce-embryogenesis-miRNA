#!/bin/bash -l

out=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/TAPIRhybrid

cd $out


sbatch -J 01 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_001.fa Fix_tapir_hybrid_001.fa.txt #Picea_mature_4334
sbatch -J 03 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_003.fa Fix_tapir_hybrid_003.fa.txt #Picea_mature_4742
sbatch -J 05 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_005.fa Fix_tapir_hybrid_005.fa.txt #Picea_mature_4190
sbatch -J 07 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_007.fa Fix_tapir_hybrid_007.fa.txt #Picea_mature_4649
sbatch -J 08 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_008.fa Fix_tapir_hybrid_008.fa.txt #Picea_mature_4626
sbatch -J 09 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_009.fa Fix_tapir_hybrid_009.fa.txt #Picea_mature_4433
sbatch -J 14 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_014.fa Fix_tapir_hybrid_014.fa.txt #Picea_mature_4586
sbatch -J 16 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_016.fa Fix_tapir_hybrid_016.fa.txt #Picea_mature_4671
sbatch -J 17 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_017.fa Fix_tapir_hybrid_017.fa.txt #Picea_mature_4740
sbatch -J 18 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_018.fa Fix_tapir_hybrid_018.fa.txt #Picea_mature_4511
sbatch -J 20 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_020.fa Fix_tapir_hybrid_020.fa.txt #Picea_mature_4655
sbatch -J 21 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_021.fa Fix_tapir_hybrid_021.fa.txt #Picea_mature_4503
sbatch -J 22 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_022.fa Fix_tapir_hybrid_022.fa.txt #Picea_mature_4050
sbatch -J 26 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_026.fa Fix_tapir_hybrid_026.fa.txt #Picea_mature_4535
sbatch -J 27 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_027.fa Fix_tapir_hybrid_027.fa.txt #Picea_mature_4697
sbatch -J 29 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_029.fa Fix_tapir_hybrid_029.fa.txt #Picea_mature_4395
sbatch -J 30 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames3500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_030.fa Fix_tapir_hybrid_030.fa.txt #Picea_mature_3902
sbatch -J 31 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames3500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_031.fa Fix_tapir_hybrid_031.fa.txt #Picea_mature_3927
sbatch -J 32 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames3000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_032.fa Fix_tapir_hybrid_032.fa.txt #Picea_mature_3061
sbatch -J 33 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_033.fa Fix_tapir_hybrid_033.fa.txt #Picea_mature_4283
sbatch -J 34 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames3500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_034.fa Fix_tapir_hybrid_034.fa.txt #Picea_mature_3520
sbatch -J 35 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_035.fa Fix_tapir_hybrid_035.fa.txt #Picea_mature_4033
sbatch -J 36 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames3500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_036.fa Fix_tapir_hybrid_036.fa.txt #Picea_mature_3611
sbatch -J 37 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_037.fa Fix_tapir_hybrid_037.fa.txt #Picea_mature_4697
sbatch -J 38 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4500.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_038.fa Fix_tapir_hybrid_038.fa.txt #Picea_mature_4815
sbatch -J 40 ${PROJECT_DIR}/pipeline/runTapirHybrid.sh ../Picea_mature_shortnames4000.fa ../PicabTx40/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_040.fa Fix_tapir_hybrid_040.fa.txt #Picea_mature_4050
