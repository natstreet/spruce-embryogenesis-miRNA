#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 4
#SBATCH --mem=50GB
#SBATCH -t 168:00:00
#SBATCH -J target

export PATH="${TOOLS_DIR}/fasta-35.4.12/bin:$PATH"
targetfinder=${TOOLS_DIR}/TargetFinder/targetfinder.pl

TMPDIR=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/TargetFinder
out=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/TargetFinder
in=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picea_mature_shortnames842.fa
#ref=${PROJECT_DIR}/data/analysis/results/targetPrediction_all/Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa

module purge
module load GCCcore/14.3.0
module load Perl/5.40.2
module load GCCcore/13.2.0
module load BioPerl/1.7.8

mkdir -p $out
#mkdir -p $TMPDIR
cd $out

# Parse FASTA and run targetfinder
awk '
    /^>/ { 
        if (name && seq) {
            print name "\t" seq
        }
        name = substr($0, 2)
        sub(/ .*/, "", name)
        seq = ""
        next
    }
    { seq = seq $0 }
    END {
        if (name && seq) {
            print name "\t" seq
        }
    }
' "$in" | while IFS=$'\t' read -r name seq; do
    echo "Processing $name $seq ..."
    $targetfinder -t 4 -s "$seq" -d ../Picab02_230926_at01_longest_no_TE_sorted_mRNA.fa -q "$name" > "$out/${name}_predicted_targets.txt"
done

# this works, cant submit but can run in bash in Rstudio ondemand (??)


# need to use relative path for ref

# # Try chunk
# for f in $(ls ../PicabTx40/*.fa)
# do
# n=$(basename ${f/Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_//})
# echo $n
# awk '
#     /^>/ {
#         if (name && seq) {
#             print name "\t" seq
#         }
#         name = substr($0, 2)
#         sub(/ .*/, "", name)
#         seq = ""
#         next
#     }
#     { seq = seq $0 }
#     END {
#         if (name && seq) {
#             print name "\t" seq
#         }
#     }
# ' "$in" | while IFS=$'\t' read -r name seq; do
#     echo "Processing $name $seq ..."
#     $targetfinder -t 4 -s "$seq" -d $f -q "$name" > "${name}_targetfinder_${n}.txt"
# done
# echo "done"
# done
# 
# awk '
#     /^>/ {
#         if (name && seq) {
#             print name "\t" seq
#         }
#         name = substr($0, 2)
#         sub(/ .*/, "", name)
#         seq = ""
#         next
#     }
#     { seq = seq $0 }
#     END {
#         if (name && seq) {
#             print name "\t" seq
#         }
#     }
# ' "$in" | while IFS=$'\t' read -r name seq; do
#     echo "Processing $name $seq ..."
#     $targetfinder -t 4 -s "$seq" -d Picab02_230926_at01_longest_no_TE_sorted_mRNA.part_001.fa -q "$name" > "${name}_targetfinder.txt"
# done
