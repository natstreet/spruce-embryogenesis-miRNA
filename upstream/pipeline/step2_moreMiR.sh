#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J blast

# folders
res=${PROJECT_DIR}/data/analysis/results/blast_moreMiR
index=${PROJECT_DIR}/data/analysis/KnownMiR
mrdp2pip=${TOOLS_DIR}/miRDeep-P2_pipeline/script
 
module purge
module load GCC/13.2.0 OpenMPI/4.1.6
module load BLAST+/2.16.0 Perl/5.38.0
module load BioPerl/1.7.8

# # Combine all db, add DB name to the seq name because some seq has duplicate name but not duplicate sequences (why...)
# cd $index
# cat miRBase_mature.fa | sed 's/^>/>miRBase_/' > all_mature.fa
# cat PmiREN_mature.fa | sed 's/^>/>PmiREN_/' >> all_mature.fa
# cat PMRD_mature.fa | sed 's/^>/>PMRD_/' >> all_mature.fa
# 
# # Index the known miRNA
# perl $mrdp2pip/fasta_U2T.pl $index/all_mature.fa $index/mature_wo_U.fa
# perl $mrdp2pip/unique_fasta_v1.3.pl $index/mature_wo_U.fa $index/mature_wo_U_uniq.fa mature_uniq
# makeblastdb -in $index/mature_wo_U_uniq.fa -dbtype 'nucl'
# 
# # Combine miRDP2 predictions
# # mkdir $res
cd $res
# #realpath */*_predictions > predictions_list.txt
# #perl $mrdp2pip/parse_miRDP2_prediction.pl predictions_list.txt Picea
# cp ${PROJECT_DIR}/data/analysis/results/Picea_mature.fa .

# # Blast
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn.txt \
#        -word_size 4 -num_alignments 1 -num_threads 1
# #set numthread to 1 to see if we can make the result more reproducible. Right now different runs give different number of known/variant

# # try more hits (num_align 20), because even with num thread 1 it still not reproducible
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_20.txt \
#        -word_size 4 -num_alignments 20 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_20.txt \
#         Picea_mature_blastn_20_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_20_parsed.txt \
#         Picea_mature_20

# # it now seems reproducible. I'll try 10 and 100 num_alignjust to benchmark
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_10.txt \
#        -word_size 4 -num_alignments 10 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_10.txt \
#         Picea_mature_blastn_10_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_10_parsed.txt \
#         Picea_mature_10
# 
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_100.txt \
#        -word_size 4 -num_alignments 100 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_100.txt \
#         Picea_mature_blastn_100_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_100_parsed.txt \
#         Picea_mature_100

# 20 is not enough, I'll try 50, 500, 1000
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_50.txt \
#        -word_size 4 -num_alignments 50 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_50.txt \
#         Picea_mature_blastn_50_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_50_parsed.txt \
#         Picea_mature_50
# 
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_500.txt \
#        -word_size 4 -num_alignments 500 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_500.txt \
#         Picea_mature_blastn_500_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_500_parsed.txt \
#         Picea_mature_500
#         
# blastn -db $index/mature_wo_U_uniq.fa \
#        -query Picea_mature.fa \
#        -out Picea_mature_blastn_1000.txt \
#        -word_size 4 -num_alignments 1000 -num_threads 1
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn_1000.txt \
#         Picea_mature_blastn_1000_parsed.txt
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_1000_parsed.txt \
#         Picea_mature_1000
# known miR saturate at 50 but variant miR saturate at 500

# try default (250)
blastn -db $index/mature_wo_U_uniq.fa \
       -query Picea_mature.fa \
       -out Picea_mature_blastn_250.txt \
       -word_size 4 -num_alignments 250 -num_threads 1
perl $mrdp2pip/general_blast_parser.pl \
        Picea_mature_blastn_250.txt \
        Picea_mature_blastn_250_parsed.txt
perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
        Picea_mature_blastn_250_parsed.txt \
        Picea_mature_250

# # Annotate if the miRDP2 predictions are known miRNAs
# ## I comment out use Bio::DB::GenBank; from the header of this file
# perl $mrdp2pip/general_blast_parser.pl \
#         Picea_mature_blastn.txt \
#         Picea_mature_blastn_parsed.txt
# 
# perl $mrdp2pip/parse_parsed_blast_known_plants.pl \
#         Picea_mature_blastn_parsed.txt \
#         Picea_mature
# 
# cut -f2 Picea_mature_known.txt | sort | uniq > Picea_mature_known_id.txt
# cut -f2 Picea_mature_variant.txt | sort | uniq > temp.txt
# 
# perl $mrdp2pip/filter_lines_by_key_words_list.pl \
#         temp.txt \
#         Picea_mature_variant_id.txt \
#         Picea_mature_known_id.txt 0
# 
# rm temp.txt
