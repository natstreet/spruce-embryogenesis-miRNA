# Notes and code when I found out that the list of known miR changes when rerun blastn

library(tidyverse)
library(here)

# known miR
knownSeq <- seqinr::read.fasta(here("data/analysis/KnownMiR/mature_wo_U_uniq.fa"), 
                               seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
  map(as.character) %>% 
  map(as.data.frame) %>% 
  bind_rows(.id = "Hit name") %>% 
  dplyr::rename(knownSequence = `.x[[i]]`) %>% 
  mutate(Description = seqinr::read.fasta(here("data/analysis/KnownMiR/mature_wo_U_uniq.fa"), 
                                          seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
           map_chr(seqinr::getAnnot) %>% 
           str_split_i(., pattern = " ",i = 2))

PmiREN <- names(seqinr::read.fasta(here("data/analysis/KnownMiR/PmiREN_mature.fa"), 
                                   seqtype = "DNA", as.string = T, forceDNAtolower = F))

PMRD <- names(seqinr::read.fasta(here("data/analysis/KnownMiR/PMRD_mature.fa"), 
                                 seqtype = "DNA", as.string = T, forceDNAtolower = F))

miRBase <- names(seqinr::read.fasta(here("data/analysis/KnownMiR/miRBase_mature.fa"), 
                                    seqtype = "DNA", as.string = T, forceDNAtolower = F))

stopifnot(length(intersect(miRBase, PmiREN)) == 0)
stopifnot(length(intersect(PmiREN, PMRD)) == 0)

forCheck <- intersect(miRBase, PMRD)
PMRDseq <- seqinr::read.fasta(here("data/analysis/KnownMiR/PMRD_mature.fa"), 
                              seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
  map(as.character) %>% 
  keep_at(forCheck)
miRBaseseq <- seqinr::read.fasta(here("data/analysis/KnownMiR/miRBase_mature.fa"), 
                                 seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
  map(as.character) %>% 
  keep_at(forCheck)

wrong <- imap(miRBaseseq, ~ (.x == PMRDseq[[.y]])) %>% 
  discard(~ .x)
PMRDseq <- PMRDseq %>%
  keep_at(names(wrong))
miRBaseseq <- miRBaseseq %>%
  keep_at(names(wrong))

mature_known <- read_lines(here("data/analysis/results/blast_moreMiR/Picea_mature_known_id.txt"))
mature_variant <- read_lines(here("data/analysis/results/blast_moreMiR/Picea_mature_variant_id.txt"))
mature_blast <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_blastn_parsed.txt"), show_col_types = F) %>%
  filter(Query_acc %in% c(mature_known,mature_variant)) %>% 
  left_join(knownSeq, by = c("Hit name", "Description")) %>% 
  select(Query_acc, Description, knownSequence) %>% 
  mutate(Description = gsub("&", "|", gsub("&$","",Description)))

sum(names(wrong) %in% unlist(str_split(mature_blast$Description, "\\|"), use.name = FALSE)) #5
# I will need to redo the annotation again... didn't expect miRNAs with the same name from different database will have different sequences

# 15 July, redo annotation : 2 less known??
# check known IDs
# the additional variant Picea_mature_87 matched to pab miR, though not vary highly match but being pab I think is's worth keeping this version
mature_known <- read_lines(here("data/analysis/results/blast_moreMiR/Picea_mature_known_id.txt"))
mature_variant <- read_lines(here("data/analysis/results/blast_moreMiR/Picea_mature_variant_id.txt"))
mature_known_old <- read_lines(here("data/analysis/results/blast_moreMiR/old/Picea_mature_known_id.txt"))
mature_variant_old <- read_lines(here("data/analysis/results/blast_moreMiR/old/Picea_mature_variant_id.txt"))

mature_known_old[!mature_known_old %in% mature_known]
# [1] "Picea_mature_2078" "Picea_mature_3428" "Picea_mature_3443" "Picea_mature_3859"
# [5] "Picea_mature_4177" "Picea_mature_4322" "Picea_mature_4336"
mature_variant_old[!mature_variant_old %in% mature_variant]
# [1] "Picea_mature_1102" "Picea_mature_1195" "Picea_mature_1464" "Picea_mature_2303"
# [5] "Picea_mature_3401" "Picea_mature_447"  "Picea_mature_765" 
mature_known[!mature_known %in% mature_known_old]
# [1] "Picea_mature_1195" "Picea_mature_1742" "Picea_mature_2303" "Picea_mature_3401"
# [5] "Picea_mature_765" 
mature_variant[!mature_variant %in% mature_variant_old]
# [1] "Picea_mature_2078" "Picea_mature_2183" "Picea_mature_2867" "Picea_mature_3066"
# [5] "Picea_mature_3146" "Picea_mature_3428" "Picea_mature_3443" "Picea_mature_3859"
# [9] "Picea_mature_4113" "Picea_mature_4177" "Picea_mature_4322" "Picea_mature_4336"
# [13] "Picea_mature_4338" "Picea_mature_4485" "Picea_mature_4494" "Picea_mature_87" 

mature_blast <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_blastn_parsed.txt"), show_col_types = F)
mature_blast_old <- read_tsv(here("data/analysis/results/blast_moreMiR/old/Picea_mature_blastn_parsed.txt"), show_col_types = F)
# Picea_mature_1195 variant->known
# Picea_mature_1742 missing->known
# Picea_mature_2303 variant->known
# Picea_mature_3401 variant->known
# Picea_mature_765  variant->known
# Picea_mature_1102 variant->missing
# Picea_mature_1464 variant->missing
# Picea_mature_447 variant->missing
# Picea_mature_2078 known->variant
# Picea_mature_3428 known->variant
# Picea_mature_3443 known->variant
# Picea_mature_3859 known->variant
# Picea_mature_4177 known->variant
# Picea_mature_4322 known->variant
# Picea_mature_4336 known->variant
# Picea_mature_2183 missing->variant
# Picea_mature_2867 missing->variant
# Picea_mature_3066 missing->variant
# Picea_mature_3146 missing->variant
# Picea_mature_4113 missing->variant
# Picea_mature_4338 missing->variant
# Picea_mature_4485 missing->variant
# Picea_mature_4494 missing->variant
# Picea_mature_87 missing->variant
# The list is not consistant when I run blastn again

# check if the extremely long names effect the result
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/experiments/Picea_mature_shortnames_20_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/experiments/Picea_mature_shortnames_20_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
# as expected, nope. good!
# also, this means that the 2 runs of -num_alignments 20 -num_threads 1 give similar and reproducible results

# check if the known list different between num_align 1 and num_align 20
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #FALSE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
length(unique(temp1$X2)) #342 known
length(unique(temp2$X2)) #370 known
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp3$X2, temp3$X4) == paste0(temp4$X2, temp4$X4)) #FALSE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #TRUE
length(unique(setdiff(temp3$X2,temp1$X2))) #385 variant
length(unique(setdiff(temp4$X2,temp2$X2))) #384 variant
all(mature_known %in% temp2$X2) #TRUE
all(mature_known_old %in% temp2$X2) #TRUE
# So the num_align 20 contains all hits that were included in num_align 1, and 28 more. And now that include all known miRNAs in all the runs we made, both old and new

# benchmark with num_align 10 and 100
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_10_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_10_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp3$X2, temp3$X4) == paste0(temp4$X2, temp4$X4)) #FALSE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #TRUE

temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #FALSE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #FALSE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_20_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp3$X2, temp3$X4) == paste0(temp4$X2, temp4$X4)) #FALSE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #FALSE
temp1$X2[!temp1$X2 %in% temp2$X2] #Picea_mature_2610
# Okay, 20 is not enough apparently

# benchmark with num_align 50, 500 and 1000
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_50_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_50_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp4$X2, temp4$X4) %in% paste0(temp3$X2, temp3$X4)) #FALSE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #TRUE

temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_100_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp4$X2, temp4$X4) %in% paste0(temp3$X2, temp3$X4)) #TRUE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #FALSE

temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_1000_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_1000_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp4$X2, temp4$X4) %in% paste0(temp3$X2, temp3$X4)) #TRUE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #TRUE
# 50 is okay in term of known miR, but variant miR saturate at 500

# try 250 (default)
temp1 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_known.txt"), show_col_types = FALSE, col_names = FALSE)
temp2 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_known.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp1$X2, temp1$X4) == paste0(temp2$X2, temp2$X4)) #TRUE
all(paste0(temp1$X2, temp1$X4) %in% paste0(temp2$X2, temp2$X4)) #TRUE
temp3 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variant.txt"), show_col_types = FALSE, col_names = FALSE)
temp4 <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_500_variant.txt"), show_col_types = FALSE, col_names = FALSE)
all(paste0(temp4$X2, temp4$X4) %in% paste0(temp3$X2, temp3$X4)) #TRUE
all(paste0(temp3$X2, temp3$X4) %in% paste0(temp4$X2, temp4$X4)) #TRUE
# 250 is fine

# Final code
mature_blast <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_blastn_250_parsed.txt"), show_col_types = F)
knownList <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_known.txt"), 
                      show_col_types = FALSE, col_names = colnames(mature_blast))
variantList <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variant.txt"), 
                        show_col_types = FALSE, col_names = colnames(mature_blast)) %>% 
  filter(!Query_acc %in% knownList$Query_acc)

knownList %>% 
  select(Query_acc, `Hit name`, Description) %>% 
  mutate(Description = gsub("&$", "", Description)) %>% 
  write_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_knownID.txt"), col_names = FALSE)

dat <- variantList %>% 
  split(., .$Query_acc) %>% 
  map(slice_max, `Query cover`*`Hit cover`*(Identity/100))

full_join(
  dat %>% map(pull, `Hit name`) %>% 
    map(paste, collapse = "|") %>% 
    map(as.data.frame) %>% 
    bind_rows(.id = "miRNA"),
  dat %>% map(pull, Description) %>% 
    map(~ gsub("&$", "", .x)) %>% 
    map(paste, collapse = "|") %>% 
    map(as.data.frame) %>% 
    bind_rows(.id = "miRNA"), 
  by = "miRNA") %>% 
  write_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variantID.txt"), col_names = FALSE)
