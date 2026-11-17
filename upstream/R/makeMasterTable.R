# table all miRNA
# - known, which known it match to
# - plant
# - DE, which transition
# - targets (predicted)
# - targets with support, opposite De (which transition), neg correlation

library(tidyverse)
library(here)
library(pwalign)

# sequence
matureSeq <- seqinr::read.fasta(here("data/analysis/results/Picea_mature.fa"), 
                                seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
  map(as.character)

# plant criteria
passed_plant_criteria_id <- readRDS(here("data/analysis/results/passed_plant_criteria.rds")) %>%
  map(function(d){
    keep_at(matureSeq, at = matureSeq %in% d$mature)
  }) %>%
  map(names)

mature_plant <- unique(unlist(passed_plant_criteria_id, use.names = F))

# known

# mature_blast <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_blastn_250_parsed.txt"), show_col_types = F)
# knownList <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_known.txt"), 
#                       show_col_types = FALSE, col_names = colnames(mature_blast))
# variantList <- read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variant.txt"), 
#                         show_col_types = FALSE, col_names = colnames(mature_blast)) %>% 
#   filter(!Query_acc %in% knownList$Query_acc)
# 
# knownList %>% 
#   select(Query_acc, `Hit name`, Description) %>% 
#   mutate(Description = gsub("&$", "", Description)) %>% 
#   write_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_knownID.txt"), col_names = FALSE)
# 
# dat <- variantList %>% 
#   split(., .$Query_acc) %>% 
#   map(slice_max, `Query cover`*`Hit cover`*(Identity/100))
# 
# full_join(
#   dat %>% map(pull, `Hit name`) %>% 
#     map(paste, collapse = "|") %>% 
#     map(as.data.frame) %>% 
#     bind_rows(.id = "miRNA"),
#   dat %>% map(pull, Description) %>% 
#     map(~ gsub("&$", "", .x)) %>% 
#     map(paste, collapse = "|") %>% 
#     map(as.data.frame) %>% 
#     bind_rows(.id = "miRNA"), 
#   by = "miRNA") %>% 
#   write_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variantID.txt"), col_names = FALSE)
# known miR
knownSeq <- seqinr::read.fasta(here("data/analysis/KnownMiR/mature_wo_U_uniq.fa"), 
                               seqtype = "DNA", as.string = T, forceDNAtolower = F) %>% 
  map(as.character)

annot <- bind_rows(
  read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_knownID.txt"), 
           show_col_types = FALSE, col_names = c("id","uniqname","description")) %>% 
    mutate(annotation = "known"),
  read_tsv(here("data/analysis/results/blast_moreMiR/Picea_mature_250_variantID.txt"), 
           show_col_types = FALSE, col_names = c("id","uniqname","description")) %>% 
    mutate(annotation = "variant")) %>% 
  rowwise() %>% 
  mutate(database = list(c("miRBase","PMRD", "PmiREN")[c(grepl("miRBase",description),grepl("PMRD",description),grepl("PmiREN",description))]) %>% 
           map_chr(paste0, collapse="|"),
         knownSeq = map_chr(uniqname, function(u){
           names <- strsplit(u, split="\\|") %>% 
             unlist(use.names = FALSE)
           keep_at(knownSeq, names) %>% 
             paste0(collapse="|")
           })) %>% 
  ungroup() %>% 
  select(id, annotation, database, description, knownSeq)

# sequence data
matureDF <- as.data.frame(matureSeq) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::rename(sequence = V1) %>% 
  mutate(P_criteria = ifelse(id %in% mature_plant, "Plant", NA)) %>% 
  left_join(annot, by = "id") %>% 
  mutate(annotation = ifelse(is.na(annotation),"novel", annotation),
         family = str_extract_all(description, "miRN?\\d+") %>% 
           map(unique) %>% 
           map_chr(paste, collapse = "|") %>% 
           gsub("NA","",.))
# write_tsv(matureDF, here("data/analysis/results/allPredicted.tsv"))

# remake family columns prioritizing miRBase IDs
matureDF <- as.data.frame(matureSeq) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::rename(sequence = V1) %>% 
  mutate(P_criteria = ifelse(id %in% mature_plant, "Plant", NA)) %>% 
  left_join(annot, by = "id") %>% 
  mutate(annotation = ifelse(is.na(annotation),"novel", annotation),
         family_miRBase = str_extract_all(description, "miRBase_[a-zA-Z]{3}-miRN?\\d+") %>% 
           map(str_extract_all,"miRN?\\d+") %>% 
           map(unique) %>% 
           map_chr(paste, collapse = "|") %>% 
           gsub("NA","",.),
         family_PmiREN = str_extract_all(description, "PmiREN_[a-zA-Z]{3}-miRN?\\d+") %>% 
           map(str_extract_all,"miRN?\\d+") %>% 
           map(unique) %>% 
           map_chr(paste, collapse = "|") %>% 
           gsub("NA","",.),
         family_PMRD = str_extract_all(description, "PMRD_[a-zA-Z]{3}-miRN?\\d+") %>% 
           map(str_extract_all,"miRN?\\d+") %>% 
           map(unique) %>% 
           map_chr(paste, collapse = "|") %>% 
           gsub("NA","",.)) %>% 
  rowwise() %>% 
  mutate(family = ifelse(! family_miRBase == "", family_miRBase, 
                         paste(setdiff(unique(c(family_PmiREN, family_PMRD)),""), collapse = "|")))

# remake family columns prioritizing miRBase IDs AND keep 3p/5p/*
extract_family <- function(description, db_prefix) {
  pattern <- paste0(db_prefix, "_([a-zA-Z]{3})-(miRN?\\d+)([a-z]?)(-3p|-5p|\\*)?")
  
  str_match_all(description, pattern) %>%
    map_chr(function(m) {
      if (nrow(m) == 0) return("")
      
      parsed <- as.data.frame(m, stringsAsFactors = FALSE) %>%
        setNames(c("full", "species", "num", "locus", "arm")) %>%
        mutate(
          locus  = ifelse(is.na(locus), "", locus),
          arm    = ifelse(is.na(arm), "", arm),
          is_pab = species %in% c("Pab", "pab")
        ) %>%
        # If a no-arm version of (num, locus) exists, drop arm-tagged versions of it
        group_by(num, locus) %>%
        filter(!(any(arm == "") & arm != "")) %>%
        ungroup() %>%
        # If a pab version of (num, arm) exists, drop non-pab versions
        group_by(num, arm) %>%
        filter(!(any(is_pab) & !is_pab)) %>%
        ungroup() %>%
        mutate(id = ifelse(is_pab,
                           paste0(species, "-", num, locus, arm),
                           paste0(num, arm)))
      
      paste(unique(parsed$id), collapse = "|")
    })
}

matureDF <- as.data.frame(matureSeq) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::rename(sequence = V1) %>% 
  mutate(P_criteria = ifelse(id %in% mature_plant, "Plant", NA)) %>% 
  left_join(annot, by = "id") %>% 
  mutate(annotation = ifelse(is.na(annotation),"novel", annotation),
         match_miRBase = extract_family(description, "miRBase"),
         match_PmiREN  = extract_family(description, "PmiREN"),
         match_PMRD    = extract_family(description, "PMRD")) %>% 
  rowwise() %>% 
  mutate(identity = ifelse(! match_miRBase == "", match_miRBase, 
                         paste(setdiff(unique(c(match_PmiREN, match_PMRD)),""), collapse = "|"))) %>% 
  mutate(identity = ifelse(identity == "" & !is.na(description), gsub("miRBase_([a-zA-Z]{3})-","",gsub("PmiREN_([a-zA-Z]{3})-","",gsub("PMRD_([a-zA-Z]{3})-","", description))), identity)) %>% 
  mutate(family = paste0(str_extract_all(identity, "miRN?\\-?\\d+") %>% map(unique) %>% map_chr(paste, collapse = "|"),
                         str_extract_all(identity, "-3p|-5p|\\*") %>% map(unique) %>% map_chr(paste, collapse = "|")))

# manual adjustment for miR319 and miR159 that have target
# Picea_mature_215	CTTGGATTGAAGGGAGCTCCA		miR159
# Picea_mature_1586	TTTGGTTTGAAGGGAGCTCTA		miR159
# Picea_mature_3170	CTTGGATTGAAGGGAGCTCCT		miR159
# Picea_mature_4366	 TTGGTTTGAAGGGAGCTCTAC	miR159
# Picea_mature_3859	CTTGGATTGAAGGGAGCTCC		miR159|miR110
# Picea_mature_4734	 TTGGATTGAAGGGAGCTCC		miR159|miR110|miR319 -> is 159
# Picea_mature_3054	 TTGGACTGAAGGGAGCTCCCT	miR319|miR159 -> is 319
# Picea_mature_3428	CTTGGACTGAAGGGAGCTCCC		miR319|miR159 -> is 319
# Picea_mature_387  CTTGGCCTGAAGGGAGCTCC    miR319 -> divergent members of the MIR159/319 superfamily with uncertain sub-family assignment
# Picea_mature_2076	CTTGGACTGAAGGGAGCTCCCT	miR319
# Picea_mature_3033	TTTGGACTGAAGGGAGCTCCC		miR319
# Picea_mature_3679	 TTGGACTGAAGGGAGCTCCC		miR319
# Picea_mature_4860	 TTGGCCTGAAGGGAGCTCC		miR319 -> divergent members of the MIR159/319 superfamily with uncertain sub-family assignment
#   check this.           ^        T=159, C=319
# 
# Picea_mature_2386	AACTGCTGGTTCATGGATTCT	miR159 -> too diif
# Picea_mature_1209	TGCTGCCGACTCATCCATCCG	miR319 -> too diff

# All annotated - manual adjustment too
# Picea_mature_43	  TTTGGCTTGAAGGGAGCTCC	miR159	Plant	variant	PmiREN_Ath-miR159c&PmiREN_Bra-miR159c&PmiREN_Esa-miR159a
# Picea_mature_118	TTTGGATTGAAGAGAGCTCCC	miR159	Plant	known	miRBase_pab-miR159d
# Picea_mature_215	CTTGGATTGAAGGGAGCTCCA	miR159	Plant	known	PMRD_psi-miR159&PmiREN_Ppe-miR159a
# Picea_mature_359	CTTGGATTGAAGAGAGCTCCC	miR159	Plant	variant	PMRD_pta-miR159b&miRBase_pta-miR159b|miRBase_pta-miR159c&PMRD_smo-miR159&PMRD_pta-miR159c&miRBase_smo-miR159
# Picea_mature_986	 TTGGCTTGAAGGGAGCTCC	miR319|miR159	NA	variant	miRBase_pab-miR319j|PmiREN_Fan-miR159d&PmiREN_Fve-miR159c&PmiREN_Fan-miR159a&PmiREN_Fan-miR159b&PmiREN_Fan-miR159c|PmiREN_Gma-miR159h&miRBase_pta-miR159a&PmiREN_Car-miR159c&PMRD_pta-miR159a|PmiREN_Mac-miR159f&PmiREN_Zma-miR159d&miRBase_pab-miR159b&miRBase_cpa-miR159b&PmiREN_Zma-miR159c&PmiREN_Pda-miR159c&PmiREN_Gbi-miR159&miRBase_mdm-miR159b&PMRD_sit-miR110-npr&PmiREN_Pab-miR159b&miRBase_pab-miR159c&miRBase_mdm-miR159a&PmiREN_Mac-miR159c|PmiREN_Ogl-miR159a&PmiREN_Lja-miR159b&PmiREN_Cla-miR159b&PmiREN_Cmo-miR159b&PmiREN_Ogl-miR159c&PmiREN_Osa-miR159c&PmiREN_Oru-miR159c&PmiREN_Osa-miR159e&PmiREN_Osa-miR159d|PmiREN_Cme-miR159b&PMRD_sof-miR319b&PMRD_zma-miR319m|PmiREN_Ath-miR159c&PmiREN_Bra-miR159c&PmiREN_Esa-miR159a
# Picea_mature_1000	 TTGGATTGAAGAGAGCTCC	miR159	NA	variant	PMRD_pta-miR159b&miRBase_pta-miR159b
# Picea_mature_1404	 TTGGATTGAAGGGAGCTCT	miR159	NA	variant	PmiREN_Tae-miR159e&PmiREN_Tae-miR159c&PmiREN_Ogl-miR159b&PmiREN_Tae-miR159d&PmiREN_Osa-miR159a&PmiREN_Oni-miR159b&PmiREN_Oru-miR159b&PmiREN_Bdi-miR159a&PMRD_bdi-miR159&miRBase_bdi-miR159a-3p|PMRD_mtr-miR319b&PmiREN_Fsy-miR159|PmiREN_Ptr-miR159b&PmiREN_Hsy-miR159i&PmiREN_Bol-miR159b&PmiREN_Asp-miR159j&PmiREN_Era-miR159c&PmiREN_Pnu-miR159d&PmiREN_Ptr-miR159a&PmiREN_Vst-miR159a&PmiREN_Aca-miR159i&PmiREN_Bol-miR159a&PmiREN_Ptr-miR159c&PmiREN_Bol-miR159c
# Picea_mature_1586	TTTGGTTTGAAGGGAGCTCTA	miR159	Plant	known	miRBase_pde-miR159
# Picea_mature_2419	TTTGGCTTGAAGGGAGCTCCC	miR159	Plant	variant	PmiREN_Gbi-miR319c|miRBase_pab-miR159d
# Picea_mature_3170	CTTGGATTGAAGGGAGCTCCT	miR159	Plant	known	PMRD_ptc-miR159d&PMRD_zma-miR159c&PMRD_zma-miR159g&PMRD_sof-miR159c&PMRD_sbi-miR159b&miRBase_sbi-miR159b&miRBase_sof-miR159c&miRBase_ptc-miR159d&PMRD_zma-miR159d&PMRD_zma-miR159f&miRBase_zma-miR159d-3p&miRBase_zma-miR159c-3p
# Picea_mature_3559	 TTGGCTTGAAGGGAGCTCCCC	miR319	Plant	variant	PmiREN_Fve-miR319b
# Picea_mature_3760	 TTGGCTTGAAGGGAGCTCCC	miR159	NA	variant	PmiREN_Fan-miR159d&PmiREN_Fve-miR159c&PmiREN_Fan-miR159a&PmiREN_Fan-miR159b&PmiREN_Fan-miR159c
# Picea_mature_3859	CTTGGATTGAAGGGAGCTCC	miR159	Plant	known	PmiREN_Mac-miR159f&PmiREN_Zma-miR159d&miRBase_pab-miR159b&miRBase_cpa-miR159b&PmiREN_Zma-miR159c&PmiREN_Pda-miR159c&PmiREN_Gbi-miR159&miRBase_mdm-miR159b&PMRD_sit-miR110-npr&PmiREN_Pab-miR159b&miRBase_pab-miR159c&miRBase_mdm-miR159a&PmiREN_Mac-miR159c
# Picea_mature_4366	 TTGGTTTGAAGGGAGCTCTAC	miR159	Plant	known	miRBase_pab-miR159a
# Picea_mature_4734	 TTGGATTGAAGGGAGCTCC	miR159	NA	variant	PmiREN_Fan-miR159d&PmiREN_Fve-miR159c&PmiREN_Fan-miR159a&PmiREN_Fan-miR159b&PmiREN_Fan-miR159c|PmiREN_Gma-miR159h&miRBase_pta-miR159a&PmiREN_Car-miR159c&PMRD_pta-miR159a|PmiREN_Mac-miR159f&PmiREN_Zma-miR159d&miRBase_pab-miR159b&miRBase_cpa-miR159b&PmiREN_Zma-miR159c&PmiREN_Pda-miR159c&PmiREN_Gbi-miR159&miRBase_mdm-miR159b&PMRD_sit-miR110-npr&PmiREN_Pab-miR159b&miRBase_pab-miR159c&miRBase_mdm-miR159a&PmiREN_Mac-miR159c|PmiREN_Ogl-miR159a&PmiREN_Lja-miR159b&PmiREN_Cla-miR159b&PmiREN_Cmo-miR159b&PmiREN_Ogl-miR159c&PmiREN_Osa-miR159c&PmiREN_Oru-miR159c&PmiREN_Osa-miR159e&PmiREN_Osa-miR159d|PmiREN_Cme-miR159b&PMRD_sof-miR319b&PMRD_zma-miR319m|PmiREN_Ath-miR159c&PmiREN_Bra-miR159c&PmiREN_Esa-miR159a

# Picea_mature_387	CTTGGCCTGAAGGGAGCTCC	miR319	Plant	known	miRBase_pab-miR319j
# Picea_mature_2066	CTTGGACTGAAGAGAGCTCC	miR319	Plant	variant	PmiREN_Ppa-miR319a&PmiREN_Ppa-miR319b&miRBase_ppt-miR319b&miRBase_ppt-miR319a&PMRD_ppt-miR319a&miRBase_sly-miR319a&PMRD_ppt-miR319b
# Picea_mature_2076	CTTGGACTGAAGGGAGCTCCCT	miR319	Plant	known	PmiREN_Fve-miR319d
# Picea_mature_3033	TTTGGACTGAAGGGAGCTCCC	miR319	Plant	known	PmiREN_Gbi-miR319c
# Picea_mature_3054	 TTGGACTGAAGGGAGCTCCCT	miR319	Plant	known	PmiREN_Ppe-miR319c&miRBase_mtr-miR319d-3p&PMRD_crb-miR159&PmiREN_Spi-miR319a&PmiREN_Mac-miR319a&PmiREN_Aof-miR319d&PMRD_nof-miR159&PmiREN_Bol-miR319a&miRBase_ath-miR319a&PmiREN_Gma-miR319h&PMRD_ath-miR319b&PmiREN_Cti-miR319a&PmiREN_Gma-miR319a&PmiREN_Ach-miR319c&PmiREN_Stu-miR319a&PmiREN_Ghi-miR319c&miRBase_mes-miR319e&PmiREN_Mac-miR319b&PmiREN_Bol-miR319b&PMRD_rco-miR319b&PmiREN_Psa-miR319f&miRBase_vvi-miR319f&PmiREN_Mes-miR319f&PmiREN_Gma-miR319b&PmiREN_Cti-miR319b&miRBase_mes-miR319b&PMRD_ape-miR159&miRBase_cca-miR319&PmiREN_Ach-miR319a&PmiREN_Sly-miR319e&PMRD_aqc-miR319&PmiREN_Bra-miR319d&miRBase_lus-miR319a&PmiREN_Ghi-miR319b&PMRD_cor-miR159&PmiREN_Mac-miR319c&PmiREN_Spi-miR319c&miRBase_bra-miR319-3p&PmiREN_Ach-miR319b&PMRD_vvi-miR319f&PmiREN_Ppe-miR319b&PmiREN_Ach-miR319g&PmiREN_Ghi-miR319a&miRBase_gma-miR319j&PMRD_aly-miR319b&miRBase_rco-miR319b&PMRD_acb-miR159&miRBase_nta-miR319b&PmiREN_Cti-miR319d&miRBase_mes-miR319c&PmiREN_Pap-miR319a&PmiREN_Aof-miR319a&PmiREN_Mac-miR319d&PmiREN_Spi-miR319d&PmiREN_Ach-miR319e&PmiREN_Spe-miR319c&PmiREN_Sly-miR319a&PmiREN_Ghi-miR319e&PMRD_fal-miR159&PmiREN_Can-miR319d&PmiREN_Pap-miR319b&PMRD_aly-miR319a&miRBase_aly-miR319a-3p&PmiREN_Sly-miR319b&PmiREN_Stu-miR319d&miRBase_rco-miR319c&PmiREN_Aof-miR319c&PMRD_ath-miR319a&PMRD_vvi-miR319c&PmiREN_Spe-miR319a&PmiREN_Gma-miR319e&PmiREN_Bra-miR319a&PmiREN_Ach-miR319d&PmiREN_Jcu-miR319b&PmiREN_Mac-miR319e&miRBase_gma-miR319k&PmiREN_Pap-miR319c&PMRD_rco-miR319a&PMRD_svi-miR159&PmiREN_Ghi-miR319d&PmiREN_Spe-miR319b&miRBase_aqc-miR319&PmiREN_Jcu-miR319a&PmiREN_Stu-miR319e&PmiREN_Bra-miR319b&miRBase_aof-miR319a&PmiREN_Bju-miR319c&PMRD_dso-miR159&PmiREN_Can-miR319f&PmiREN_Cme-miR319a&PmiREN_Psa-miR319a&PmiREN_Lsa-miR319b&PmiREN_Stu-miR319f&PmiREN_Cru-miR319b&PMRD_bvl-miR159&miRBase_sly-miR319b&miRBase_tae-miR319&PmiREN_Lsa-miR319a&PmiREN_Cme-miR319b&PmiREN_Psa-miR319g&PmiREN_Cmo-miR319c&PmiREN_Mes-miR319g&PmiREN_Cmo-miR319b&PMRD_mma-miR159&PmiREN_Ghi-miR319f&miRBase_mdm-miR319b-3p&PMRD_ech-miR159&PmiREN_Car-miR319b&PmiREN_Mtr-miR319d&PmiREN_Psa-miR319c&miRBase_aly-miR319b-3p&PmiREN_Mes-miR319c&PmiREN_Cmo-miR319a&miRBase_vvi-miR319c&PmiREN_Cru-miR319c&miRBase_nta-miR319a&PMRD_vvi-miR319b&PmiREN_Ach-miR319f&PmiREN_Lsa-miR319c&PmiREN_Car-miR319h&PMRD_cmi-miR159&PmiREN_Car-miR319a&PmiREN_Ath-miR319b&miRBase_aau-miR319&miRBase_rco-miR319a&PMRD_rco-miR319c&PmiREN_Mtr-miR319c&PmiREN_Mes-miR319d&PMRD_ahi-miR159&miRBase_mes-miR319d&PmiREN_Esa-miR319b&miRBase_hbr-miR319&PmiREN_Ath-miR319a&PMRD_aha-miR159&miRBase_gma-miR319h&PmiREN_Ghi-miR319k&PmiREN_Nnu-miR319e&miRBase_mes-miR319a&PmiREN_Esa-miR319c&miRBase_vvi-miR319b&miRBase_vun-miR319a&miRBase_stu-miR319a-3p&PmiREN_Aly-miR319a&PmiREN_Mes-miR319e&PMRD_tar-miR159&PmiREN_Mtr-miR319b&miRBase_mdm-miR319a&miRBase_ath-miR319b&PMRD_rin-miR159&miRBase_amg-miR319&PmiREN_Cmo-miR319d&PmiREN_Aly-miR319b&miRBase_gma-miR319m
# Picea_mature_3428	CTTGGACTGAAGGGAGCTCCC	miR319|miR159	Plant	known	miRBase_ppt-miR319c&PmiREN_Mac-miR319h&PmiREN_Pvu-miR319d&miRBase_csi-miR159c-3p&PmiREN_Pvu-miR319a&PmiREN_Sha-miR319d&miRBase_ppt-miR319e&PmiREN_Sha-miR319b&PmiREN_Mac-miR319i&PmiREN_Mac-miR319f&PMRD_ppt-miR319d&PmiREN_Can-miR319a&miRBase_ppt-miR319d-3p&PMRD_ppt-miR319e&PmiREN_Pvu-miR319e&PmiREN_Mac-miR319k&PmiREN_Mtr-miR319a&PmiREN_Sha-miR319c&PMRD_ppt-miR319c&PmiREN_Mac-miR319j&PmiREN_Pab-miR319a&PmiREN_Ppa-miR319e&PmiREN_Ppa-miR319c&PmiREN_Gbi-miR319a
# Picea_mature_3548	 TTGGACTGAAGGGAGCTCAC	miR319	Plant	variant	miRBase_pta-miR319&PMRD_pta-miR319
# Picea_mature_3579	 TTGGACTGAAGGGAGCTCCCC	miR319	Plant	known	PmiREN_Fve-miR319b
# Picea_mature_3679	 TTGGACTGAAGGGAGCTCCC	miR319	Plant	known	PmiREN_Bna-miR319e&PmiREN_Ptr-miR319b&PmiREN_Lja-miR319e&miRBase_pab-miR319k&PmiREN_Hsy-miR319n&PMRD_vun-miR319a.1&PmiREN_Vvi-miR319d&PmiREN_Qro-miR319e&PmiREN_Hsy-miR319m&PmiREN_Ptr-miR319a&miRBase_pab-miR319l&PMRD_mtr-miR319&PMRD_ctr-miR319&PmiREN_Bna-miR319d&PmiREN_Hsy-miR319b&PmiREN_Qro-miR319d&miRBase_pab-miR319c&PmiREN_Vun-miR319g&PmiREN_Vun-miR319a&PMRD_ptc-miR319a&PmiREN_Vun-miR319h&PmiREN_Bju-miR319k&PmiREN_Hsy-miR319a&miRBase_atr-miR319c&miRBase_ptc-miR319b&PmiREN_Ptr-miR319c&PmiREN_Gma-miR319j&PmiREN_Bna-miR319c&PmiREN_Vvi-miR319a&PmiREN_Bju-miR319i&PmiREN_Lja-miR319c&miRBase_pab-miR319i&PmiREN_Hsy-miR319s&PmiREN_Can-miR319c&PmiREN_Bju-miR319f&PmiREN_Hsy-miR319t&PmiREN_Qro-miR319c&PmiREN_Vvi-miR319b&PmiREN_Ppe-miR319a&PmiREN_Gma-miR319m&PmiREN_Fex-miR319f&PmiREN_Nta-miR319f&miRBase_mtr-miR319a-3p&PMRD_ptc-miR319b&PmiREN_Fex-miR319i&miRBase_ptc-miR319c&PmiREN_Mac-miR319g&PmiREN_Atr-miR319f&PmiREN_Cit-miR319b&PmiREN_Csi-miR319b&miRBase_gma-miR319b&PmiREN_Vvi-miR319c&miRBase_gma-miR319e&PmiREN_Atr-miR319i&PmiREN_Lja-miR319g&miRBase_pab-miR319b&PmiREN_Bna-miR319b&miRBase_cme-miR319b&PmiREN_Lja-miR319b&PmiREN_Pda-miR319c&miRBase_ppe-miR319a&miRBase_pab-miR319f&miRBase_cme-miR319a&PmiREN_Nta-miR319d&PmiREN_Pab-miR319d&PmiREN_Fex-miR319d&PMRD_ptc-miR319c&PMRD_ltu-miR319&miRBase_pab-miR319a&PmiREN_Hbr-miR319b&PmiREN_Gma-miR319k&PmiREN_Atr-miR319e&miRBase_pab-miR319m&PMRD_gma-miR319a&PmiREN_Hbr-miR319a&PmiREN_Bju-miR319d&miRBase_pab-miR319g&miRBase_atr-miR319a&miRBase_gma-miR319a&PmiREN_Bju-miR319l&PmiREN_Pab-miR319e&PmiREN_Fex-miR319e&PmiREN_Nta-miR319e&PmiREN_Fsy-miR319d&PmiREN_Fan-miR319d&PmiREN_Vun-miR319f&PmiREN_Pda-miR319a&miRBase_lus-miR319b&miRBase_mtr-miR319b-3p&PMRD_ptc-miR319d&PmiREN_Hbr-miR319c&PmiREN_Pda-miR319b&PmiREN_Nnu-miR319d&PmiREN_Bju-miR319e&PmiREN_Nta-miR319l&PmiREN_Hsy-miR319o&PMRD_tae-miR319a&PmiREN_Cla-miR319c&PmiREN_Fsy-miR319c&PmiREN_Pda-miR319e&miRBase_ptc-miR319a&PmiREN_Nta-miR319a&PmiREN_Bju-miR319g&PmiREN_Atr-miR319j&PmiREN_Nta-miR319h&PmiREN_Nta-miR319g&PmiREN_Atr-miR319c&PmiREN_Bju-miR319h&PmiREN_Hbr-miR319d&miRBase_pab-miR319n&PmiREN_Pab-miR319b&PMRD_vun-miR319a.2&miRBase_ctr-miR319&PmiREN_Nnu-miR319c&PmiREN_Nta-miR319b&PmiREN_Pab-miR319c&PmiREN_Fex-miR319c&PmiREN_Peu-miR319&PmiREN_Nnu-miR319b&PmiREN_Pda-miR319d&PmiREN_Nta-miR319c&PmiREN_Bju-miR319j&PmiREN_Sme-miR319a&PmiREN_Cla-miR319a&miRBase_ptc-miR319d&PmiREN_Pvu-miR319f&PmiREN_Mtr-miR319e&PmiREN_Hsy-miR319u&PmiREN_Hsy-miR319r&PmiREN_Hbr-miR319e&PmiREN_Hsy-miR319w&PmiREN_Nnu-miR319a&PmiREN_Nta-miR319j&PmiREN_Atr-miR319a&PmiREN_Sme-miR319b&PmiREN_Ccl-miR319a&PMRD_gma-miR319b
# Picea_mature_4680	ATTGGACTGAAGGGAGCTCCC	miR319|miR159	Plant	variant	PmiREN_Bna-miR319e&PmiREN_Ptr-miR319b&PmiREN_Lja-miR319e&miRBase_pab-miR319k&PmiREN_Hsy-miR319n&PMRD_vun-miR319a.1&PmiREN_Vvi-miR319d&PmiREN_Qro-miR319e&PmiREN_Hsy-miR319m&PmiREN_Ptr-miR319a&miRBase_pab-miR319l&PMRD_mtr-miR319&PMRD_ctr-miR319&PmiREN_Bna-miR319d&PmiREN_Hsy-miR319b&PmiREN_Qro-miR319d&miRBase_pab-miR319c&PmiREN_Vun-miR319g&PmiREN_Vun-miR319a&PMRD_ptc-miR319a&PmiREN_Vun-miR319h&PmiREN_Bju-miR319k&PmiREN_Hsy-miR319a&miRBase_atr-miR319c&miRBase_ptc-miR319b&PmiREN_Ptr-miR319c&PmiREN_Gma-miR319j&PmiREN_Bna-miR319c&PmiREN_Vvi-miR319a&PmiREN_Bju-miR319i&PmiREN_Lja-miR319c&miRBase_pab-miR319i&PmiREN_Hsy-miR319s&PmiREN_Can-miR319c&PmiREN_Bju-miR319f&PmiREN_Hsy-miR319t&PmiREN_Qro-miR319c&PmiREN_Vvi-miR319b&PmiREN_Ppe-miR319a&PmiREN_Gma-miR319m&PmiREN_Fex-miR319f&PmiREN_Nta-miR319f&miRBase_mtr-miR319a-3p&PMRD_ptc-miR319b&PmiREN_Fex-miR319i&miRBase_ptc-miR319c&PmiREN_Mac-miR319g&PmiREN_Atr-miR319f&PmiREN_Cit-miR319b&PmiREN_Csi-miR319b&miRBase_gma-miR319b&PmiREN_Vvi-miR319c&miRBase_gma-miR319e&PmiREN_Atr-miR319i&PmiREN_Lja-miR319g&miRBase_pab-miR319b&PmiREN_Bna-miR319b&miRBase_cme-miR319b&PmiREN_Lja-miR319b&PmiREN_Pda-miR319c&miRBase_ppe-miR319a&miRBase_pab-miR319f&miRBase_cme-miR319a&PmiREN_Nta-miR319d&PmiREN_Pab-miR319d&PmiREN_Fex-miR319d&PMRD_ptc-miR319c&PMRD_ltu-miR319&miRBase_pab-miR319a&PmiREN_Hbr-miR319b&PmiREN_Gma-miR319k&PmiREN_Atr-miR319e&miRBase_pab-miR319m&PMRD_gma-miR319a&PmiREN_Hbr-miR319a&PmiREN_Bju-miR319d&miRBase_pab-miR319g&miRBase_atr-miR319a&miRBase_gma-miR319a&PmiREN_Bju-miR319l&PmiREN_Pab-miR319e&PmiREN_Fex-miR319e&PmiREN_Nta-miR319e&PmiREN_Fsy-miR319d&PmiREN_Fan-miR319d&PmiREN_Vun-miR319f&PmiREN_Pda-miR319a&miRBase_lus-miR319b&miRBase_mtr-miR319b-3p&PMRD_ptc-miR319d&PmiREN_Hbr-miR319c&PmiREN_Pda-miR319b&PmiREN_Nnu-miR319d&PmiREN_Bju-miR319e&PmiREN_Nta-miR319l&PmiREN_Hsy-miR319o&PMRD_tae-miR319a&PmiREN_Cla-miR319c&PmiREN_Fsy-miR319c&PmiREN_Pda-miR319e&miRBase_ptc-miR319a&PmiREN_Nta-miR319a&PmiREN_Bju-miR319g&PmiREN_Atr-miR319j&PmiREN_Nta-miR319h&PmiREN_Nta-miR319g&PmiREN_Atr-miR319c&PmiREN_Bju-miR319h&PmiREN_Hbr-miR319d&miRBase_pab-miR319n&PmiREN_Pab-miR319b&PMRD_vun-miR319a.2&miRBase_ctr-miR319&PmiREN_Nnu-miR319c&PmiREN_Nta-miR319b&PmiREN_Pab-miR319c&PmiREN_Fex-miR319c&PmiREN_Peu-miR319&PmiREN_Nnu-miR319b&PmiREN_Pda-miR319d&PmiREN_Nta-miR319c&PmiREN_Bju-miR319j&PmiREN_Sme-miR319a&PmiREN_Cla-miR319a&miRBase_ptc-miR319d&PmiREN_Pvu-miR319f&PmiREN_Mtr-miR319e&PmiREN_Hsy-miR319u&PmiREN_Hsy-miR319r&PmiREN_Hbr-miR319e&PmiREN_Hsy-miR319w&PmiREN_Nnu-miR319a&PmiREN_Nta-miR319j&PmiREN_Atr-miR319a&PmiREN_Sme-miR319b&PmiREN_Ccl-miR319a&PMRD_gma-miR319b|PmiREN_Hsy-miR319q&PmiREN_Rco-miR319a&PmiREN_Stu-miR319b&miRBase_mes-miR319f&miRBase_cpa-miR319&PmiREN_Ghi-miR319i&PmiREN_Hsy-miR319p&PmiREN_Sme-miR319d|miRBase_fve-miR159c
# Picea_mature_4860	 TTGGCCTGAAGGGAGCTCC	miR319	NA	variant	miRBase_pab-miR319j

# Picea_mature_918	GAGCTTCCTTAAGTCCAAAG	miR159	NA	variant	miRBase_vca-miR159-5p
# Picea_mature_1209	TGCTGCCGACTCATCCATCCG	miR319	NA	variant	PMRD_smo-miR319&miRBase_smo-miR319
# Picea_mature_2386	AACTGCTGGTTCATGGATTCT	miR159	Plant	known	miRBase_pab-miR159i&PmiREN_Pab-miR159d&miRBase_pab-miR159g&PmiREN_Pab-miR159c&miRBase_pab-miR159h&PmiREN_Pab-miR159f
# Picea_mature_3694	GAGTTCTCTTCAGTCCACTC	miR319	NA	variant	miRBase_zma-miR319a-5p&miRBase_zma-miR319c-5p

# test <- matureDF %>% 
#   filter(id %in% c(
#     "Picea_mature_918",
#     "Picea_mature_1209",
#     "Picea_mature_2386",
#     "Picea_mature_3694"
#   ))

# manual edit
matureDF <- matureDF %>% 
  mutate(family = case_when(
    id %in% c("Picea_mature_986","Picea_mature_3559") ~ "miR159",
    id %in% c("Picea_mature_3428","Picea_mature_4680") ~ "miR319",
    id == "Picea_mature_1209" ~ "",
    id == "Picea_mature_2386" ~ "miR159-phased",
    .default = family)) %>% 
  mutate(identity = case_when(
    id == "Picea_mature_986" ~ gsub("pab-miR319j\\|", "", identity),
    id == "Picea_mature_3559" ~ "miR159",
    id == "Picea_mature_3428" ~ gsub("miR159-3p\\|", "", identity),
    id == "Picea_mature_4680" ~ gsub("\\|miR159", "", identity),
    id == "Picea_mature_1209" ~ "",
    .default = identity))

# sanity check
# Compute pairwise %identity for all sequence pairs within a family
family_identity <- function(seqs) {
  if (length(seqs) < 2) return(NA_real_)  # singleton family, nothing to compare
  combn(seqs, 2, function(pair) {
    aln <- pairwiseAlignment(pair[1], pair[2], type = "global-local")
    pid(aln, type = "PID2")   # % identity over aligned (non-gap) positions
  })
}

family_check <- matureDF %>%
  filter(family != "") %>%
  group_by(family) %>%
  summarise(n_members = n(), seqs = list(sequence), .groups = "drop") %>%
  mutate(
    pid_values = map(seqs, family_identity),
    min_pid    = map_dbl(pid_values, ~ ifelse(all(is.na(.x)), NA, min(.x))),
    mean_pid   = map_dbl(pid_values, ~ ifelse(all(is.na(.x)), NA, mean(.x)))
  ) %>%
  select(family, n_members, min_pid, mean_pid) %>%
  arrange(min_pid)
saveRDS(family_check, here("data/analysis/results/family_check.rds"))

# Axtell and Meyers (2018) figure 3, all, with Claude's comment
conservedFam <- family_check %>% 
  filter(family %in% c(
    "miR156", # core deeply-conserved in land plant
    "miR159", # core deeply-conserved in land plant
    "miR160", # core deeply-conserved in land plant
    "miR162",
    "miR164",
    "miR166", # core deeply-conserved in land plant
    "miR167", # core deeply-conserved in land plant
    "miR168", # core deeply-conserved in land plant
    "miR169", # core deeply-conserved in land plant
    "miR171", # core deeply-conserved in land plant
    "miR172", # core deeply-conserved in land plant
    "miR319", # core deeply-conserved in land plant
    "miR390", # core deeply-conserved in land plant
    "miR393", # core deeply-conserved in land plant
    "miR394", # core deeply-conserved in land plant
    "miR395", # core deeply-conserved in land plant
    "miR396", # core deeply-conserved in land plant
    "miR397", # core deeply-conserved in land plant
    "miR398", # core deeply-conserved in land plant
    "miR399", # core deeply-conserved in land plant
    "miR403",
    "miR408", # core deeply-conserved in land plant
    "miR477",
    "miR479",
    "miR482",
    "miR529",
    "miR530",
    "miR535",
    "miR827",
    "miR1446",
    "miR1863",
    "miR2111",
    "miR2118",
    "miR2275",
    "miR2950",
    "miR4376"))

# check fam than min pid < 80
# miR399 40.0
# miR482 42.8
# miR408 53.3
# miR169 71.4
# miR395 71.4
# miR529 71.4
# miR164 75.0
# miR159 76.2
test <- matureDF %>% 
  filter(family %in% c(
    "miR399",
    "miR482",
    "miR408",
    "miR169",
    "miR395",
    "miR529",
    "miR164",
    "miR159"
  )) %>% 
  split(., .$family)

# $miR159
# "TTTGGCTTGAAGGGAGCTCC"  
# "TTTGGATTGAAGAGAGCTCCC" 
# "CTTGGATTGAAGGGAGCTCCA" 
# "CTTGGATTGAAGAGAGCTCCC"
#  "TTGGCTTGAAGGGAGCTCC"   
#  "TTGGATTGAAGAGAGCTCC"   
# "TTTGGTTTGAAGGGAGCTCTA" 
# "TTTGGCTTGAAGGGAGCTCCC"
# "CTTGGATTGAAGGGAGCTCCT" 
#  "TTGGCTTGAAGGGAGCTCCCC" 
#  "TTGGCTTGAAGGGAGCTCCC"  
# "CTTGGATTGAAGGGAGCTCC" 
#  "TTGGTTTGAAGGGAGCTCTAC" 
#  "TTGGATTGAAGGGAGCTCC"  
# 
# $miR164 ?????
#  "GGAGAAGCAGGGCACGTGCG"  
# "TGGAGAAGCAGGGCACGTGCG" 
# "CATGTGCTCTTTTTCTCCATC" 
# "CACGTGCTCCCCTTCTCCAA" 
# "CACGTGCTCCCCTTCTCCAAC"
# 
# $miR169 the TGATTTGCC ones are variants, TGACTTGCC are mostly known
# "CTCAGCCAAGAATGACTTGCCCG" 
#   "CAGCCAAGAATGACTTGCCCG"   
#  "TCAGCCAAGAATGACTTGCCC"   
#   "CAGCCAAGGATGACTTGCCGGC" 
#   "CAGCCAAGGATGACTTGCCTA"   
#   "GAGCCAAGGATGATTTGCCGA"  
#   "GAGACAAGGATGATTTGCCGA"   
#   "CAGCCAAGAATGATTTGCCCGC"  
#    "AGCCAAGGATGATTTGCCGA"   
# 
# $miR395 the ones that have A and T in the middle are both known and variants from miRBase
# "GCTGAAGTGTTTGGGGGAATTC" 
#  "ATGAAGTGTTTGGAGGAACTT"  
#  "CAGAAGTGTTTGGGGGGACTC"  
#  "CAGAAGTGTTTGGGGGGACT"  
#  "CTGAAGTGTTTGGAGGAACTT"  
#  "CAGAAGTGTTTGGGTGGACTC"  
#  "CTGAAGTGTTTGGGGGAACTC" 
# 
# $miR399 ????
# "GGCAGCTCTCCTTTGGCAGG"  
# "CGCCAAAGGAGAGTTGCTCAG" 
# "TGCCAAAGGAGAAGTGCCCTG" 
#  "GCCAAAGGAGAGTTGCCTT"  
# "CGCCAAAGGAGAGTTGCCCTG" 
# "CGCCAAAGGAGAGTTGCGCTG" 
# "CGCCAAAGGAGAGTTGCCT"   
# "TGCCAAAGGAGAGTTGCTCTG"
# " GCCAAAGGAGAGTTGCCTTC"  
# "CGCCAAAGGAGAGTTGCCTT"  
# "TGCCAAAGGAGAATTGCCCTG"
# 
# $miR408 ????
# "TGCACTGCCTCTTCCCTGGCT" 
# "ACAGGGAAGAGTTAGGGCATG"
# 
# $miR482 a diverse family targeting NBS-LRR disease-resistance genes
# "TCTTTCCACTTCTACCCATTTC" 
# "TCTTTCCTACTCCTCCCAAGCC" 
# "TCTTCCCTATTCCTCCCATTCC" 
# "CTTTCCACTTCTACCCATTTCC"
# "TTTTTCCCACTCCTCCCAAGCC" 
# "TTCCCTATTCCTCCCATTCCT"  
# "TTGAGAAACTGTGAGCCAAATC" 
# "TTTTTCCTACACCTCCCAACCC"
# "TTTCCTACTCCTCCCAAGCCC"  
# "TTCCCTATTCCTCCCATTCC"   
# "TTTCCAACGCCTCCCATGCCT"  
# "CTTTCCAACGCCGGCCATTCCA"
# "TTTTCCCGGCTCCGCCCATTCC" 
# "AATGGGTAGATTGGAATGACC"  
# "TTTTTCCTACACCTCCCAACC"  
# "TCTTTTCAGTTCTACCCATTTC"
# "CTTCCCCACCCCTCCCATT"    
# "TTGAGAGACTGTGAGCCAAATC" 
# "TCTTTCCAACGCCGGCCATTCC" 
# "TCTTCCCAACCCCTCCCATGCC"
# "TTCCCTATTCCTCCCATTCCTA" 
# "TCTTCCCTAAACCTCCCAAACC" 
# "TTTCCAACGCCTCCCATGCCTA" 
# "TGAGAGACTGTGAGCCAAATCA"
# "TTTCCTACTCCTCCCAAGCCCA" 
# "TTTTTCCAACTCCACCTAGACC" 
# "TGAGAAACTGTGAGCCAAATCA" 
# "CTTTCCACTTCTACCCATTTC" 
# "TCTTTCCTACTCCTCCCATTCC" 
# "TCTTCCCTACTCCTCCCATTCC"
# 
# $miR529 the two GAGATA are known and variants from PmiREN, the rest GAGAGA are known and variants from either miRBase or PmiREN
# "TAGAAGAGAGATAGTACAACC" 
#  "AGAAGAGAGATAGTACAACCA" 
# "CAGAAGAGAGAGAGCACAGCC" 
#  "CGAAGAGAGAGAGCACCTCCT"
#  "TGAAGAGAGAGAGCACAACCC" 
#  "AGAAGAGAGAGAGCACAGCCA" 
#  "TGAAGAGAGAGAGCACAACC"  
#  "AGAAGAGAGATAGCACAACCT"
#  "AGAAGAGAGAGAGTACAGCCC" 
#  "AGAAGAGAGAGAGCACAGCCT" 
#  "TGAAGAGAGAGAGCACAGCC"  
#  "AGAAGAGAGAGAGCACCGCCT"
#   "GAAGAGAGAGAGCACAGCC"   
#  "TGAAGAGAGAGAGCACCGCCT" 
#  "CGAAGAGAGAGAGCACCGCCT"

# check fam than min pid 80-100
# miR156 80.0
# miR398 80.9
# miR393 83.3
# miR396 83.3
# miR535 85.7
# miR160 89.5
# miR319 89.5
# miR171 90.0
# miR2111 90.0
# miR397 90.0
# miR162 90.5
# miR166 94.7
# miR172 94.7
# miR394 94.7
# miR167 100
# miR168 100
# miR2950 100

test <- matureDF %>% 
  filter(family %in% c(
    "miR156",
    "miR398",
    "miR393",
    "miR396",
    "miR535",
    "miR160",
    "miR319",
    "miR171",
    "miR2111",
    "miR397",
    "miR162",
    "miR166",
    "miR172",
    "miR394",
    "miR167",
    "miR168",
    "miR2950"
  )) %>% 
  split(., .$family)

# $miR156 both GAG and TAG are known ana variants in miRBase
#   "TGACAGGAGATAGAGAGCAC"   
#   "TGACAGAAGAGAGTGGGCAC"   
#   "TGACAGAAGAGAGTGAGCAC"   
#   "TGACAGAAGAGAGAGAGCAC"  
# "GCTGACAGGAGATAGAGAGCAC" 
#   "TGACAGAAGATAGTGGGCAC"   
#   CTGACAGGAGATAGAGAGCAC" 
# 
# $miR160 ???? T in the middle
#   "TGCCTGGCTCTCTGCATGCCA" 
# "TCTGCCTGGCTCCCTGTATGC" 
#  "CTGCCTGGCTCCCTGTATGCC" 
#  "GTGCCTGGCTCCCTGTATGCC"
# "TGTGCCTGGCTCCCTGTATGC" 
#   "TGCCTGGCTCCCTGTATGCC" 
# 
# $miR162
# "TCGATAAGCCTCTGCATCCAG" 
# "TCGATAAACCTCTGCATCCGG"
# 
# $miR166
# "TCGGACCAGGCTTCATTCCCC" 
# "CCGGACCAGGCTTCATTCCCC" 
# "TCGGACCAGGCTTCATTCC"  
# 
# $miR167
# "TGAAGCTGCCAGCATGATCTG"  
# "TGAAGCTGCCAGCATGATCTGG"
# 
# $miR168
# "TCGATTGGTGCAGATCGGG"   
# "TCGATTGGTGCAGATCGGGA"  
# "TCGATTGGTGCAGATCGGGAA"
# 
# $miR171 ??? seed position slided
# "TGATTGAGCCGTGCCAATATC"  
#      "GAGCCGTGCCAATATCACAT"   
#    "TTGAGCCGTGCCAATATCGCAC" 
#     "TGAGCCGTGCCAATATCGCA"   
#    "TTGAGCCGTGCCAATATCGC"  
#    "TTGAGCCGTGCCAATATCGCA"  
#     "TGAGCCGTGCCAATATCACAT" 
# 
# $miR172
#  "GAATCCTGATGATGCTGCATC" 
#   "AATCCTGATGATGCTGCAT"   
# "AGAATCCTGATGATGCTGCAT" 
# "AGAATCTTGATGATGCTGCAT"
# 
# $miR2111 ????
# "AACCCTAGAATGCAGATTATC" 
# "TAATCTGCATCCTGGGGTTTG"
# 
# $miR2950
# "TTCCATCTCTTGCACACTGAG" 
# "TTCCATCTCTTGCACACTGA" 
# 
# $miR319
# "CTTGGCCTGAAGGGAGCTCC"   
# "CTTGGACTGAAGAGAGCTCC"   
# "CTTGGACTGAAGGGAGCTCCCT" 
# "TTTGGACTGAAGGGAGCTCCC" 
#  "TTGGACTGAAGGGAGCTCCCT"  
# "CTTGGACTGAAGGGAGCTCCC"  
#  "TTGGACTGAAGGGAGCTCAC"   
#  "TTGGACTGAAGGGAGCTCCCC" 
#  "TTGGACTGAAGGGAGCTCCC"   
# "ATTGGACTGAAGGGAGCTCCC"  
#  "TTGGCCTGAAGGGAGCTCC"   
# 
# $miR393 ???
# "ATCATGCAATCCCTTCGGAAT"   
# "TCCAAAGGGATTGCATTGATTC"  
# "TCCAAAGGGATTGCATTGATT"  
# "TCCAAAGGGATTGCATTGATTCA"
# 
# $miR394
#  "TGGCATTCCGTCCACCTCCA"  
# "CTGGCATTCTGTCCACCTCCA" 
# "CTGGCATTCCGTCCACCTCC"  
# " TGGCATTCTGTCCACCTCCA" 
# "CTGGCATTCCGTCCACCTCCA"
# 
# $miR396
#   "TCCCACTGCTTTCTTGAACTT"  
#   "TCCAGCAGCTTTCTTGAACT"   
#   "TCCAGCAGCTTTCTTGAACTT"  
# "TTTCCACATGCTTTCTTGAACT"
#   "TCCACATGCTTTCTTGAACT"   
#   "TTCCACAGCTTTCTTGAACTA"  
#   "TCCCACTGCTTTCTTGAACT"   
#  "TTCCACATGCTTTCTTGAACTT"
#    "CCAGCAGCTTTCTTGAACT"    
#    "CCCACTGCTTTCTTGAACT"    
#   "TCCAGCAGCTTTCTTGAAC"    
#  "TTCCCACTGCTTTCTTGAACT" 
#  "TTCCACATGCTTTCTTGAACT"  
#   "TCCCACTGCTTTCTTGAAC"   
# 
# $miR397
# "TATTGAGTGCAGCGTTGAC"   
#  "ATTGAGTGCAGCGTTGACGA"  
# "TATTGAGTGCAGCGTTCACGA" 
# "CATTGAGTGCAGCGTTGACGA"
# "TATTGAGTGCAGCGTTGACGA" 
# "TATTGAGTGCAGCGTTCACG"  
# "TATTGAGTGCAGCGTTGACGG" 
# "TATTGAGTGTAGCGTTGACGA"
# "TATTGAGTGCAGCGTTGACG" 
# 
# $miR398
# "CGTGTTCTCAGGTCACCCCTG" 
# "CGTGTTCCCAGGTCGCCTCAG" 
# "CGTGTTCCCAGGTCGCCCCTG"
# 
# $miR535
# "TGACAATGAGAGAGAGCACGCTT" 
# "TGACAATGAGAGAGAGCACGCT"  
# "TGACAACGAGAGAGAGCACGCT" 
# "TGACAGAGAGAGAGAGCACGC"   
# "TGACAACTAGAGAGAGCACGC"   
# "TGACAACGAGAGAGAGCACGC"  

# Quastionable
test <- matureDF %>% 
  filter(family %in% c(
    "miR164",
    "miR399",
    "miR408",
    "miR160",
    "miR171",
    "miR2111",
    "miR393"
  )) %>% 
  split(., .$family)

# $miR164 ???
#  "GGAGAAGCAGGGCACGTGCG"  Picea_mature_1230  variant  match to many species miR164
# "TGGAGAAGCAGGGCACGTGCG"  Picea_mature_1595  known    match to many species miR164
# "CATGTGCTCTTTTTCTCCATC"  Picea_mature_1251  known    miRBase_pab-miR164b
# "CACGTGCTCCCCTTCTCCAA"   Picea_mature_2626  variant  miRBase_pab-miR164a&miRBase_fve-miR164a-3p
# "CACGTGCTCCCCTTCTCCAAC"  Picea_mature_4703  known    miRBase_pab-miR164a&miRBase_fve-miR164a-3p
# 
# $miR399 ???
# "GGCAGCTCTCCTTTGGCAGG"  Picea_mature_923  known  miRBase_pab-miR399f&miRBase_pab-miR399g
# "CGCCAAAGGAGAGTTGCTCAG" 
# "TGCCAAAGGAGAAGTGCCCTG" 
#  "GCCAAAGGAGAGTTGCCTT"  
# "CGCCAAAGGAGAGTTGCCCTG" 
# "CGCCAAAGGAGAGTTGCGCTG" 
# "CGCCAAAGGAGAGTTGCCT"   
# "TGCCAAAGGAGAGTTGCTCTG"
# " GCCAAAGGAGAGTTGCCTTC"  
# "CGCCAAAGGAGAGTTGCCTT"  
# "TGCCAAAGGAGAATTGCCCTG"
# 
# $miR408 ???
# "TGCACTGCCTCTTCCCTGGCT"  Picea_mature_2303  known  match to many species miR408
# "ACAGGGAAGAGTTAGGGCATG"  Picea_mature_4138  known  miRBase_pab-miR408
# 
# $miR160 ??? T in the middle
#   "TGCCTGGCTCTCTGCATGCCA"  Picea_mature_1694  variant(known seq has CCC, not CTC)  miREN_Zos-miR160a&PmiREN_Nta-miR160d&miRBase_ptc-miR160h&PMRD_ptc-miR160h&PmiREN_Gma-miR160h&miRBase_nta-miR160d&PmiREN_Mac-miR160a
# "TCTGCCTGGCTCCCTGTATGC" 
#  "CTGCCTGGCTCCCTGTATGCC" 
#  "GTGCCTGGCTCCCTGTATGCC"
# "TGTGCCTGGCTCCCTGTATGC" 
#   "TGCCTGGCTCCCTGTATGCC" 
# 
# $miR171 ??? seed position slided
# "TGATTGAGCCGTGCCAATATC"      Picea_mature_2839  known    match to many species miR171
#      "GAGCCGTGCCAATATCACAT"  Picea_mature_782   variant  PmiREN_Smo-miR171e&PMRD_smo-miR171b&PmiREN_Smo-miR171b&miRBase_smo-miR171b
#    "TTGAGCCGTGCCAATATCGCAC"  Picea_mature_1545  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#     "TGAGCCGTGCCAATATCGCA"   Picea_mature_1643  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#    "TTGAGCCGTGCCAATATCGC"    Picea_mature_1763  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#    "TTGAGCCGTGCCAATATCGCA"   Picea_mature_3003  known    PmiREN_Pab-miR171c&miRBase_pab-miR171a
#     "TGAGCCGTGCCAATATCACAT"  Picea_mature_4024  known    PmiREN_Smo-miR171e&PMRD_smo-miR171b&PmiREN_Smo-miR171b&miRBase_smo-miR171b
# 
# $miR2111 ???
# "AACCCTAGAATGCAGATTATC"  Picea_mature_1142  variant  PmiREN_Mtr-miR2111y
# "TAATCTGCATCCTGGGGTTTG"  Picea_mature_3134  known    miRBase_pab-miR2111&PmiREN_Pab-miR2111a
# 
# $miR393 ???
# "ATCATGCAATCCCTTCGGAAT"   Picea_mature_290  variant  miRBase_sly-miR393&miRBase_stu-miR393-3p
# "TCCAAAGGGATTGCATTGATTC"  
# "TCCAAAGGGATTGCATTGATT"  
# "TCCAAAGGGATTGCATTGATTCA"

# Redo family assignment with both 5p/3p and * nomenclature
extract_family <- function(description, db_prefix) {
  pattern <- paste0(db_prefix, "_([a-zA-Z]{3})-(miRN?\\d+)([a-z]?)(-3p|-5p|\\*)?")
  
  matches <- str_match_all(description, pattern)
  
  map(matches, function(m) {
    if (nrow(m) == 0) return(tibble(family = character(0), identity = character(0), arm = character(0)))
    
    parsed <- as.data.frame(m, stringsAsFactors = FALSE) %>%
      setNames(c("full", "species", "num", "locus", "arm")) %>%
      mutate(
        locus  = ifelse(is.na(locus), "", locus),
        arm    = ifelse(is.na(arm), "", arm),
        is_pab = species %in% c("Pab", "pab"),
        style  = case_when(
          arm %in% c("-3p", "-5p") ~ "explicit",
          arm == "*"               ~ "star",
          TRUE                     ~ "none"
        )
      )
    
    # collapse non-pab duplicates when a pab entry shares (num, arm)
    parsed <- parsed %>%
      group_by(num, arm) %>%
      filter(!(any(is_pab) & !is_pab)) %>%
      ungroup()
    
    # for each locus, record arm information
    parsed <- parsed %>%
      group_by(num, locus) %>%
      mutate(has_explicit = any(style == "explicit")) %>%
      ungroup() %>%
      # priority: keep the explicit -3p/-5p entry when one exists for this locus;
      # only fall back to bare/star if no explicit version was found
      filter(!(has_explicit & style != "explicit")) %>%
      mutate(
        family   = num,
        identity = ifelse(is_pab,
                          paste0(species, "-", num, locus, arm),
                          paste0(num, arm)),
        arm      = gsub("^-", "", arm)   # "-3p" -> "3p", "-5p" -> "5p", "*" stays "*", "" stays ""
      ) %>%
      distinct(identity, .keep_all = TRUE) %>%
      select(family, identity, arm)
    
    parsed
  })
}

matureDF <- as.data.frame(matureSeq) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::rename(sequence = V1) %>% 
  mutate(P_criteria = ifelse(id %in% mature_plant, "Plant", NA)) %>% 
  left_join(annot, by = "id") %>% 
  mutate(
    parsed_miRBase = extract_family(description, "miRBase"),
    family_miRBase   = map_chr(parsed_miRBase, ~ paste(unique(.x$family),   collapse = "|")),
    identity_miRBase = map_chr(parsed_miRBase, ~ paste(unique(.x$identity), collapse = "|")),
    arm_miRBase      = map_chr(parsed_miRBase, ~ paste(unique(setdiff(.x$arm,"")),      collapse = "|")),
    
    parsed_PmiREN = extract_family(description, "PmiREN"),
    family_PmiREN   = map_chr(parsed_PmiREN, ~ paste(unique(.x$family),   collapse = "|")),
    identity_PmiREN = map_chr(parsed_PmiREN, ~ paste(unique(.x$identity), collapse = "|")),
    arm_PmiREN      = map_chr(parsed_PmiREN, ~ paste(unique(setdiff(.x$arm,"")),      collapse = "|")),
    
    parsed_PMRD = extract_family(description, "PMRD"),
    family_PMRD   = map_chr(parsed_PMRD, ~ paste(unique(.x$family),   collapse = "|")),
    identity_PMRD = map_chr(parsed_PMRD, ~ paste(unique(.x$identity), collapse = "|")),
    arm_PMRD      = map_chr(parsed_PMRD, ~ paste(unique(setdiff(.x$arm,"")),      collapse = "|"))
  ) %>%
  select(-starts_with("parsed_")) %>%
  rowwise() %>%
  mutate(
    family   = ifelse(family_miRBase != "", family_miRBase,
                      paste(setdiff(unique(c(family_PmiREN, family_PMRD)), ""), collapse = "|")),
    identity = ifelse(identity_miRBase != "", identity_miRBase,
                      paste(setdiff(unique(c(identity_PmiREN, identity_PMRD)), ""), collapse = "|")),
    arm      = ifelse(arm_miRBase != "", arm_miRBase,
                      paste(setdiff(unique(c(arm_PmiREN, arm_PMRD)), ""), collapse = "|"))
  )

# manual adjustment for miR319 and miR159
matureDF <- matureDF %>% 
  mutate(family = case_when(
    id %in% c("Picea_mature_986","Picea_mature_3559") ~ "miR159",
    id %in% c("Picea_mature_3428","Picea_mature_4680") ~ "miR319",
    id == "Picea_mature_1209" ~ "",
    .default = family)) %>% 
  mutate(identity = case_when(
    id == "Picea_mature_986" ~ gsub("pab-miR319j\\|", "", identity),
    id == "Picea_mature_3559" ~ "miR159",
    id == "Picea_mature_3428" ~ gsub("miR159-3p\\|", "", identity),
    id == "Picea_mature_4680" ~ gsub("\\|miR159", "", identity),
    id == "Picea_mature_1209" ~ "",
    .default = identity)) %>% 
  mutate(arm = ifelse(id == "Picea_mature_2386", "phased", arm))

# sanity check
# --- 0. Reshape matureDF so each row is one (family, identity, arm) membership ---
# family/identity/arm were built in lockstep (same source, same order), so
# separate_rows can split them together safely.
sanity_df <- matureDF %>%
  filter(!family %in% c("","NA")) %>%
  separate_rows(family, arm, sep = "\\|")

# --- 1. Pairwise stats helper, with coverage on BOTH sides ---
# Switched to type = "overlap": it allows free end-gaps on both sequences
# (no penalty for unaligned prefix/suffix on either side), so coverage is
# symmetric regardless of which sequence is passed first — "global-local"
# would force pair[1] to align in full, biasing cov1 toward ~100% by
# construction and making the two sides not comparable.
pairwise_stats <- function(seqs, ids, min_coverage = 0.5) {
  n <- length(seqs)
  if (n < 2) return(tibble(id1 = character(0), id2 = character(0), pid = numeric(0),
                           cov1 = numeric(0), cov2 = numeric(0), pass_coverage = logical(0)))
  pairs <- combn(n, 2)
  map_dfr(seq_len(ncol(pairs)), function(i) {
    a <- pairs[1, i]; b <- pairs[2, i]
    s1 <- seqs[a]; s2 <- seqs[b]
    aln <- pairwiseAlignment(s1, s2, type = "overlap")
    aligned_len <- nmatch(aln) + nmismatch(aln)   # non-gap aligned columns
    cov1 <- aligned_len / nchar(s1)
    cov2 <- aligned_len / nchar(s2)
    tibble(
      id1 = ids[a], id2 = ids[b],
      pid  = pid(aln, type = "PID2"),
      cov1 = cov1, cov2 = cov2,
      pass_coverage = cov1 > min_coverage & cov2 > min_coverage
    )
  })
}

# --- 2. Within-arm check: each (family, arm) group compared separately ---
family_check <- sanity_df %>%
  group_by(family, arm) %>%
  summarise(n_members = n(), seqs = list(sequence), ids = list(id), .groups = "drop") %>%
  mutate(
    pairwise       = map2(seqs, ids, pairwise_stats),
    valid_pairwise = map(pairwise, ~ filter(.x, pass_coverage)),
    n_total_pairs  = map_int(pairwise, nrow),
    n_valid_pairs  = map_int(valid_pairwise, nrow),
    low_cov_pairs  = n_total_pairs - n_valid_pairs,     # pairs excluded for poor overlap
    min_pid        = map_dbl(valid_pairwise, ~ ifelse(nrow(.x) == 0, NA_real_, min(.x$pid))),
    mean_pid       = map_dbl(valid_pairwise, ~ ifelse(nrow(.x) == 0, NA_real_, mean(.x$pid)))
  ) %>%
  select(family, arm, n_members, n_total_pairs, n_valid_pairs, low_cov_pairs, min_pid, mean_pid) %>%
  arrange(min_pid)
saveRDS(family_check, here("data/analysis/results/family_check.rds"))

# Families/arms worth a manual look:
family_check %>% filter(n_members > 1, min_pid < 70 | is.na(min_pid))
test <- matureDF %>% 
  filter(family %in% c("miR11408", "miR11487", "miR2111"))
# they are likely different arms, but no arm info was indicated in database
# if they turn up in target step, we will try identifying them
# but if not, we will keep them in supp table like this

# --- 3. Cross-arm check: compare sequences with DIFFERENT arm labels ---
# within the same family (e.g. "" vs "5p", "*" vs "3p", or "5p" vs "3p").
cross_arm_check <- sanity_df %>%
  group_by(family) %>%
  filter(n_distinct(arm) > 1) %>%          # only families with >1 arm label present
  summarise(seqs = list(sequence), ids = list(id), arms = list(arm), .groups = "drop") %>%
  mutate(
    pairwise = pmap(list(seqs, ids, arms), function(s, i, a) {
      arm_lookup <- setNames(a, i)
      pairwise_stats(s, i) %>%
        mutate(arm1 = arm_lookup[id1], arm2 = arm_lookup[id2]) %>%
        filter(arm1 != arm2)               # keep only cross-labeled pairs
    })
  ) %>%
  select(family, pairwise) %>%
  unnest(pairwise)

# Flag: high identity + good coverage but conflicting arm labels
# -> either "no-suffix"/"*" corresponds to a known arm, or a possible
#    genuine mislabeling (e.g. same sequence called both 5p and 3p).
cross_arm_flags <- cross_arm_check %>%
  filter(pass_coverage, pid > 70) %>%
  arrange(desc(pid))

cross_arm_flags
# only miR950 has 3p/5p conflict with pid > 70
test <- matureDF %>% 
  filter(family == "miR950")
# Picea_mature_1259 TGACATCTGGGCCACGATGGTT ""   
# Picea_mature_1264 TCACGTCTGGGCCGCGATGGTT ""   
# Picea_mature_2719 TAACATCTGGGCCACGATGGTT ""   
# Picea_mature_3476 TTACGTCTGGGCCGCGATGGTT ""   
# Picea_mature_4074 TGACGTCTGGTCCACGATGGTT ""   
# Picea_mature_4162 TCACATCTGGGCCACGATG    "5p" 
# Picea_mature_4526 TCACGTCAGGGCCACGATGGTT ""   
# Picea_mature_4655 TGACGTCTGGTCCACGATG    ""   
# Picea_mature_4780 TCACATCTGGGCCACGATGGTT "5p" 
# Picea_mature_2103 AGAACCACCGTGGCCCTGACG  ""   
# Picea_mature_1015 AGAACCACCGTGGCCCTGATG  ""   
# Picea_mature_1005 TCTGGGCCCCGGTGGTTTATGA "3p" 
# reverse comp 1005 TCATAAACCACCGGGGCCCAGA
# they just have many G and C, but the 5p and 3p are not similar

# Axtell and Meyers (2018) figure 3, all, with Claude's comment
conservedFam <- family_check %>% 
  filter(family %in% c(
    "miR156", # core deeply-conserved in land plant
    "miR159", # core deeply-conserved in land plant
    "miR160", # core deeply-conserved in land plant
    "miR162",
    "miR164",
    "miR166", # core deeply-conserved in land plant
    "miR167", # core deeply-conserved in land plant
    "miR168", # core deeply-conserved in land plant
    "miR169", # core deeply-conserved in land plant
    "miR171", # core deeply-conserved in land plant
    "miR172", # core deeply-conserved in land plant
    "miR319", # core deeply-conserved in land plant
    "miR390", # core deeply-conserved in land plant
    "miR393", # core deeply-conserved in land plant
    "miR394", # core deeply-conserved in land plant
    "miR395", # core deeply-conserved in land plant
    "miR396", # core deeply-conserved in land plant
    "miR397", # core deeply-conserved in land plant
    "miR398", # core deeply-conserved in land plant
    "miR399", # core deeply-conserved in land plant
    "miR403",
    "miR408", # core deeply-conserved in land plant
    "miR477",
    "miR479",
    "miR482",
    "miR529",
    "miR530",
    "miR535",
    "miR827",
    "miR1446",
    "miR1863",
    "miR2111",
    "miR2118",
    "miR2275",
    "miR2950",
    "miR4376"))

# lets check the ones we have problem with before 
filter(conservedFam, family %in% c("miR164", "miR399", "miR408", "miR160", "miR171", "miR2111", "miR393")) %>% arrange(family)

# $miR164 good!
#  "GGAGAAGCAGGGCACGTGCG"  Picea_mature_1230  variant  match to many species miR164 -> 5p in Arabidopsis and bra
# "TGGAGAAGCAGGGCACGTGCG"  Picea_mature_1595  known    match to many species miR164 -> 5p in Arabidopsis
# "CATGTGCTCTTTTTCTCCATC"  Picea_mature_1251  known    miRBase_pab-miR164b -> 3p?
# "CACGTGCTCCCCTTCTCCAA"   Picea_mature_2626  variant  miRBase_pab-miR164a&miRBase_fve-miR164a-3p -> 3p
# "CACGTGCTCCCCTTCTCCAAC"  Picea_mature_4703  known    miRBase_pab-miR164a&miRBase_fve-miR164a-3p -> 3p
# 
# $miR399 good!
# "GGCAGCTCTCCTTTGGCAGG"  Picea_mature_923  known  miRBase_pab-miR399f&miRBase_pab-miR399g -> 5p
# "CGCCAAAGGAGAGTTGCTCAG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "TGCCAAAGGAGAAGTGCCCTG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
#  "GCCAAAGGAGAGTTGCCTT"  -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "CGCCAAAGGAGAGTTGCCCTG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "CGCCAAAGGAGAGTTGCGCTG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "CGCCAAAGGAGAGTTGCCT"   -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "TGCCAAAGGAGAGTTGCTCTG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# " GCCAAAGGAGAGTTGCCTTC" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "CGCCAAAGGAGAGTTGCCTT"  -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# "TGCCAAAGGAGAATTGCCCTG" -> 3p ath-miR399a-3p TGCCAAAGGAGATTTGCCCTG
# 
# $miR408 good! 
# "TGCACTGCCTCTTCCCTGGCT"  Picea_mature_2303  known  match to many species miR408 -> 3p in stu and vca
# "ACAGGGAAGAGTTAGGGCATG"  Picea_mature_4138  known  miRBase_pab-miR408 -> 5p canonical ACAGGGAACAAGCAGAGCATG Arabidopsis 5p

# $miR160 good! note as central-mismatch variant
#   "TGCCTGGCTCTCTGCATGCCA"  Picea_mature_1694  variant(known seq has CCC, not CTC)  miREN_Zos-miR160a&PmiREN_Nta-miR160d&miRBase_ptc-miR160h&PMRD_ptc-miR160h&PmiREN_Gma-miR160h&miRBase_nta-miR160d&PmiREN_Mac-miR160a
# "TCTGCCTGGCTCCCTGTATGC" 
#  "CTGCCTGGCTCCCTGTATGCC" 
#  "GTGCCTGGCTCCCTGTATGCC"
# "TGTGCCTGGCTCCCTGTATGC" 
#   "TGCCTGGCTCCCTGTATGCC" 
# 
# $miR171 all 3p, note as seed-shifted
# "TGATTGAGCCGTGCCAATATC"      Picea_mature_2839  known    match to many species miR171 -> 3p in ptc, zma, csi, osa, and several other species
#      "GAGCCGTGCCAATATCACAT"  Picea_mature_782   variant  PmiREN_Smo-miR171e&PMRD_smo-miR171b&PmiREN_Smo-miR171b&miRBase_smo-miR171b
#    "TTGAGCCGTGCCAATATCGCAC"  Picea_mature_1545  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#     "TGAGCCGTGCCAATATCGCA"   Picea_mature_1643  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#    "TTGAGCCGTGCCAATATCGC"    Picea_mature_1763  variant  PmiREN_Pab-miR171c&miRBase_pab-miR171a
#    "TTGAGCCGTGCCAATATCGCA"   Picea_mature_3003  known    PmiREN_Pab-miR171c&miRBase_pab-miR171a
#     "TGAGCCGTGCCAATATCACAT"  Picea_mature_4024  known    PmiREN_Smo-miR171e&PMRD_smo-miR171b&PmiREN_Smo-miR171b&miRBase_smo-miR171b
# 
# $miR2111 good! 
# "AACCCTAGAATGCAGATTATC"  Picea_mature_1142  variant  PmiREN_Mtr-miR2111y -> * version CCATTAGGCGTAGGGCTCCTG, CCATTTGACATAGGGCTCCTA in Arabidopsis, 3p GTCCTCGGGATGCGGATTACC
#. GATAATCTGCATTCTAGGGTT reverse complement 1142
#   "TAATCTGCATCCTGGGGTTTG"  Picea_mature_3134  known    miRBase_pab-miR2111&PmiREN_Pab-miR2111a -> canonical TAATCTGCATCCTGAGGTTTA, 5p TAATCTGCATCCTGAGGTTTA
# 
# $miR393 good!
# "ATCATGCAATCCCTTCGGAAT"   Picea_mature_290  variant  miRBase_sly-miR393&miRBase_stu-miR393-3p -> 3p ATCATGCTATCTCTTTGGATT
# "TCCAAAGGGATTGCATTGATTC"   -> 5p, canonical TCCAAAGGGATCGCATTGATCT, 5p TCCAAAGGGATCGCATTGATCC 
# "TCCAAAGGGATTGCATTGATT"   -> 5p, canonical TCCAAAGGGATCGCATTGATCT, 5p TCCAAAGGGATCGCATTGATCC 
# "TCCAAAGGGATTGCATTGATTCA" -> 5p, canonical TCCAAAGGGATCGCATTGATCT, 5p TCCAAAGGGATCGCATTGATCC 

# add the arm information back to the ones without, only if pid > 90
cross_arm_flags <- cross_arm_check %>%
  filter(pass_coverage, pid > 90) %>%
  arrange(desc(pid))

# build an arm-inference lookup from the flagged pairs
arm_inference <- cross_arm_flags %>%
  mutate(
    empty_id  = ifelse(arm1 == "", id1, ifelse(arm2 == "", id2, NA)),
    known_arm = ifelse(arm1 == "", arm2, ifelse(arm2 == "", arm1, NA))
  ) %>%
  filter(!is.na(empty_id))

# if an "" id happens to pair highly with BOTH 5p and 3p in different rows,
# that's genuinely ambiguous -- don't auto-fill those, flag for manual review
arm_inference_ambiguous <- arm_inference %>%
  group_by(empty_id) %>%
  filter(n_distinct(known_arm) > 1) %>%
  ungroup()
if (nrow(arm_inference_ambiguous) > 0) {
  message("Ambiguous -- not auto-filled, please check manually:")
  print(arm_inference_ambiguous %>% distinct(empty_id, known_arm))
}

arm_inference_clean <- arm_inference %>%
  distinct(empty_id, known_arm) %>%
  group_by(empty_id) %>%
  filter(n() == 1) %>%    # drops ids that appear in arm_inference_ambiguous
  ungroup()

# collapse back to matureDF's original
matureDF <- matureDF %>%
  left_join(arm_inference_clean, by = c("id" = "empty_id")) %>%
  mutate(arm      = coalesce(known_arm, arm)) %>%
  select(-known_arm)

# manual edit
test <- filter(matureDF, family %in% c(
  # "miR164"
  # "miR399"
  # "miR408"
  # "miR160"
  "miR171" #need further assessment
  # "miR2111"
  # "miR393"
  # "miR11408" #gymnosperm niche
  # "miR11487" #conifer niche
  ))

## miR171
# Picea_mature_2839 TGATTGAGCCGTGCCAATATC  3p 
# Picea_mature_782       GAGCCGTGCCAATATCACAT   3p   
# Picea_mature_1545    TTGAGCCGTGCCAATATCGCAC 3p   
# Picea_mature_1643     TGAGCCGTGCCAATATCGCA   3p   
# Picea_mature_1763    TTGAGCCGTGCCAATATCGC   3p   
# Picea_mature_3003    TTGAGCCGTGCCAATATCGCA  3p   
# Picea_mature_3471    TTGAGCCGCGCCAATATCACT  3p|5p
# Picea_mature_4024     TGAGCCGTGCCAATATCACAT  3p   
# ath-miR171c-3p:      TTGAGCCGTGCCAATATCACG
# ath-miR171c-5p:   AGATATTGGTGCGGTTCAATC
# Picea_mature_4855 AGATGCTGGTGCGGTTCAATC  5p   
# Picea_mature_285  AGATATTGGGGCGGTTCAATC  5p   
# Picea_mature_1470 AGATATTGGAGCGGTTCAAC   5p 

matureDF <- matureDF %>% 
  mutate(arm = case_when(
    id == "Picea_mature_1251" ~ "3p",
    id == "Picea_mature_923" ~ "5p",
    id %in% c("Picea_mature_192","Picea_mature_3955") ~ "3p",
    id == "Picea_mature_4138" ~ "5p",
    id == "Picea_mature_3134" ~ "5p",
    id == "Picea_mature_1142" ~ "3p",
    id %in% c("Picea_mature_1501","Picea_mature_2959","Picea_mature_4411") ~ "5p",
    id == "Picea_mature_3471" ~ "3p",
    .default = arm))
write_tsv(matureDF, here("data/analysis/results/allPredicted.tsv"))

# check family again
sanity_df <- matureDF %>%
  filter(!family %in% c("","NA")) %>%
  separate_rows(family, arm, sep = "\\|")
family_check <- sanity_df %>%
  group_by(family, arm) %>%
  summarise(n_members = n(), seqs = list(sequence), ids = list(id), .groups = "drop") %>%
  mutate(
    pairwise       = map2(seqs, ids, pairwise_stats),
    valid_pairwise = map(pairwise, ~ filter(.x, pass_coverage)),
    n_total_pairs  = map_int(pairwise, nrow),
    n_valid_pairs  = map_int(valid_pairwise, nrow),
    low_cov_pairs  = n_total_pairs - n_valid_pairs,     # pairs excluded for poor overlap
    min_pid        = map_dbl(valid_pairwise, ~ ifelse(nrow(.x) == 0, NA_real_, min(.x$pid))),
    mean_pid       = map_dbl(valid_pairwise, ~ ifelse(nrow(.x) == 0, NA_real_, mean(.x$pid)))
  ) %>%
  select(family, arm, n_members, n_total_pairs, n_valid_pairs, low_cov_pairs, min_pid, mean_pid) %>%
  arrange(min_pid)
saveRDS(family_check, here("data/analysis/results/family_check.rds"))








# Replicate check
repChecked <- read_lines(here("data/analysis/results/repChecked.txt"))

# expression data
samples <- bind_rows(
  read_tsv(here("doc/metadata_SE.tsv"), show_col_types=F) %>%
    dplyr::filter(! grepl("-unk",stage)) %>% 
    dplyr::rename(sample = NGI_ID,
                  sampleName = User_ID) %>%
    mutate(temp = gsub("PabK14-03","",sampleName),
           Stage = gsub("-unk","",stage)) %>%
    separate_wider_delim(temp, delim = "_", names = c("SampleType","rep")) %>%
    mutate(Dataset = "SE",
           Type = ifelse(stage %in% c("S7","S8"),"SDL","Sem"),
           StageType = paste0(Type, gsub("S","_",Stage)),
           Replicate = substr(rep,3,3),
           SampleName = paste0(Stage,"-",SampleType,"-",Replicate)) %>% 
    dplyr::select(-stage, -Mreads, -`>=Q30`, -Status, -rep, -replicate),
  read_csv(here("doc/metadataZE_mRNA_miRNA.csv"), show_col_types = FALSE) %>% 
    dplyr::filter(! is.na(miR_NGI_ID)) %>% 
    mutate(Dataset = "ZE",
           Type = gsub("ZE","Zem",Type),
           StageType = paste0(Type, "_", Stage),
           Stage = paste0("Z",Stage)) %>% 
    dplyr::select(-SampleID, -lincRNA_ID) %>% 
    dplyr::rename(sample = miR_NGI_ID,
                  sampleName = miR_User_ID,
                  description = SampleDescription)) %>% 
  mutate(StageType = factor(StageType,
                            levels = c(
                              paste0("SD_",1:3),
                              paste0("Zem_",c(4:7,"8-9",10)),
                              paste0("FG_",c(4:7,"8-9",10)),
                              paste0("Sem_",1:6),
                              paste0("SDL_",7:8))),
         Stage = factor(Stage, levels = c(paste0("Z",c(1:7,"8-9",10)), paste0("S",1:8))),
         Type = factor(Type, level = c("SD","Zem","FG","Sem","SDL"))) %>% 
  dplyr::filter(! sampleName %in% c("PabK14-03SM_SR1","PabK14-03RO_SR2","Pab738_2_FMG_MIR1"))

countAll <- read_tsv(here("data/analysis/results/bowtie/count_DE.txt"), show_col_types = F) %>%
  filter(! Gene_id == "Unaligned") %>%
  column_to_rownames("Gene_id") %>%
  setNames(gsub("_count","",basename(colnames(.)))) %>%
  as.matrix() %>% 
  .[repChecked, samples$sample] # Remove the lowQual samples, remove miR that do not have rep support

countSE <- countAll[, dplyr::filter(samples, Dataset == "SE")$sample] %>%
  .[rowSums(.) > 0,]
countZE <- countAll[, dplyr::filter(samples, Dataset == "ZE")$sample] %>%
  .[rowSums(.) > 0,]
countSem <- countSE[, dplyr::filter(samples, Type == "Sem")$sample] %>%
  .[rowSums(.) > 0,]
countSDL <- countSE[, dplyr::filter(samples, Type == "SDL")$sample] %>%
  .[rowSums(.) > 0,]
countSD <- countZE[, dplyr::filter(samples, Type == "SD")$sample] %>%
  .[rowSums(.) > 0,]
countZem <- countZE[, dplyr::filter(samples, Type == "Zem")$sample] %>%
  .[rowSums(.) > 0,]
countFG <- countZE[, dplyr::filter(samples, Type == "FG")$sample] %>%
  .[rowSums(.) > 0,]

# DE data
res <- c(
  readRDS(here("data/analysis/results/flipComp_1stTier_finalSE_ddsResults.rds"))[["names"]] %>% 
   transpose(),
  readRDS(here("data/analysis/results/flipComp_1stTier_finalZE_ddsResults.rds"))[["names"]] %>% 
   transpose())

# target data
targetPairsFinal <- readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairsFinal_withPval.rds")) %>% 
  select(miRNA,target) %>% 
  split(., .$miRNA) %>% 
  map(pull, target) %>% 
  map(paste, collapse="|") %>% 
  map(as.data.frame) %>% 
  bind_rows(.id = "id") %>% 
  dplyr::rename(target_bySeq = `.x[[i]]`)
targetPairsFinal_wSupport <- readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairsFinal_withPval.rds")) %>% 
  filter(!(is.na(selComp) & is.na(corSE) & is.na(corSDZE) & is.na(corSDFG))) %>% 
  select(miRNA,target) %>% 
  split(., .$miRNA) %>% 
  map(pull, target) %>% 
  map(paste, collapse="|") %>% 
  map(as.data.frame) %>% 
  bind_rows(.id = "id") %>% 
  dplyr::rename(target_withSupport = `.x[[i]]`)
targetPairsDE<- readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairs_finalMetadata_DEdirections.rds")) %>% 
  rowwise() %>% 
  mutate(miRNA_up_mRNA_dn = ifelse(!is.na(miRNA_up_mRNA_dn), paste0(miRNA_up_mRNA_dn, "_miRNA_up_mRNA_dn"), miRNA_up_mRNA_dn),
         miRNA_dn_mRNA_up = ifelse(!is.na(miRNA_dn_mRNA_up), paste0(miRNA_dn_mRNA_up, "_miRNA_dn_mRNA_up"), miRNA_dn_mRNA_up),
         miRNA_up_mRNA_up = ifelse(!is.na(miRNA_up_mRNA_up), paste0(miRNA_up_mRNA_up, "_miRNA_up_mRNA_up"), miRNA_up_mRNA_up),
         miRNA_dn_mRNA_dn = ifelse(!is.na(miRNA_dn_mRNA_dn), paste0(miRNA_dn_mRNA_dn, "_miRNA_dn_mRNA_dn"), miRNA_dn_mRNA_dn)) %>% 
  mutate(oppositeDE =  paste(c(na.omit(miRNA_up_mRNA_dn), na.omit(miRNA_dn_mRNA_up)), collapse = ","),
         sameDE =  paste(c(na.omit(miRNA_up_mRNA_up), na.omit(miRNA_dn_mRNA_dn)), collapse = ",")) %>% 
  ungroup() %>% 
  dplyr::select(miRNA,target, oppositeDE, sameDE)

# master table
dat <- matureDF %>% 
  filter(id %in% rownames(countAll[rowSums(countAll) > 0,])) %>% 
  mutate(expressSem = ifelse(id %in% rownames(countSem), "express", NA),
         expressSDL = ifelse(id %in% rownames(countSDL), "express", NA),
         expressSD = ifelse(id %in% rownames(countSD), "express", NA),
         expressFG = ifelse(id %in% rownames(countFG), "express", NA),
         expressZem = ifelse(id %in% rownames(countZem), "express", NA)) %>% 
  mutate(Sem2_vs_Sem1 = ifelse(id %in% res[["SE2_vs_SE1"]][["up"]], "up", ifelse(id %in% res[["SE2_vs_SE1"]][["dn"]], "down", NA)),
         Sem3_vs_Sem2 = ifelse(id %in% res[["SE3_vs_SE2"]][["up"]], "up", ifelse(id %in% res[["SE3_vs_SE2"]][["dn"]], "down", NA)),
         Sem4_vs_Sem3 = ifelse(id %in% res[["SE4_vs_SE3"]][["up"]], "up", ifelse(id %in% res[["SE4_vs_SE3"]][["dn"]], "down", NA)),
         Sem5_vs_Sem4 = ifelse(id %in% res[["SE5_vs_SE4"]][["up"]], "up", ifelse(id %in% res[["SE5_vs_SE4"]][["dn"]], "down", NA)),
         Sem6_vs_Sem5 = ifelse(id %in% res[["SE6_vs_SE5"]][["up"]], "up", ifelse(id %in% res[["SE6_vs_SE5"]][["dn"]], "down", NA)),
         SD2_vs_SD1 = ifelse(id %in% res[["SD2_vs_SD1"]][["up"]], "up", ifelse(id %in% res[["SD2_vs_SD1"]][["dn"]], "down", NA)),
         SD3_vs_SD2 = ifelse(id %in% res[["SD3_vs_SD2"]][["up"]], "up", ifelse(id %in% res[["SD3_vs_SD2"]][["dn"]], "down", NA)),
         Zem4_vs_SD3 = ifelse(id %in% res[["ZE4_vs_SD3"]][["up"]], "up", ifelse(id %in% res[["ZE4_vs_SD3"]][["dn"]], "down", NA)),
         Zem5_vs_Zem4 = ifelse(id %in% res[["ZE5_vs_ZE4"]][["up"]], "up", ifelse(id %in% res[["ZE5_vs_ZE4"]][["dn"]], "down", NA)),
         Zem6_vs_Zem5 = ifelse(id %in% res[["ZE6_vs_ZE5"]][["up"]], "up", ifelse(id %in% res[["ZE6_vs_ZE5"]][["dn"]], "down", NA)),
         Zem7_vs_Zem6 = ifelse(id %in% res[["ZE7_vs_ZE6"]][["up"]], "up", ifelse(id %in% res[["ZE7_vs_ZE6"]][["dn"]], "down", NA)),
         `Zem8-9_vs_Zem7` = ifelse(id %in% res[["ZE8-9_vs_ZE7"]][["up"]], "up", ifelse(id %in% res[["ZE8-9_vs_ZE7"]][["dn"]], "down", NA)),
         `Zem10_vs_Zem8-9`= ifelse(id %in% res[["ZE10_vs_ZE8-9"]][["up"]], "up", ifelse(id %in% res[["ZE10_vs_ZE8-9"]][["dn"]], "down", NA)),
         FG4_vs_SD3 = ifelse(id %in% res[["FG4_vs_SD3"]][["up"]], "up", ifelse(id %in% res[["FG4_vs_SD3"]][["dn"]], "down", NA)),
         FG5_vs_FG4 = ifelse(id %in% res[["FG5_vs_FG4"]][["up"]], "up", ifelse(id %in% res[["FG5_vs_FG4"]][["dn"]], "down", NA)),
         FG6_vs_FG5 = ifelse(id %in% res[["FG6_vs_FG5"]][["up"]], "up", ifelse(id %in% res[["FG6_vs_FG5"]][["dn"]], "down", NA)),
         FG7_vs_FG6 = ifelse(id %in% res[["FG7_vs_FG6"]][["up"]], "up", ifelse(id %in% res[["FG7_vs_FG6"]][["dn"]], "down", NA)),
         `FG8-9_vs_FG7` = ifelse(id %in% res[["FG8-9_vs_FG7"]][["up"]], "up", ifelse(id %in% res[["FG8-9_vs_FG7"]][["dn"]], "down", NA)),
         `FG10_vs_FG8-9`= ifelse(id %in% res[["FG10_vs_FG8-9"]][["up"]], "up", ifelse(id %in% res[["FG10_vs_FG8-9"]][["dn"]], "down", NA))) %>% 
  left_join(targetPairsFinal, by = "id") %>% 
  left_join(targetPairsFinal_wSupport, by = "id") %>% 
  mutate(family = gsub("\\|miRN[0-9]+$","",gsub("^miRN[0-9]+\\|","",family)))
dat[is.na(dat)] <- ""
write_csv(dat,here("data/analysis/results/masterTable.csv"))

# recheck
# test <- dat %>% 
#   filter(annotation %in% c("known", "variant") | !is.na(P_criteria))
# sum(!is.na(test$expressSem))
# sum(!is.na(test$SE7_vs_SE6))

# target table
readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairsFinal_withPval.rds")) %>% 
  filter(!(is.na(selComp) & is.na(corSE) & is.na(corSDZE) & is.na(corSDFG))) %>% 
  select(-isSE, -isZE) %>% 
  mutate(selComp = gsub(",","|", selComp)) %>% 
  dplyr::rename(oppositeDE = selComp) %>% 
  write_csv(here("data/analysis/results/targetPairs_withSupport.csv"))

# # filter target table
# mAnnot <- read_tsv("reference/Picea-abies/v2.0/annotation/Picab02_230926_at01_longest_representative_annotations_merged_sorted.tsv.gz", col_names = T, show_col_types = F) %>%
#   mutate(target = map_chr(id, str_split_i, "\\.", 1)) %>% 
#   select(-id, -eggnog_seed_ortholog, -eggnog_evalue, -eggnog_score, -interpro_evalue, -interpro_panther_evalue) %>% 
#   split(., .$target) %>% 
#   map( ~ summarise(.x, across(everything(), ~ paste0(unique(na.omit(.x)), collapse = ",")))) %>% 
#   bind_rows()
# write_tsv(mAnnot, here("data/analysis/PicabAnnot_collapsedPanther.tsv"))

mAnnot <- read_tsv(here("data/analysis/PicabAnnot_collapsedPanther.tsv"), show_col_types = FALSE)

# remove no annot
# dat <- readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairsFinal_withPval.rds")) %>% 
#   filter(!(is.na(selComp) & is.na(corSE) & is.na(corSDZE) & is.na(corSDFG))) %>% 
#   dplyr::rename(oppositeDE = selComp) %>% 
#   rowwise() %>% 
#   mutate(negCorCoeff = paste(c("Sem","SDZem","SDFG")[c(!is.na(corSE),!is.na(corSDZE), !is.na(corSDFG))], collapse = ","),
#          oppositeDE = gsub("SE","Sem", gsub("ZE","Zem", oppositeDE))) %>% 
#   ungroup() %>% 
#   select(-NumberOfPrograms,-isSE, -isZE, -pvalSE,-padjSE,-pvalSDZE, -padjSDZE, -pvalSDFG, -padjSDFG, -corSE,-corSDZE,-corSDFG) %>% 
#   left_join(matureDF, by = c("miRNA" = "id")) %>% 
#   select(-database) %>% 
#   inner_join(mAnnot, by = "target") %>% 
#   select(-eggnog_go, -eggnog_KEGG_ko, -eggnog_KEGG_Pathway, -eggnog_KEGG_Module, -eggnog_KEGG_Reaction, -eggnog_KEGG_rclass, -interpro_pfam, -interpro_id, -interpro_go, -interpro_panther) %>% 
#   filter(! if_all(c(eggnog_description,eggnog_preferred_name,eggnog_pfams,interpro_description,interpro_panther_description), ~ .x == "")) %>% 
#   rename(miRNA_seq = sequence,
#          known_miRNA_id = description,
#          known_miRNA_seq = knownSeq)
# dat[is.na(dat)] <- ""

# keep no annot
dat <- readRDS(here("data/analysis/results/flipComp_1stTier_TargetPairsFinal_withPval.rds")) %>% 
  filter(!(is.na(selComp) & is.na(corSE) & is.na(corSDZE) & is.na(corSDFG))) %>% 
  left_join(targetPairsDE,  by = c("miRNA","target")) %>% 
  rowwise() %>% 
  mutate(negCorCoeff = paste(c("Sem","SDZem","SDFG")[c(!is.na(corSE),!is.na(corSDZE), !is.na(corSDFG))], collapse = ","),
         oppositeDE = gsub("SE","Sem", gsub("ZE","Zem", oppositeDE))) %>% 
  ungroup() %>% 
  select(-NumberOfPrograms,-isSE, -isZE, -pvalSE,-padjSE,-pvalSDZE, -padjSDZE, -pvalSDFG, -padjSDFG, -corSE,-corSDZE,-corSDFG) %>%
  left_join(matureDF, by = c("miRNA" = "id")) %>% 
  select(-database) %>% 
  left_join(mAnnot, by = "target") %>% 
  select(-eggnog_go, -eggnog_KEGG_ko, -eggnog_KEGG_Pathway, -eggnog_KEGG_Module, -eggnog_KEGG_Reaction, -eggnog_KEGG_rclass, -interpro_pfam, -interpro_id, -interpro_go, -interpro_panther) %>% 
  rename(miRNA_seq = sequence,
         known_miRNA_id = description,
         known_miRNA_seq = knownSeq)
dat[is.na(dat)] <- ""

write_tsv(dat, here("data/analysis/results/targetPairs_withSupport_annot.tsv"))



