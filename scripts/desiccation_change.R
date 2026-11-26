## Change of the 147 embryo-enriched (Zem-enriched) miRNAs at the somatic desiccation treatment (S5 -> S6) and over
## the period of zygotic desiccation in the cone (Zem6 -> Zem8-9; stages 7 to 9 are the desiccated seed in
## Stojkovic et al. 2024, which used the same samples), with the steps Zem6 -> Zem7, Zem7 -> Zem8-9 and
## Zem8-9 -> Zem10 (desiccated seed in the cone to dry seed of the previous year), and the embryo signature index
## per stage.
## Change = stage-mean vst of the later stage minus that of the earlier stage, for each miRNA.
## Inputs: data/miRNA_vst_ZE.rds, data/miRNA_vst_SE.rds, data/metadata*_mRNA_miRNA.csv,
##   results/Zem_vs_FG_enrichment.tsv, results/embryo_signature_index.tsv
## Output: results/desiccation_change.tsv
vZ <- readRDS("data/miRNA_vst_ZE.rds"); vS <- readRDS("data/miRNA_vst_SE.rds")
if (!is.matrix(vZ)) vZ <- SummarizedExperiment::assay(vZ)
if (!is.matrix(vS)) vS <- SummarizedExperiment::assay(vS)
mz <- read.csv("data/metadataZE_mRNA_miRNA.csv"); ms <- read.csv("data/metadataSE_mRNA_miRNA.csv")
gz <- with(mz, setNames(ifelse(Type == "SD", paste0("SD", Stage), paste0(ifelse(Type %in% c("ZE", "Zem"), "Zem", "FG"), Stage)), miR_User_ID))
gs <- with(ms, setNames(paste0("S", Stage), miR_User_ID))
smean <- function(v, g, st) { lib <- colnames(v)[colnames(v) %in% names(g) & g[colnames(v)] == st]
  stopifnot(length(lib) > 0); rowMeans(v[, lib, drop = FALSE]) }
enr <- read.delim("results/Zem_vs_FG_enrichment.tsv")
emb <- enr$miR[enr$enriched_deseq == "Zem"]; stopifnot(length(emb) == 147, all(emb %in% rownames(vZ)), all(emb %in% rownames(vS)))
steps <- list("S5->S6" = c("S", "S5", "S6"), "Zem6->Zem7" = c("Z", "Zem6", "Zem7"), "Zem7->Zem8-9" = c("Z", "Zem7", "Zem8-9"),
              "Zem6->Zem8-9" = c("Z", "Zem6", "Zem8-9"), "Zem8-9->Zem10" = c("Z", "Zem8-9", "Zem10"))
out <- do.call(rbind, lapply(names(steps), function(nm) { s <- steps[[nm]]
  v <- if (s[1] == "S") vS else vZ; g <- if (s[1] == "S") gs else gz
  d <- smean(v, g, s[3])[emb] - smean(v, g, s[2])[emb]
  data.frame(step = nm, n_miRNAs = length(d), median_change = round(median(d), 2), n_decreased = sum(d < 0)) }))
idx <- read.delim("results/embryo_signature_index.tsv")
im <- aggregate(index ~ stage, idx[idx$tissue %in% c("Zem", "Sem") | grepl("^S[0-9]", idx$stage), ], mean)
write.table(out, "results/desiccation_change.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out); print(im[im$stage %in% c(paste0("Zem", c(4:7, "8-9", 10)), paste0("S", 1:6)), ], digits = 3)
