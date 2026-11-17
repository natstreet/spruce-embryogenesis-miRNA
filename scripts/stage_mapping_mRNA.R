## Stage mapping of the mRNA profiles (Figure 5, lower panel).
## Inputs: data/mRNA_vst_ZE.rds, data/mRNA_vst_SE.rds (column names carry _S<n>),
## data/metadataZE_mRNA_miRNA.csv, data/metadataSE_mRNA_miRNA.csv. Output: results/stage_mapping_mRNA.tsv
## Rule: average vst per stage; keep genes whose mean of stage means is > 1 in
## each dataset (26,598 genes); score each gene by the variance of its centred
## stage means in ZE (15 groups) plus that in SE (S1-S6); take the top 5,000;
## centre each gene within each dataset; Spearman correlation of each SE stage
## with each ZE group.
mze <- readRDS("data/mRNA_vst_ZE.rds"); mse <- readRDS("data/mRNA_vst_SE.rds")
mz <- read.csv("data/metadataZE_mRNA_miRNA.csv"); ms <- read.csv("data/metadataSE_mRNA_miRNA.csv")
colnames(mse) <- sub("_S\\d+$", "", colnames(mse))
mz$grp <- ifelse(mz$Type == "SD", paste0("SD", mz$Stage),
                 paste0(ifelse(mz$Type == "ZE", "Zem", "FG"), mz$Stage))
ms$grp <- paste0("S", ms$Stage)
ZE_ORDER <- c("SD1", "SD2", "SD3", paste0("Zem", c(4:7, "8-9", 10)),
              paste0("FG", c(4:7, "8-9", 10)))
SEL <- paste0("S", 1:6)
smean <- function(X, g, lev) sapply(lev, function(l) rowMeans(X[, g == l, drop = FALSE]))
Zm <- smean(mze, mz$grp[match(colnames(mze), mz$SampleID)], ZE_ORDER)
Sm <- smean(mse, ms$grp[match(colnames(mse), ms$SampleID)], SEL)
expressed <- rowMeans(Zm) > 1 & rowMeans(Sm) > 1
score <- apply(Zm - rowMeans(Zm), 1, var) + apply(Sm - rowMeans(Sm), 1, var)
score[!expressed] <- -Inf
keep <- order(score, decreasing = TRUE)[1:5000]
M <- cor(Sm[keep, ] - rowMeans(Sm[keep, ]), Zm[keep, ] - rowMeans(Zm[keep, ]),
         method = "spearman")
write.table(data.frame(SE = rownames(M), round(M, 3), check.names = FALSE),
            "results/stage_mapping_mRNA.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("expressed in both:", sum(expressed), "\n")
