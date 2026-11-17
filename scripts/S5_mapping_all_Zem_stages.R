## Supplementary Figure S5B: correlation of S5 with every zygotic embryo stage (Zem4 to Zem10)
## for each feature set. Same rule as Figure 5: stage means of vst,
## centred per feature within each dataset over all its groups, Spearman correlation.
## Gene sets:
##   all_expressed: mean of stage means > 1 in each dataset
##   targets_focal: predicted targets (26,812 pairs) of the 961 retained miRNAs
##   targets_supported: targets in the 718 pairs with expression support
##   targets_top10pct_contrib: targets of the 96 miRNAs contributing most to
##       rho(S5, Zem10) - rho(S5, Zem6)
##   targets_Zem10_DE: targets of the miRNAs DE at Zem8-9 -> Zem10
ZE_ORDER <- c("SD1", "SD2", "SD3", paste0("Zem", c(4:7, "8-9", 10)), paste0("FG", c(4:7, "8-9", 10)))
ZEM <- paste0("Zem", c(4:7, "8-9", 10)); SEL <- paste0("S", 1:6)
mz <- read.csv("data/metadataZE_mRNA_miRNA.csv"); ms <- read.csv("data/metadataSE_mRNA_miRNA.csv")
mz$grp <- ifelse(mz$Type == "SD", paste0("SD", mz$Stage), paste0(ifelse(mz$Type == "ZE", "Zem", "FG"), mz$Stage))
ms$grp <- paste0("S", ms$Stage)
smean <- function(X, g, lev) sapply(lev, function(l) rowMeans(X[, g == l, drop = FALSE]))
ctr <- function(M) M - rowMeans(M)
## miRNA
a <- readRDS("data/miRNA_vst_ZE.rds"); b <- readRDS("data/miRNA_vst_SE.rds")
D1 <- read.delim("results/Zem_vs_FG_enrichment.tsv"); a <- a[rownames(a) %in% D1$miR, ]; b <- b[rownames(a), ]
cZm <- ctr(smean(a, mz$grp[match(colnames(a), mz$miR_User_ID)], ZE_ORDER))
cSm <- ctr(smean(b, ms$grp[match(colnames(b), ms$miR_User_ID)], SEL))
rk <- function(x) { r <- rank(x); r - mean(r) }
ctb <- function(z) rk(cSm[, "S5"]) * rk(cZm[, z]) / sqrt(sum(rk(cSm[, "S5"])^2) * sum(rk(cZm[, z])^2))
Dm <- ctb("Zem10") - ctb("Zem6"); o <- order(Dm, decreasing = TRUE)
keep20 <- rownames(cSm)[-o[1:20]]
## mRNA
mze <- readRDS("data/mRNA_vst_ZE.rds"); mse <- readRDS("data/mRNA_vst_SE.rds")
colnames(mse) <- sub("_S\\d+$", "", colnames(mse))
Zg <- smean(mze, mz$grp[match(colnames(mze), mz$SampleID)], ZE_ORDER)
Sg <- smean(mse, ms$grp[match(colnames(mse), ms$SampleID)], SEL)
stopifnot(identical(rownames(Zg), rownames(Sg)))
g_all <- rownames(Zg)[rowMeans(Zg) > 1 & rowMeans(Sg) > 1]
cZg <- ctr(Zg[g_all, ]); cSg <- ctr(Sg[g_all, ])
tp <- read.delim("data/predicted_target_pairs.tsv", stringsAsFactors = FALSE); tp$target <- sub("\\.mRNA\\.\\d+$", "", tp$target)
sup <- read.delim("results/supported_target_pairs_annotated.tsv"); sup$target <- sub("\\.mRNA\\.\\d+$", "", sup$target)
z10 <- read.csv("results/miRNA_DE/ZE10_vs_ZE8-9_miRresults.csv"); z10DE <- z10$miR[!is.na(z10$padj) & z10$padj < 0.01]
sets <- list(all_expressed = g_all,
             targets_focal = intersect(unique(tp$target[tp$miRNA %in% D1$miR]), g_all),
             targets_supported = intersect(unique(sup$target), g_all),
             targets_top10pct_contrib = intersect(unique(tp$target[tp$miRNA %in% rownames(cSm)[o][1:96]]), g_all),
             targets_Zem10_DE = intersect(unique(tp$target[tp$miRNA %in% z10DE]), g_all))
rowm <- function(cs, cz, f, layer, set) data.frame(layer = layer, set = set, n = length(f),
          t(sapply(ZEM, function(z) cor(cs[f, "S5"], cz[f, z], method = "spearman"))), check.names = FALSE)
TS <- rbind(rowm(cSm, cZm, rownames(cSm), "miRNA", "focal_miRNAs"),
            rowm(cSm, cZm, keep20, "miRNA", "without_top20_contributors"),
            do.call(rbind, lapply(names(sets), function(s) rowm(cSg, cZg, sets[[s]], "mRNA", s))))
TS[, ZEM] <- round(TS[, ZEM], 4)
write.table(TS, "results/S5_mapping_all_Zem_stages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(TS, digits = 2)
