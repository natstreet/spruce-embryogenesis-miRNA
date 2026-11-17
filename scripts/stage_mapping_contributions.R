## Decomposition of the miRNA stage mapping (Figure 5, Supplementary Figure S5B) into the
## contributions of single miRNAs and of the enrichment classes.
## Same mapping rule as Figure 5: stage means of vst, centred per miRNA within each dataset,
## Spearman correlation of each SE stage with each ZE group. A Spearman correlation is the
## sum over features of the products of centred ranks, divided by the norm of the two rank
## vectors, so each miRNA's product is its contribution to the correlation.
## Output: results/stage_mapping_contributions.txt
## (1) S5 mapping for each S5 replicate library
## (2) miRNAs contributing most to rho(S5, Zem10) - rho(S5, Zem7) and - rho(S5, Zem8-9):
##     their mean centred values in the zygotic and somatic stages
## (3) rho(S6, FG8-9) and rho(S6, Zem10) split by enrichment class, and the mean centred
##     values of each class at S4, S5 and S6
ZE_ORDER <- c("SD1", "SD2", "SD3", paste0("Zem", c(4:7, "8-9", 10)), paste0("FG", c(4:7, "8-9", 10)))
SEL <- paste0("S", 1:6)
mz <- read.csv("data/metadataZE_mRNA_miRNA.csv"); ms <- read.csv("data/metadataSE_mRNA_miRNA.csv")
mz$grp <- ifelse(mz$Type == "SD", paste0("SD", mz$Stage), paste0(ifelse(mz$Type == "ZE", "Zem", "FG"), mz$Stage))
ms$grp <- paste0("S", ms$Stage)
smean <- function(X, g, lev) sapply(lev, function(l) rowMeans(X[, g == l, drop = FALSE]))
a <- readRDS("data/miRNA_vst_ZE.rds"); b <- readRDS("data/miRNA_vst_SE.rds")
cls <- read.delim("results/Zem_vs_FG_enrichment.tsv")
a <- a[rownames(a) %in% cls$miR, ]; b <- b[rownames(a), ]
gz <- mz$grp[match(colnames(a), mz$miR_User_ID)]; gs <- ms$grp[match(colnames(b), ms$miR_User_ID)]
Zm <- smean(a, gz, ZE_ORDER); Sm <- smean(b, gs, SEL)
cZ <- Zm - rowMeans(Zm); cS <- Sm - rowMeans(Sm)
rk <- function(x) { r <- rank(x); r - mean(r) }
contrib <- function(s, z) { p <- rk(cS[, s]) * rk(cZ[, z]); p / sqrt(sum(rk(cS[, s])^2) * sum(rk(cZ[, z])^2)) }
cl <- cls$enriched_deseq[match(rownames(cZ), cls$miR)]   # Zem, FG or ns

sink("results/stage_mapping_contributions.txt")
cat("miRNAs:", nrow(cZ), "\n")
cat("\n## (1) Spearman correlation of each S5 replicate (centred on the SE stage means) with the zygotic stage means\n")
R <- b[, gs == "S5", drop = FALSE] - rowMeans(Sm)
print(round(cor(R, cZ[, c("Zem4", "Zem5", "Zem6", "Zem7", "Zem8-9", "Zem10")], method = "spearman"), 2))

cat("\n## (2) miRNAs contributing most to the difference between S5 mapping on Zem10 and on Zem7 or Zem8-9\n")
for (z in c("Zem7", "Zem8-9")) {
  dlt <- contrib("S5", "Zem10") - contrib("S5", z)
  top <- order(dlt, decreasing = TRUE)[1:50]
  cat("\nrho(S5, Zem10) - rho(S5,", z, ") =", round(sum(dlt), 3),
      "| share from the 50 largest contributors:", round(sum(dlt[top]) / sum(dlt), 2), "\n")
  cat("mean centred vst of these 50 miRNAs:\n")
  print(round(colMeans(cbind(cZ[top, c("Zem5", "Zem6", "Zem7", "Zem8-9", "Zem10")], cS[top, c("S4", "S5", "S6")])), 2))
  cat("enrichment class of these 50 miRNAs:\n"); print(table(cl[top]))
}

cat("\n## (3) S6 mapping split by enrichment class (sum of contributions)\n")
for (pr in list(c("S6", "FG8-9"), c("S6", "Zem10"))) {
  p <- contrib(pr[1], pr[2])
  cat("\nrho(", pr[1], ",", pr[2], ") =", round(sum(p), 3), "; by class:\n"); print(round(tapply(p, cl, sum), 3))
}
cat("\nmean centred vst by enrichment class at S4, S5 and S6:\n")
print(round(sapply(c("S4", "S5", "S6"), function(s) tapply(cS[, s], cl, mean)), 2))
sink()
cat(readLines("results/stage_mapping_contributions.txt"), sep = "\n")
