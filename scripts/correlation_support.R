## Correlation support for predicted miRNA-target pairs, using the stage-to-stage miRNA DE sets.
## Pearson correlation of vst values within each series (Sem S1-S6; SD+Zem; SD+FG). In each series every
## miRNA DE in either dataset is correlated with every predicted target gene, and BH is applied over all
## these combinations. A pair is kept in a series only if its miRNA is DE in the same dataset (ZE or SE).
## Inputs: data/predicted_target_pairs.tsv, results/miRNA_DE, data/miRNA_vst_*.rds, data/mRNA_vst_*.rds
##   and data/metadata*_mRNA_miRNA.csv.
## Output: results/target_pair_correlations.tsv
tpF <- read.delim("data/predicted_target_pairs.tsv", stringsAsFactors = FALSE)
vmZ <- readRDS("data/miRNA_vst_ZE.rds"); vmS <- readRDS("data/miRNA_vst_SE.rds")
vgZ <- readRDS("data/mRNA_vst_ZE.rds"); vgS <- readRDS("data/mRNA_vst_SE.rds")
mdZ <- read.csv("data/metadataZE_mRNA_miRNA.csv"); mdS <- read.csv("data/metadataSE_mRNA_miRNA.csv")
cS <- sub("_S[0-9]+$", "", colnames(vgS))
mapS <- mdS[mdS$miR_User_ID %in% colnames(vmS) & mdS$SampleID %in% cS & mdS$Stage <= 6, ]
mapZ <- mdZ[mdZ$miR_User_ID %in% colnames(vmZ) & mdZ$SampleID %in% colnames(vgZ), ]
series <- list(SDZE = mapZ[mapZ$Type %in% c("SD","ZE","Zem"), ], SDFG = mapZ[mapZ$Type %in% c("SD","FG"), ], SE = mapS)
zs <- function(X) { X <- X - rowMeans(X); X / sqrt(rowSums(X^2)) }
Z <- lapply(setNames(names(series), names(series)), function(s) { m <- series[[s]]
  if (s == "SE") list(mi = zs(vmS[, m$miR_User_ID]), g = zs(vgS[, match(m$SampleID, cS)]), n = nrow(m))
  else list(mi = zs(vmZ[, m$miR_User_ID]), g = zs(vgZ[, m$SampleID]), n = nrow(m)) })
pfromr <- function(r, n) { t <- r * sqrt((n - 2) / (1 - r^2)); 2 * pt(-abs(t), n - 2) }
comps <- setdiff(sub("_miRresults.csv", "", list.files("results/miRNA_DE", "_miRresults.csv$")), c("SE6_vs_SE7","SE7_vs_SE8"))
de <- function(cc) { x <- read.csv(file.path("results/miRNA_DE", paste0(cc, "_miRresults.csv"))); setdiff(x$miR[!is.na(x$padj) & x$padj < 0.01 & abs(x$log2FoldChange) >= 0.5], "Unaligned") }
D <- lapply(setNames(comps, comps), de)
deSE <- unique(unlist(D[grep("^SE", comps)])); deZE <- unique(unlist(D[grep("^SE", comps, invert = TRUE)])); deU <- union(deSE, deZE)
keepset <- list(SE = deSE, SDZE = deZE, SDFG = deZE)
o <- data.frame(miRNA = tpF$miRNA, target = tpF$target)
for (s in c("SE","SDZE","SDFG")) { z <- Z[[s]]; m <- intersect(deU, rownames(z$mi)); g <- intersect(unique(tpF$target), rownames(z$g))
  R <- z$mi[m, ] %*% t(z$g[g, ]); PA <- matrix(p.adjust(pfromr(R, z$n), "BH"), nrow(R))
  i <- cbind(match(tpF$miRNA, m), match(tpF$target, g)); ok <- !is.na(i[, 1]) & !is.na(i[, 2]) & tpF$miRNA %in% keepset[[s]]
  r <- pa <- rep(NA_real_, nrow(tpF)); r[ok] <- R[i[ok, , drop = FALSE]]; pa[ok] <- PA[i[ok, , drop = FALSE]]
  o[[paste0("cor", s)]] <- r; o[[paste0("padj", s)]] <- pa
  cat(s, "negatively correlated pairs (padj < 0.05):", sum(!is.na(r) & r < 0 & pa < 0.05), "\n") }
write.table(o, "results/target_pair_correlations.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
