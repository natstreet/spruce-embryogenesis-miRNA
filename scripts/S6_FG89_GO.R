## Which functions underlie the mRNA similarity of S6 (after the desiccation treatment) to the late FG (FG8-9)?
## Genes: the 5,000 genes of Figure 5. Per-gene contribution to Spearman rho(S6, FG8-9) = product of centred ranks.
## Test sets: top 10% contributors (500) split into genes high in both and genes low in both.
## Background: the 5,000 genes. topGO parent-child, Fisher, BH <= 0.05 (as in the Methods).
suppressPackageStartupMessages({library(topGO); library(dplyr)})
ZE_ORDER <- c("SD1","SD2","SD3",paste0("Zem",c(4:7,"8-9",10)),paste0("FG",c(4:7,"8-9",10)))
mz <- read.csv("data/metadataZE_mRNA_miRNA.csv"); ms <- read.csv("data/metadataSE_mRNA_miRNA.csv")
mz$grp <- ifelse(mz$Type=="SD",paste0("SD",mz$Stage),paste0(ifelse(mz$Type=="ZE","Zem","FG"),mz$Stage)); ms$grp <- paste0("S",ms$Stage)
sm <- function(X,g,l) sapply(l,function(x) rowMeans(X[,g==x,drop=FALSE]))
mze <- readRDS("data/mRNA_vst_ZE.rds"); mse <- readRDS("data/mRNA_vst_SE.rds")
colnames(mse) <- sub("_S\\d+$","",colnames(mse))
Zm <- sm(mze, mz$grp[match(colnames(mze), mz$SampleID)], ZE_ORDER); Sm <- sm(mse, ms$grp[match(colnames(mse), ms$SampleID)], paste0("S",1:6))
expressed <- rowMeans(Zm) > 1 & rowMeans(Sm) > 1
score <- apply(Zm - rowMeans(Zm), 1, var) + apply(Sm - rowMeans(Sm), 1, var); score[!expressed] <- -Inf
keep <- order(score, decreasing = TRUE)[1:5000]
cZ <- (Zm - rowMeans(Zm))[keep, ]; cS <- (Sm - rowMeans(Sm))[keep, ]
rk <- function(x){r<-rank(x); r-mean(r)}
res <- list()
for (z in c("FG8-9", "Zem8-9")) {
  p <- rk(cS[,"S6"]) * rk(cZ[,z]); p <- p / sqrt(sum(rk(cS[,"S6"])^2) * sum(rk(cZ[,z])^2))
  top <- order(p, decreasing = TRUE)[1:500]
  res[[paste(z,"high")]] <- rownames(cZ)[top][cS[top,"S6"] > 0]
  res[[paste(z,"low")]]  <- rownames(cZ)[top][cS[top,"S6"] < 0]
  cat(z, "rho", round(sum(p),3), "| top 500 share", round(sum(p[top])/sum(p),2), "| high", length(res[[paste(z,"high")]]), "low", length(res[[paste(z,"low")]]), "\n")
}
ga <- read.delim(gzfile("data/GO_annotation_gopher_20261002.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
ga$V1 <- sub("\\.mRNA\\.\\d+$", "", ga$V1)
g2go <- setNames(strsplit(ga$V2, "\\|"), ga$V1)
bg <- rownames(cZ)
out <- list()
for (nm in names(res)) for (o in c("BP","MF","CC")) {
  al <- factor(as.integer(bg %in% res[[nm]])); names(al) <- bg
  gd <- suppressMessages(new("topGOdata", ontology = o, allGenes = al, annot = annFUN.gene2GO, gene2GO = g2go))
  r <- suppressMessages(runTest(gd, algorithm = "parentchild", statistic = "fisher"))
  bh <- p.adjust(score(r), "BH"); n <- sum(bh <= 0.05)
  if (n > 0) { t <- GenTable(gd, r, topNodes = n, numChar = 80); t <- t[, 1:5]; t$P <- signif(score(r)[t$GO.ID], 2); t$BH <- signif(bh[t$GO.ID], 2); t$set <- nm; t$ont <- o; out[[paste(nm,o)]] <- t }
}
tab <- do.call(rbind, out)
if (is.null(tab)) cat("no terms at BH <= 0.05\n") else {
  write.table(tab, "results/S6_FG89_GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  print(tab[, c("set","ont","Term","Annotated","Significant","Expected","BH")], row.names = FALSE) }
