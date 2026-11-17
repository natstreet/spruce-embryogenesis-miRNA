## Stage-to-stage differential expression with limma voomWithQualityWeights, sample weights
## estimated per stage (var.group = stage), for mRNA and miRNA in both series.
## mRNA: TMM, filterByExpr, eBayes, BH < 0.01 and |log2FC| > 0.5.
## miRNA: library sizes = DESeq2 size factors (with the Unaligned row) x geometric mean of library sizes,
##        filterByExpr, treat against |log2FC| = 0.5, BH < 0.01.
suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR)})
## Inputs: data/mRNA_counts_ZE.tsv.gz, data/mRNA_counts_SE.tsv.gz and data/metadata*_mRNA_miRNA.csv (mRNA),
##   data/miRNA_counts_DESeqDataSet.rds (miRNA, with the Unaligned row).
## Outputs: results/miRNA_DE/*_miRresults.csv, results/mRNA_DE_voomQWB.tsv.gz, results/mRNA_DE_names_voomQWB.rds, results/stage_weights.rds
MI <- "data/miRNA_counts_DESeqDataSet.rds"; OUTD <- "results"
dir.create(file.path(OUTD, "miRNA_DE"), recursive = TRUE, showWarnings = FALSE)
lv <- function(x) sub("^([A-Z]+)(.*)$", "\\1_\\2", x)
comps <- list(ZE = c("SD2_vs_SD1","SD3_vs_SD2","ZE4_vs_SD3","ZE5_vs_ZE4","ZE6_vs_ZE5","ZE7_vs_ZE6","ZE8-9_vs_ZE7","ZE10_vs_ZE8-9",
                     "FG4_vs_SD3","FG5_vs_FG4","FG6_vs_FG5","FG7_vs_FG6","FG8-9_vs_FG7","FG10_vs_FG8-9"),
              SE = paste0("SE", 2:6, "_vs_SE", 1:5))
mk <- function(x) make.names(x)
fitq <- function(cnt, grp, libsize, cc, layer) {
  grp <- droplevels(factor(grp)); y <- DGEList(cnt, group = grp, lib.size = libsize)
  if (layer == "mRNA") y <- calcNormFactors(y)
  keep <- filterByExpr(y, group = grp); y <- y[keep, , keep.lib.sizes = TRUE]
  levels(grp) <- mk(levels(grp)); design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  cm <- makeContrasts(contrasts = sapply(cc, function(z) { p <- strsplit(z, "_vs_")[[1]]; paste0(mk(lv(p[1])), "-", mk(lv(p[2]))) }), levels = design)
  colnames(cm) <- cc
  v <- voomWithQualityWeights(y, design, var.group = grp)
  f <- contrasts.fit(lmFit(v, design), cm)
  f <- if (layer == "mRNA") eBayes(f) else treat(f, fc = 2^0.5)
  tabs <- lapply(setNames(cc, cc), function(z) if (layer == "mRNA") topTable(f, coef = z, number = Inf, sort.by = "none") else topTreat(f, coef = z, number = Inf, sort.by = "none"))
  w <- tapply(v$targets$sample.weights, grp, mean)
  list(tabs = tabs, keep = keep, weights = w, grp = grp) }
call <- function(tt, layer) { g <- rownames(tt); k <- tt$adj.P.Val < 0.01
  if (layer == "mRNA") k <- k & abs(tt$logFC) > 0.5 & grepl("^PA_", g) else k <- k & g != "Unaligned"
  list(up = g[k & tt$logFC > 0], dn = g[k & tt$logFC < 0]) }
## mRNA
NM <- list(up = list(), dn = list()); RES <- list(); WT <- list()
for (s in c("ZE","SE")) {
  cnt <- as.matrix(read.delim(gzfile(sprintf("data/mRNA_counts_%s.tsv.gz", s)), row.names = 1, check.names = FALSE))
  md <- read.csv(sprintf("data/metadata%s_mRNA_miRNA.csv", s)); md <- md[match(colnames(cnt), md$SampleID), ]
  stopifnot(!anyNA(md$SampleID)); stage <- paste0(md$Type, "_", md$Stage)   # e.g. SD_1, ZE_4, FG_8-9, SE_1
  r <- fitq(cnt, stage, colSums(cnt), comps[[s]], "mRNA"); WT[[paste0("mRNA_", s)]] <- r$weights
  for (cc in comps[[s]]) { x <- call(r$tabs[[cc]], "mRNA"); NM$up[[cc]] <- x$up; NM$dn[[cc]] <- x$dn; RES[[cc]] <- r$tabs[[cc]] } }
saveRDS(NM, file.path(OUTD, "mRNA_DE_names_voomQWB.rds")); g <- gzfile(file.path(OUTD, "mRNA_DE_voomQWB.tsv.gz"), "w")
write.table(do.call(rbind, lapply(names(RES), function(cc) data.frame(comparison = cc, gene = rownames(RES[[cc]]),
  signif(RES[[cc]], 6), check.names = FALSE))), g, sep = "\t", quote = FALSE, row.names = FALSE); close(g)
## miRNA
dd <- readRDS(MI)
for (s in c("ZE","SE")) { d <- dd[, dd$Dataset == s & !dd$StageType %in% c("SE_7","SE_8")]
  d <- estimateSizeFactors(d); cnt <- counts(d); ls <- sizeFactors(d) * exp(mean(log(colSums(cnt))))
  nc <- counts(d, normalized = TRUE); cnt <- cnt[rownames(cnt) != "Unaligned", ]
  r <- fitq(cnt, d$StageType, ls, comps[[s]], "miRNA"); WT[[paste0("miRNA_", s)]] <- r$weights
  for (cc in comps[[s]]) { tt <- r$tabs[[cc]]; p <- strsplit(cc, "_vs_")[[1]]; inb <- d$StageType %in% lv(p)
    o <- data.frame(miR = rownames(cnt), baseMean = rowMeans(nc[rownames(cnt), inb]), log2FoldChange = NA_real_, lfcSE = NA_real_, stat = NA_real_, pvalue = NA_real_, padj = NA_real_)
    i <- match(rownames(tt), o$miR)
    o$log2FoldChange[i] <- tt$logFC; o$stat[i] <- tt$t; o$lfcSE[i] <- abs(tt$logFC / tt$t); o$pvalue[i] <- tt$P.Value; o$padj[i] <- tt$adj.P.Val
    o <- o[order(is.na(o$padj), o$padj), ]
    write.csv(o, file.path(OUTD, "miRNA_DE", paste0(cc, "_miRresults.csv")), row.names = FALSE) } }
saveRDS(WT, file.path(OUTD, "stage_weights.rds"))
