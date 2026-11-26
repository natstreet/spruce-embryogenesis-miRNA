## Supplementary Table S11: for each miRNA DE in only one series, the change at the matching developmental
## phase of the other series (same limma model, eBayes without a fold-change threshold; P < 0.05, and BH < 0.05).
## Inputs: data/miRNA_counts_DESeqDataSet.rds, results/miRNA_DE/*_miRresults.csv, results/miRNA_class.tsv
## Output: results/TableS11_consistent_direction.tsv
suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR)})
lv <- function(x) sub("^([A-Z]+)(.*)$", "\\1_\\2", x)
mk <- function(x) make.names(x)
pretty <- function(z) { p <- strsplit(z, "_vs_")[[1]]; x <- paste0(p[2], "\u2192", p[1]); x <- gsub("ZE", "Zem", x); gsub("SE", "S", x) }
comps <- list(ZE = c("SD2_vs_SD1","SD3_vs_SD2","ZE4_vs_SD3","ZE5_vs_ZE4","ZE6_vs_ZE5","ZE7_vs_ZE6","ZE8-9_vs_ZE7",
                     "ZE10_vs_ZE8-9","FG4_vs_SD3","FG5_vs_FG4","FG6_vs_FG5","FG7_vs_FG6","FG8-9_vs_FG7","FG10_vs_FG8-9"),
              SE = paste0("SE", 2:6, "_vs_SE", 1:5))
## developmental phases used to match transitions between the series. Zygotic desiccation is the period from
## the last stage before desiccation (Zem6) to the desiccated seed in the cone (Zem8-9; stages 7 to 9 are the
## desiccated, dormant seed in Stojkovic et al. 2024, which used the same samples). Zem8-9 -> Zem10 (desiccated seed
## in the cone to dry seed of the previous year, from other mother trees) has no somatic counterpart and is not matched.
phase <- c(SD2_vs_SD1 = "pre", SD3_vs_SD2 = "pre", SE2_vs_SE1 = "pre",
           ZE5_vs_ZE4 = "mat", ZE6_vs_ZE5 = "mat",
           SE3_vs_SE2 = "mat", SE4_vs_SE3 = "mat", SE5_vs_SE4 = "mat",
           ZE7_vs_ZE6 = "des", "ZE8-9_vs_ZE7" = "des", SE6_vs_SE5 = "des")

## published DE calls
de <- do.call(rbind, lapply(unlist(comps), function(cc) {
  x <- read.csv(file.path("results/miRNA_DE", paste0(cc, "_miRresults.csv")))
  x <- x[x$miR != "Unaligned" & !is.na(x$padj) & x$padj < 0.01 & abs(x$log2FoldChange) >= 0.5, ]
  data.frame(miR = x$miR, comp = rep(cc, nrow(x)), dir = sign(x$log2FoldChange),
             series = rep(ifelse(startsWith(cc, "SE"), "SE", "ZE"), nrow(x))) }))
ser <- tapply(de$series, de$miR, function(s) paste(sort(unique(s)), collapse = "+"))
only <- list(ZE = names(ser)[ser == "ZE"], SE = names(ser)[ser == "SE"])

## the published model per series, tested with eBayes instead of treat
dd <- readRDS("data/miRNA_counts_DESeqDataSet.rds"); ALT <- list(); TESTED <- list()
for (s in c("ZE", "SE")) {
  d <- dd[, dd$Dataset == s & !dd$StageType %in% c("SE_7", "SE_8")]
  d <- estimateSizeFactors(d); cnt <- counts(d); ls <- sizeFactors(d) * exp(mean(log(colSums(cnt))))
  cnt <- cnt[rownames(cnt) != "Unaligned", ]
  grp <- droplevels(factor(d$StageType)); y <- DGEList(cnt, group = grp, lib.size = ls)
  keep <- filterByExpr(y, group = grp); y <- y[keep, , keep.lib.sizes = TRUE]
  levels(grp) <- mk(levels(grp)); design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  cc <- comps[[s]]
  cm <- makeContrasts(contrasts = sapply(cc, function(z) { p <- strsplit(z, "_vs_")[[1]]; paste0(mk(lv(p[1])), "-", mk(lv(p[2]))) }), levels = design)
  colnames(cm) <- cc
  v <- voomWithQualityWeights(y, design, var.group = grp)
  f <- eBayes(contrasts.fit(lmFit(v, design), cm))
  TESTED[[s]] <- rownames(y)
  for (z in cc) { tt <- topTable(f, coef = z, number = Inf, sort.by = "none")
    ALT[[z]] <- data.frame(miR = rownames(tt), comp = z, logFC = tt$logFC, P = tt$P.Value, BH = tt$adj.P.Val) } }
alt <- do.call(rbind, ALT)

## classify each miRNA
one <- function(m, own) {
  other <- setdiff(c("ZE", "SE"), own)
  ev <- de[de$miR == m & de$series == own, ]
  ev <- ev[order(match(ev$comp, comps[[own]])), ]
  de_txt <- paste(sprintf("%s %s", vapply(ev$comp, pretty, ""), ifelse(ev$dir > 0, "up", "down")), collapse = "; ")
  evm <- ev[ev$comp %in% names(phase), ]
  oc <- unique(unlist(lapply(evm$comp, function(z) names(phase)[phase == phase[[z]] & (startsWith(names(phase), "SE") == (other == "SE"))])))
  oc <- oc[order(match(oc, comps[[other]]))]
  res <- c(main = NA_character_, strict = NA_character_); oth_txt <- ""
  if (nrow(evm) == 0) res[] <- "not assessable"
  else if (!m %in% TESTED[[other]]) res[] <- "not tested"
  else {
    a <- alt[alt$miR == m & alt$comp %in% oc, ]; a <- a[order(match(a$comp, comps[[other]])), ]
    oth_txt <- paste(sprintf("%s %.2f (P = %.2g)", vapply(a$comp, pretty, ""), a$logFC, a$P), collapse = "; ")
    for (tier in c("main", "strict")) { same <- opp <- FALSE
      for (i in seq_len(nrow(evm))) {
        b <- a[a$comp %in% names(phase)[phase == phase[[evm$comp[i]]]], ]
        sig <- if (tier == "main") b$P < 0.05 else b$BH < 0.05
        same <- same || any(sig & sign(b$logFC) == evm$dir[i]); opp <- opp || any(sig & sign(b$logFC) == -evm$dir[i]) }
      res[tier] <- if (same && opp) "both directions" else if (same) "same direction" else if (opp) "opposite direction" else "no change" } }
  data.frame(miRNA = m, DE_only_in = ifelse(own == "ZE", "zygotic series", "somatic series"),
             DE_transitions = de_txt, matched_transitions_other_series = oth_txt,
             result_P_0.05 = res[["main"]], result_BH_0.05 = res[["strict"]]) }
out <- do.call(rbind, c(lapply(only$ZE, one, own = "ZE"), lapply(only$SE, one, own = "SE")))
cls <- read.delim("results/miRNA_class.tsv"); out$miRNA_class <- cls$class[match(out$miRNA, cls$miR)]
out <- out[, c("miRNA", "miRNA_class", "DE_only_in", "DE_transitions", "matched_transitions_other_series",
               "result_P_0.05", "result_BH_0.05")]
out <- out[order(out$DE_only_in != "zygotic series", as.integer(sub("Picea_mature_", "", out$miRNA))), ]
write.table(out, "results/TableS11_consistent_direction.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(table(out$result_P_0.05, out$DE_only_in)); print(table(out$result_BH_0.05, out$DE_only_in))
