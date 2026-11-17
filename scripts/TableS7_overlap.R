## Supplementary Table S7: hypergeometric overlap of DE miRNAs between the two somatic culture transitions and
## each zygotic transition with at least one DE miRNA, with three background sets and BH across the 18 pairs.
suppressPackageStartupMessages({library(DESeq2)})

dd <- readRDS("data/miRNA_counts_DESeqDataSet.rds")
cn <- counts(dd); cn <- cn[rownames(cn) != "Unaligned", ]
z <- dd$Dataset == "ZE"; s <- dd$Dataset == "SE" & !dd$StageType %in% c("SE_7","SE_8"); s2 <- dd$Dataset == "SE"
for (k in list(c("SE1-6", "x"), c("SE all","y"))) { ss <- if (k[2]=="x") s else s2; cat(k[1], "expressed in both:", sum(rowSums(cn[, z]) > 0 & rowSums(cn[, ss]) > 0), "\n") }
mdir <- "results/miRNA_DE"
ld <- function(cc) { x <- read.csv(sprintf("%s/%s_miRresults.csv", mdir, cc)); k <- !is.na(x$padj) & x$padj < 0.01 & abs(x$log2FoldChange) >= 0.5
  list(up = setdiff(x$miR[k & x$log2FoldChange > 0], "Unaligned"), dn = setdiff(x$miR[k & x$log2FoldChange < 0], "Unaligned")) }
zc <- c("SD2_vs_SD1","SD3_vs_SD2","ZE4_vs_SD3","ZE5_vs_ZE4","ZE6_vs_ZE5","ZE7_vs_ZE6","ZE8-9_vs_ZE7","ZE10_vs_ZE8-9","FG4_vs_SD3","FG5_vs_FG4","FG6_vs_FG5","FG7_vs_FG6","FG8-9_vs_FG7","FG10_vs_FG8-9")
sc <- paste0("SE", 2:6, "_vs_SE", 1:5)
L <- lapply(setNames(c(zc, sc), c(zc, sc)), ld); de <- function(cc) c(L[[cc]]$up, L[[cc]]$dn)
U <- unique(unlist(lapply(L, unlist))); EB <- rownames(cn)[rowSums(cn[, z]) > 0 & rowSums(cn[, s]) > 0]
cat("union DE", length(U), " zygotic", length(unique(unlist(lapply(L[zc], unlist)))), " somatic", length(unique(unlist(lapply(L[sc], unlist)))), " both", length(intersect(unique(unlist(lapply(L[zc], unlist))), unique(unlist(lapply(L[sc], unlist))))), "\n")
lab <- function(cc) { p <- strsplit(cc, "_vs_")[[1]]; f <- function(x) sub("^ZE", "Zem", sub("^SE", "S", x)); paste0(f(p[2]), "\u2192", f(p[1])) }
zt <- zc[sapply(zc, function(x) length(de(x)) > 0)]
g <- expand.grid(ze = zt, se = c("SE3_vs_SE2","SE6_vs_SE5"), stringsAsFactors = FALSE)[, c("se","ze")]
hp <- function(a, b, N, univ) { a <- intersect(a, univ); b <- intersect(b, univ); k <- length(intersect(a, b)); phyper(k - 1, length(b), N - length(b), length(a), lower.tail = FALSE) }
o <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) { a <- de(g$se[i]); b <- de(g$ze[i]); sh <- intersect(a, b)
  sd <- length(intersect(L[[g$se[i]]]$up, L[[g$ze[i]]]$up)) + length(intersect(L[[g$se[i]]]$dn, L[[g$ze[i]]]$dn))
  data.frame(SE_transition = lab(g$se[i]), ZE_transition = lab(g$ze[i]), `DE in SE transition` = length(a), `DE in ZE transition` = length(b), overlap = length(sh),
    N_retained = 961, p_retained = hp(a, b, 961, rownames(cn)), N_expressed_both = length(EB), p_expressed_both = hp(a, b, length(EB), EB),
    N_union_DE = length(U), p_union_DE = hp(a, b, length(U), U), `same direction` = sd, `opposite direction` = length(sh) - sd,
    `shared miRNAs` = paste(sh, collapse = ";"), check.names = FALSE) }))
o$BH_retained <- p.adjust(o$p_retained, "BH"); o$BH_expressed_both <- p.adjust(o$p_expressed_both, "BH"); o$BH_union_DE <- p.adjust(o$p_union_DE, "BH")
readr::write_tsv(o, "results/TableS7_overlap.tsv")   # UTF-8 output whatever the locale (transition arrows)
print(o[, c(1:5, 7, 15, 9, 16, 11, 17, 12, 13)], digits = 3)
