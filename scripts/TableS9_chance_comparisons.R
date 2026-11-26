## Supplementary Table S9 | expression support of predicted miRNA-target pairs compared with chance.
## For a set of pairs, five counts are computed:
##   opp       pairs in which the miRNA and the target were DE in opposite directions in at least one comparison
##   same      pairs in which both were DE in the same direction in at least one comparison
##   negcor    pairs with a significant negative correlation in at least one series
##   poscor    pairs with a significant positive correlation in at least one series
##   supported pairs with opposite DE or a significant negative correlation
## and two fractions: opp_frac = opp / (opp + same) and neg_frac = negcor / (negcor + poscor).
## DE and correlation follow scripts/DE_stages_voomQWB.R and scripts/correlation_support.R: miRNAs DE at BH-adjusted
## p < 0.01 and |log2 fold change| >= 0.5; Pearson correlation of vst values within each series (S1-S6; SD+Zem; SD+FG);
## a correlation counts when its p-value is at or below the largest p-value with BH-adjusted p < 0.05 in that series in
## correlation_support.R, and only when the miRNA is DE in the dataset of the series.
## Comparison pairs:
##   (1) targets shuffled among the pairs of the set (target_support_chance_shuffled_targets.tsv);
##   (2) targets replaced by genes drawn at random from the genes expressed in all three series
##       (target_support_chance_random_genes.tsv).
## Pair sets: all predicted pairs; pairs predicted by all four programs; conserved families paired with their expected
## target class (keyword rules below, applied to the eggNOG, InterPro and PANTHER descriptions of the target, taken
## from the first annotation row of each gene in the P. abies 2.0 annotation: data/target_gene_descriptions.tsv.gz);
## conserved transcription-factor modules; and single modules.
## miRNA families: the 'family' column of data/predicted_miRNA_catalogue.tsv. A miRNA with more than one family label
## (for example 'miR156|miR157') is treated as a member of each of its families.
## P = (k + 1) / (n + 1), where k is the number of the n random sets with a value at least as large as observed.
## Usage, from the deposit root: Rscript scripts/TableS9_chance_comparisons.R
## Outputs: results/target_support_chance_shuffled_targets.tsv, results/target_support_chance_random_genes.tsv,
##   results/target_support_by_family.tsv (the five counts for the pairs of each conserved family with its expected target class)
suppressPackageStartupMessages(library(data.table))

FAMILY_TABLE <- Sys.getenv("FAMILY_TABLE", "data/predicted_miRNA_catalogue.tsv")   # id, family
DESCRIPTIONS <- Sys.getenv("DESCRIPTIONS", "data/target_gene_descriptions.tsv.gz")
OUT_PREFIX   <- Sys.getenv("OUT_PREFIX", "results/target_support_")

tpF <- read.delim("data/predicted_target_pairs.tsv", stringsAsFactors = FALSE)
vmZ <- readRDS("data/miRNA_vst_ZE.rds"); vmS <- readRDS("data/miRNA_vst_SE.rds")
vgZ <- readRDS("data/mRNA_vst_ZE.rds"); vgS <- readRDS("data/mRNA_vst_SE.rds")
mdZ <- read.csv("data/metadataZE_mRNA_miRNA.csv"); mdS <- read.csv("data/metadataSE_mRNA_miRNA.csv")
cS <- sub("_S[0-9]+$", "", colnames(vgS))
mapS <- mdS[mdS$miR_User_ID %in% colnames(vmS) & mdS$SampleID %in% cS & mdS$Stage <= 6, ]
mapZ <- mdZ[mdZ$miR_User_ID %in% colnames(vmZ) & mdZ$SampleID %in% colnames(vgZ), ]
series <- list(SDZE = mapZ[mapZ$Type %in% c("SD", "ZE", "Zem"), ], SDFG = mapZ[mapZ$Type %in% c("SD", "FG"), ], SE = mapS)
M <- lapply(setNames(names(series), names(series)), function(s) { m <- series[[s]]
  if (s == "SE") list(mi = vmS[, m$miR_User_ID], g = vgS[, match(m$SampleID, cS)])
  else list(mi = vmZ[, m$miR_User_ID], g = vgZ[, m$SampleID]) })
zs <- function(X) { X <- X - rowMeans(X); X / sqrt(rowSums(X^2)) }
Z <- lapply(M, function(m) list(mi = zs(m$mi), g = zs(m$g), n = ncol(m$mi)))
pfromr <- function(r, n) { t <- r * sqrt((n - 2) / (1 - r^2)); 2 * pt(-abs(t), n - 2) }

## DE sets
comps <- setdiff(sub("_miRresults.csv", "", list.files("results/miRNA_DE", "_miRresults.csv$")), c("SE6_vs_SE7", "SE7_vs_SE8"))
mDE <- lapply(setNames(comps, comps), function(cc) {
  x <- read.csv(file.path("results/miRNA_DE", paste0(cc, "_miRresults.csv")))
  x <- x[!is.na(x$padj) & x$padj < 0.01 & abs(x$log2FoldChange) >= 0.5 & x$miR != "Unaligned", ]
  list(up = x$miR[x$log2FoldChange > 0], dn = x$miR[x$log2FoldChange < 0]) })
NM <- readRDS("results/mRNA_DE_names_voomQWB.rds")
gDE <- lapply(setNames(comps, comps), function(cc) list(up = NM$up[[cc]], dn = NM$dn[[cc]]))
deSE <- unique(unlist(mDE[grep("^SE", comps)])); deZE <- unique(unlist(mDE[grep("^SE", comps, invert = TRUE)]))
deU <- union(deSE, deZE)
keepset <- list(SE = deSE, SDZE = deZE, SDFG = deZE)

## correlation cut-off per series: the largest p-value with BH-adjusted p < 0.05 over all DE miRNA x target combinations
pcut <- sapply(names(Z), function(s) { z <- Z[[s]]
  m <- intersect(deU, rownames(z$mi)); g <- intersect(unique(tpF$target), rownames(z$g))
  P <- pfromr(z$mi[m, ] %*% t(z$g[g, ]), z$n); PA <- p.adjust(P, "BH"); max(P[!is.na(PA) & PA < 0.05]) })
cat("correlation p-value cut-offs:", paste(names(pcut), signif(pcut, 3), collapse = "; "), "\n")

paircor <- function(s, mi, g) { z <- Z[[s]]; ok <- mi %in% rownames(z$mi) & g %in% rownames(z$g)
  r <- rep(NA_real_, length(mi)); r[ok] <- rowSums(z$mi[mi[ok], , drop = FALSE] * z$g[g[ok], , drop = FALSE]); r }
stat2 <- function(pairs) {
  opp <- same <- neg <- pos <- rep(FALSE, nrow(pairs))
  for (cc in comps) { a <- mDE[[cc]]; b <- gDE[[cc]]
    mu <- pairs$miRNA %in% a$up; md <- pairs$miRNA %in% a$dn; gu <- pairs$target %in% b$up; gd <- pairs$target %in% b$dn
    opp <- opp | (mu & gd) | (md & gu); same <- same | (mu & gu) | (md & gd) }
  for (s in names(Z)) { r <- paircor(s, pairs$miRNA, pairs$target); p <- pfromr(r, Z[[s]]$n)
    k <- pairs$miRNA %in% keepset[[s]] & !is.na(r) & p <= pcut[[s]]
    neg <- neg | (k & r < 0); pos <- pos | (k & r > 0) }
  s <- c(opp = sum(opp), same = sum(same), negcor = sum(neg), poscor = sum(pos), supported = sum(opp | neg))
  c(s, opp_frac = unname(s["opp"] / max(1, s["opp"] + s["same"])), neg_frac = unname(s["negcor"] / max(1, s["negcor"] + s["poscor"]))) }

## families and expected target classes
fam <- read.delim(FAMILY_TABLE, stringsAsFactors = FALSE, colClasses = "character")
famv <- setNames(fam$family, fam[[1]]); famv[famv %in% c("", "NA")] <- NA
desc <- read.delim(gzfile(DESCRIPTIONS), quote = "", stringsAsFactors = FALSE)
blob <- setNames(tolower(paste(desc$eggnog_description, desc$eggnog_preferred_name, desc$interpro_description,
                               desc$interpro_panther_description)), desc$gene)
rules <- c(miR156 = "squamosa|\\bspl", miR159 = "myb", miR319 = "\\btcp",
  miR160 = "auxin response factor|\\barf", miR167 = "auxin response factor|\\barf",
  miR164 = "\\bnac\\b|nac domain", miR165 = "hd-zip|homeobox", miR166 = "hd-zip|homeobox",
  miR171 = "scarecrow|\\bgras|\\bscl", miR172 = "apetala|\\bap2",
  miR396 = "growth-regulating|growth regulating|\\bgrf", miR397 = "laccase",
  miR398 = "superoxide dismutase|copper chaperone", miR408 = "plantacyanin|laccase|cupredoxin|blue copper",
  miR482 = "nb-arc|lrr|leucine|resistance|nbs|nucleotide binding|nucleotide-binding",
  miR2118 = "nb-arc|lrr|leucine|resistance|nbs|nucleotide binding|nucleotide-binding",
  miR394 = "f-box", miR393 = "f-box|transport inhibitor", miR162 = "dicer|\\bdcl",
  miR168 = "argonaute|\\bago", miR399 = "ubiquitin-conjugating|pho2", miR395 = "sulfurylase|sulfate")
tfmods <- c("miR156", "miR159", "miR160", "miR164", "miR165", "miR166", "miR167", "miR171", "miR172", "miR319", "miR396")

tp <- data.frame(miRNA = tpF$miRNA, target = tpF$target, nprog = tpF$NumberOfPrograms, stringsAsFactors = FALSE)
fl <- lapply(famv[tp$miRNA], function(f) if (is.na(f)) character(0) else strsplit(f, "|", fixed = TRUE)[[1]])
## the families of a pair whose expected target class matches the target's descriptions
tp_canon_fams <- mapply(function(fs, g) { b <- blob[g]; if (is.na(b)) return(character(0))
  fs[fs %in% names(rules)][vapply(fs[fs %in% names(rules)], function(f) grepl(rules[[f]], b), logical(1))] }, fl, tp$target, SIMPLIFY = FALSE)
canon <- lengths(tp_canon_fams) > 0
inmod <- function(fams) vapply(tp_canon_fams, function(x) any(x %in% fams), logical(1))
cat("pairs: all", nrow(tp), "| four programs", sum(tp$nprog == 4), "| conserved family, expected class", sum(canon),
    "| TF modules", sum(inmod(tfmods)), "| TF modules excl. miR159", sum(inmod(setdiff(tfmods, "miR159"))), "\n")

## counts per conserved family with its expected target class
fams_all <- unique(unlist(tp_canon_fams))
byfam <- do.call(rbind, lapply(fams_all, function(f) { o <- stat2(tp[inmod(f), c("miRNA", "target")])
  data.frame(family = f, pairs = sum(inmod(f)), t(o[c("opp", "same", "negcor", "poscor", "supported")])) }))
byfam <- byfam[order(-byfam$pairs, byfam$family), ]
fwrite(byfam, paste0(OUT_PREFIX, "by_family.tsv"), sep = "\t")

summ <- function(q0, o, nl, label, extra = list()) {
  nl <- nl[complete.cases(nl), , drop = FALSE]
  d <- data.frame(subset = label, n_pairs = nrow(q0), stat = names(o), observed = round(o, 3),
             null_median = round(apply(nl, 2, median), 3),
             null_lo = round(apply(nl, 2, quantile, 0.025), 3), null_hi = round(apply(nl, 2, quantile, 0.975), 3),
             p_upper = round((colSums(nl >= matrix(o, nrow(nl), length(o), byrow = TRUE)) + 1) / (nrow(nl) + 1), 3),
             row.names = NULL)
  if (length(extra)) d <- cbind(d[1], as.data.frame(extra), d[-1])
  d }
runsub <- function(sel, label, B = 500) { q0 <- tp[sel, c("miRNA", "target")]; o <- stat2(q0)
  nl <- t(replicate(B, { q <- q0; q$target <- sample(q$target); stat2(q) })); summ(q0, o, nl, label) }
pool <- Reduce(intersect, lapply(M, function(m) rownames(m$g)[apply(m$g, 1, var) > 0]))
cat("random-gene pool", length(pool), "\n")
runrand <- function(sel, label, B = 1000) { q0 <- tp[sel, c("miRNA", "target")]; o <- stat2(q0)
  nl <- t(replicate(B, { q <- q0; q$target <- sample(pool, nrow(q), replace = TRUE); stat2(q) }))
  summ(q0, o, nl, label, list(null = "random expressed genes")) }

set.seed(7)
R <- rbind(runsub(rep(TRUE, nrow(tp)), "all predicted", 1000),
           runsub(tp$nprog == 4, "four programs"),
           runsub(canon, "conserved family, expected target class"),
           runsub(canon & !inmod(c("miR482", "miR2118", "miR397")), "conserved modules excl. miR482 and miR397"))
set.seed(11)
R <- rbind(R, runsub(inmod(tfmods), "conserved transcription-factor modules", 2000))
fwrite(R, paste0(OUT_PREFIX, "chance_shuffled_targets.tsv"), sep = "\t")

set.seed(3)
RR <- rbind(runrand(inmod("miR159"), "miR159, MYB targets"),
            runrand(inmod("miR397"), "miR397, laccase targets"),
            runrand(inmod(setdiff(tfmods, "miR159")), "TF modules excl. miR159"),
            runrand(inmod(c("miR482", "miR2118")), "miR482, NB-LRR targets"),
            runrand(!canon, "all other predicted pairs", 300))
fwrite(RR, paste0(OUT_PREFIX, "chance_random_genes.tsv"), sep = "\t")
cat("done\n")
