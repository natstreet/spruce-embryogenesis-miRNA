## GO enrichment among the target genes of supported pairs: topGO parent-child algorithm, Fisher's exact
## test, terms kept at BH <= 0.05. Series level (Supplementary Tables S4 and S5, Figure 4) and per miRNA
## family and arm (Supplementary Table S10).
## Background genes: vst > 0 in at least two libraries of data/mRNA_vst_SE.rds and data/mRNA_vst_ZE.rds.
## Usage, from the deposit root: Rscript scripts/GO_enrichment.R
## Outputs: results/GO/GO_terms.tsv, results/GO/GO_results.rds, results/GO_terms_Figure4.tsv
suppressPackageStartupMessages({library(tidyverse); library(matrixStats); library(topGO); library(GO.db)})
invisible(Sys.setlocale("LC_COLLATE", "C"))   # locale-independent gene order within each term
out <- "results/GO"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
cat("topGO", as.character(packageVersion("topGO")), "GO.db", as.character(packageVersion("GO.db")), "\n")

goannot <- read_tsv("data/GO_annotation_gopher_20261002.tsv.gz", col_names = FALSE, show_col_types = FALSE) %>%
  mutate(X2 = strsplit(X2, "\\|")) %>% split(., .$X1) %>% map(pull, X2) %>% map(unlist)

vstmSE <- readRDS("data/mRNA_vst_SE.rds")
vstmZE <- readRDS("data/mRNA_vst_ZE.rds")
fl <- function(m) m[rowSums(m) > 0, ] %>% .[rowSds(.) > 0, ]
bSE <- fl(vstmSE); bZE <- fl(vstmZE)
backgroundSE <- rownames(bSE)[rowSums(bSE > 0) > 1]
backgroundZE <- rownames(bZE)[rowSums(bZE > 0) > 1]
rname <- intersect(rownames(vstmSE), rownames(vstmZE))
vstm <- fl(cbind(vstmSE[rname, ], vstmZE[rname, ]))
background <- rownames(vstm)[rowSums(vstm > 0) > 1]
cat("background SE", length(backgroundSE), "ZE", length(backgroundZE), "combined", length(background), "\n")

ap <- read_tsv("data/predicted_miRNA_catalogue.tsv", show_col_types = FALSE, col_types = cols(.default = "c"))
fam <- ap %>% transmute(miRNA = id, family = na_if(family, ""), arm = na_if(arm, ""))
key <- function(f, ar) gsub("-NA$", "", paste0(ifelse(is.na(f), "novel", f), "-", ar))

  q <- read_tsv("results/target_pair_support_flags.tsv", show_col_types = FALSE)   # supported pairs
resSE <- q %>% filter(isSE) %>% pull(target) %>% unique()
resZE <- q %>% filter(isZE) %>% pull(target) %>% unique()
sup <- q %>% left_join(fam, by = "miRNA") %>% mutate(k = key(family, arm)) %>% dplyr::select(miRNA, target, k)
q <- read_tsv("results/target_pair_correlations.tsv", show_col_types = FALSE) %>% dplyr::select(miRNA, target)  # all predicted pairs
de <- list.files("results/miRNA_DE", "_miRresults.csv$", full.names = TRUE) %>%
  map(read_csv, show_col_types = FALSE) %>% map(~ .x$miR[!is.na(.x$padj) & .x$padj < 0.01]) %>% unlist() %>% unique()
pred <- q %>% filter(miRNA %in% de) %>% left_join(fam, by = "miRNA") %>% filter(!is.na(family)) %>%
  transmute(miRNA, target, k = key(family, arm))
cat("set SE", length(resSE), "ZE", length(resZE), "supported pairs", nrow(sup), "predicted pairs (DE, annotated)", nrow(pred), "\n")

## Same steps as the topGO() helper of the UPSC bioinformatics facility (terms returned = those with BH <= 0.05; its 'FDR'
## column is GenTable's raw p-value), with the BH value of each returned term added.
topGO_bh <- function(set, bg, annotation, algorithm, alpha = 0.05, ontology = c("BP", "CC", "MF")) {
  if (!any(bg %in% set)) return(setNames(vector("list", length(ontology)), ontology))
  allGenes <- factor(as.integer(bg %in% set)); names(allGenes) <- bg
  lst <- map(ontology, function(o) {
    GOdata <- new("topGOdata", ontology = o, allGenes = allGenes, annot = annFUN.gene2GO, gene2GO = annotation)
    allGO <- genesInTerm(GOdata)
    results <- runTest(GOdata, algorithm = algorithm, statistic = "fisher")
    bh <- p.adjust(score(results), method = "BH")
    n <- sum(bh <= alpha)
    if (n == 0) return(NULL)
    as_tibble(GenTable(GOdata, results, topNodes = n, numChar = 1000)) %>%
      rename_with(function(sel) "FDR", .cols = last_col()) %>%
      mutate(P = score(results)[GO.ID], BH = bh[GO.ID],
             siggenes = map_chr(GO.ID, ~ paste(intersect(allGO[[.x]], set), collapse = "|")))
  })
  names(lst) <- ontology; lst
}
run <- function(s, bg, alg) suppressMessages(topGO_bh(s, bg, goannot, algorithm = alg))
ds <- list(Targets_SE = run(resSE, backgroundSE, "parentchild"),
           Targets_ZE = run(resZE, backgroundZE, "parentchild"))
pmap1 <- function(L) parallel::mclapply(L, run, background, "parentchild", mc.cores = 12)
fw <- sup %>% split(.$k) %>% map(~ unique(.x$target)) %>% pmap1()
pw <- pred %>% split(.$k) %>% map(~ unique(.x$target)) %>% pmap1()
stopifnot(!any(map_lgl(c(fw, pw), ~ inherits(.x, "try-error"))))

## BH over all terms scored, kept alongside the raw P that the helper labels FDR
saveRDS(list(dataset = ds, family = fw, predicted = pw,
             sets = list(SE = resSE, ZE = resZE, family = map(split(sup, sup$k), ~ unique(.x$target)),
                         predicted = map(split(pred, pred$k), ~ unique(.x$target))),
             background = list(SE = backgroundSE, ZE = backgroundZE, combined = background)),
        file.path(out, "GO_results.rds"))
flat <- function(L) imap_dfr(L, function(e, nm) imap_dfr(compact(e), ~ as_tibble(.x) %>% mutate(set = nm, ont = .y, .before = 1)))
write_tsv(bind_rows(flat(ds) %>% mutate(level = "dataset"), flat(fw) %>% mutate(level = "family"),
                    flat(pw) %>% mutate(level = "predicted")) %>% dplyr::select(-any_of("allgenes")),
          file.path(out, "GO_terms.tsv"))
## Figure 4: series-level terms with fewer than 5,000 annotated genes
f4 <- flat(ds) %>% filter(as.numeric(Annotated) < 5000) %>%
  transmute(cat = ont, id = GO.ID, term = Term, bh = BH, dataset = sub("^Targets_", "", set))
write_tsv(f4, "results/GO_terms_Figure4.tsv")
cat("done\n")
