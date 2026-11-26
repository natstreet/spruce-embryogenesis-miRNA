## Supplementary Figure S5 | the embryo miRNA signature across development, and
## whether the predicted targets of miRNAs share the advanced miRNA state of S5
source("scripts/fig_style.R")
suppressPackageStartupMessages(library(patchwork))

## ---- A  embryo signature index ---------------------------------------------
## Per library: mean centred vst of the Zem-enriched miRNAs minus that of the
## FG-enriched miRNAs (results/embryo_signature_index.tsv). Centring is within each dataset (stage means), so
## values compare stages within a series, not levels between series.
e6 <- read.delim("results/embryo_signature_index.tsv", stringsAsFactors = FALSE)
ZST <- c(4:7, "8-9", 10)
e6 <- e6 %>% mutate(pos = factor(ifelse(series == "zygotic", sub("^(Zem|FG)", "", stage), stage),
                                 levels = c(ZST, paste0("S", 1:6))),
                    series = factor(ifelse(series == "zygotic", "Zygotic series", "Somatic series"),
                                    levels = c("Zygotic series", "Somatic series")))
mean6 <- e6 %>% group_by(series, tissue, pos) %>% summarise(index = mean(index), .groups = "drop")
labs_z <- c("4", "5", "6", "7", "8-9", "10"); labs_s <- paste0("S", 1:6)
pA <- ggplot(e6, aes(pos, index, colour = tissue)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(data = mean6, aes(group = tissue), linewidth = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  facet_wrap(~ series, scales = "free_x") +
  scale_colour_manual(values = COL_TISSUE[c("Zem", "FG", "Sem")], breaks = c("Zem", "FG", "Sem"), name = NULL) +
  labs(x = "developmental stage", y = "embryo signature index",
       title = NULL) +
  theme_ms() + theme(legend.position = "bottom")

## ---- B  S5 mapping on miRNAs and on target gene sets, every zygotic embryo stage --
## results/S5_mapping_all_Zem_stages.tsv, written by scripts/S5_mapping_all_Zem_stages.R
ts <- read.delim("results/S5_mapping_all_Zem_stages.tsv", stringsAsFactors = FALSE, check.names = FALSE)
SET_LAB <- c(focal_miRNAs = "miRNA, all retained (961)", without_top20_contributors = "miRNA, 20 largest contributors removed (941)",
             all_expressed = "mRNA, all expressed genes (26,598)", targets_focal = "mRNA, targets of retained miRNAs (3,272)",
             targets_supported = "mRNA, targets in supported pairs (474)", targets_top10pct_contrib = "mRNA, targets of the 96 largest contributors (475)",
             targets_Zem10_DE = "mRNA, targets of miRNAs DE at Zem8-9\u2192Zem10 (216)")
ZEM <- c("Zem4", "Zem5", "Zem6", "Zem7", "Zem8-9", "Zem10")
tl <- ts %>% mutate(lab = factor(SET_LAB[set], levels = rev(SET_LAB[set]))) %>%
  select(lab, all_of(ZEM)) %>% pivot_longer(-lab, names_to = "zstage", values_to = "rho") %>%
  mutate(zstage = factor(zstage, levels = ZEM)) %>%
  group_by(lab) %>% mutate(best = rho == max(rho)) %>% ungroup()
lim <- max(abs(tl$rho))
pB <- ggplot(tl, aes(zstage, lab, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = filter(tl, best), fill = NA, colour = "black", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = COL_TISSUE[["FG"]], midpoint = 0, limits = c(-lim, lim),
                       name = "\u03c1 with S5") +
  labs(x = "zygotic embryo stage", y = NULL, title = NULL) +
  theme_ms() + theme(panel.grid = element_blank())

## ---- C  miRNAs that decrease from Zem8-9 to Zem10, in both series ------------
## The 26 miRNAs DE and lower at Zem10 than at Zem8-9. Values are vst stage means
## centred on the mean of each miRNA within its own dataset.
ZE_ORDER <- c("SD1", "SD2", "SD3", ZEM, paste0("FG", c(4:7, "8-9", 10)))
mzC <- read.csv("data/metadataZE_mRNA_miRNA.csv"); msC <- read.csv("data/metadataSE_mRNA_miRNA.csv")
mzC$grp <- ifelse(mzC$Type == "SD", paste0("SD", mzC$Stage), paste0(ifelse(mzC$Type == "ZE", "Zem", "FG"), mzC$Stage))
msC$grp <- paste0("S", msC$Stage)
vZ <- readRDS("data/miRNA_vst_ZE.rds"); vS <- readRDS("data/miRNA_vst_SE.rds")
d1 <- read.delim("results/Zem_vs_FG_enrichment.tsv"); vZ <- vZ[rownames(vZ) %in% d1$miR, ]; vS <- vS[rownames(vZ), ]
smC <- function(X, g, l) sapply(l, function(x) rowMeans(X[, g == x, drop = FALSE]))
cZC <- smC(vZ, mzC$grp[match(colnames(vZ), mzC$miR_User_ID)], ZE_ORDER); cZC <- cZC - rowMeans(cZC)
cSC <- smC(vS, msC$grp[match(colnames(vS), msC$miR_User_ID)], paste0("S", 1:6)); cSC <- cSC - rowMeans(cSC)
z10 <- read.csv("results/miRNA_DE/ZE10_vs_ZE8-9_miRresults.csv")
dn <- z10$miR[!is.na(z10$padj) & z10$padj < 0.01 & z10$log2FoldChange < 0]
cl <- rbind(data.frame(miR = rep(dn, length(ZEM)), pos = rep(ZEM, each = length(dn)), v = as.vector(cZC[dn, ZEM]), series = "Zygotic embryo"),
            data.frame(miR = rep(dn, 6), pos = rep(paste0("S", 1:6), each = length(dn)), v = as.vector(cSC[dn, ]), series = "Somatic embryo")) %>%
  mutate(series = factor(series, levels = c("Zygotic embryo", "Somatic embryo")),
         pos = factor(pos, levels = c(ZEM, paste0("S", 1:6))))
cm <- cl %>% group_by(series, pos) %>% summarise(v = median(v), .groups = "drop")
cat("Figure S5C:", length(dn), "miRNAs; S5 below SE mean:", sum(cSC[dn, "S5"] < 0), "\n")
pC <- ggplot(cl, aes(pos, v)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(group = miR), colour = "grey75", linewidth = 0.25) +
  geom_line(data = cm, aes(group = series), colour = "black", linewidth = 0.7) +
  geom_point(data = cm, size = 1.2) +
  facet_wrap(~ series, scales = "free_x") +
  labs(x = "developmental stage", y = "relative abundance (vst, centred)", title = NULL) +
  theme_ms()

figS6 <- (free(pA) / pB / free(pC)) + plot_layout(heights = c(1, 1.05, 0.9)) + plot_annotation(tag_levels = "A") & tag_theme()
cat("S6 index means:\n"); print(as.data.frame(mean6 %>% mutate(index = round(index, 2))), row.names = FALSE)
save_fig(figS6, "FigureS5", MM_DOUBLE, 205)
