## Supplementary Figure S7 | miRNA-mRNA target pairs with supporting evidence
source("scripts/fig_style.R")
suppressPackageStartupMessages({
  library(ggVennDiagram); library(igraph); library(ggraph); library(tidygraph)
  library(ggrepel)
})
## Negative correlation is assessed within each dataset (SE libraries, zygotic
## embryo libraries, female gametophyte libraries), not across pooled libraries,
## which would mix genotypes and sequencing batches.
negc <- function(r, p) !is.na(r) & r < 0 & !is.na(p) & p < 0.05
## Opposite-DE comparisons per pair (miRNA and mRNA both DE, BH < 0.01 and
## abs(lfc) >= 0.5, in opposite directions in the same comparison) and the
## isSE / isZE flags are read from results/target_pair_support_flags.tsv.
flags <- read_tsv("results/target_pair_support_flags.tsv",
                  show_col_types = FALSE) %>%
  transmute(miRNA, target, selComp_new = opposite_DE_comparisons,
            isSE_new = isSE, isZE_new = isZE)
## Within-dataset correlations (scripts/correlation_support.R): BH over all DE miRNAs
## x all predicted target genes in each series, kept only where the miRNA is DE in that
## series.
corr <- read_tsv("results/target_pair_correlations.tsv", show_col_types = FALSE)
tp <- read_tsv("data/predicted_target_pairs.tsv", show_col_types = FALSE) %>%
  select(-any_of(c("corSE","padjSE","corSDZE","padjSDZE","corSDFG","padjSDFG"))) %>%
  left_join(corr, by = c("miRNA", "target")) %>%
  left_join(flags, by = c("miRNA", "target")) %>%
  mutate(selComp = selComp_new,
         isSE = !is.na(isSE_new) & isSE_new,
         isZE = !is.na(isZE_new) & isZE_new) %>%
  select(-selComp_new, -isSE_new, -isZE_new) %>%
  mutate(oppDE = !is.na(selComp),
         negCorr = negc(corSE, padjSE) | negc(corSDZE, padjSDZE) |
                   negc(corSDFG, padjSDFG))
stopifnot(!any(tp$negCorr & !(tp$isSE | tp$isZE)))
ann <- read_tsv("results/supported_target_pairs_annotated.tsv", show_col_types = FALSE) %>%
  select(miRNA, target, family, arm, annotation) %>% distinct(miRNA, target, .keep_all = TRUE)
## family assignment as in Supplementary Table S8 (Picea_mature_1209 has no family)
famS8 <- read_tsv("results/supported_pair_family_TableS8.tsv", show_col_types = FALSE)
ann <- ann %>% select(-family) %>% left_join(famS8, by = c("miRNA", "target"))

## miRNA direction per transition, in developmental direction
de <- list.files("results/miRNA_DE", pattern = "_miRresults.csv$", full.names = TRUE)
names(de) <- sub("_miRresults.csv$", "", basename(de))
## the comparisons are named later_vs_earlier, so no sign change is applied
mdir <- imap_dfr(de, ~ read.csv(.x) %>% mutate(comp = .y)) %>%
  filter(miR != "Unaligned") %>%
  transmute(miRNA = miR, comp, mlfc = log2FoldChange)

pretty_comp <- function(x) {
  p <- strsplit(as.character(x), "_vs_")
  x <- vapply(p, function(z) paste0(z[2], "\u2192", z[1]), character(1))
  gsub("SE", "S", gsub("ZE", "Zem", x))
}
COMP_ORDER <- c("SD2_vs_SD1", "SD3_vs_SD2", "ZE4_vs_SD3", "ZE5_vs_ZE4",
                "ZE6_vs_ZE5", "ZE7_vs_ZE6", "ZE8-9_vs_ZE7", "ZE10_vs_ZE8-9",
                "FG4_vs_SD3", "FG5_vs_FG4", "FG6_vs_FG5", "FG7_vs_FG6",
                "FG8-9_vs_FG7", "FG10_vs_FG8-9",
                paste0("SE", 2:6, "_vs_SE", 1:5))

## ---- A  four-way overlap of the evidence types ---------------------------
f3a <- list(
  "opposite DE"          = which(tp$oppDE),
  "negative correlation" = which(tp$negCorr),
  "in SE"                = which(tp$isSE),
  "in ZE"                = which(tp$isZE))
cat("Fig6A:", paste(names(f3a), lengths(f3a), sep = "="), "\n")
pA <- ggVennDiagram(f3a, label = "count", label_alpha = 0,
                    label_size = SZ_TICK / .pt, set_size = SZ_ANN / .pt,
                    edge_size = 0.3) +
  scale_fill_gradient(low = "white", high = "#9ecae1", guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = NULL) +
  theme_void_ms()

## ---- B, C  opposite-DE pairs per transition, by miRNA direction ----------
pairs_long <- tp %>% filter(oppDE) %>%
  separate_longer_delim(selComp, delim = ",") %>%
  mutate(comp = trimws(selComp)) %>%
  left_join(mdir, by = c("miRNA", "comp")) %>%
  mutate(mdirection = ifelse(mlfc > 0, "miRNA higher in later stage",
                             "miRNA lower in later stage"))
stopifnot(!any(is.na(pairs_long$mlfc)))
pair_panel <- function(cs, title) {
  d <- pairs_long %>% filter(comp %in% cs) %>% count(comp, mdirection) %>%
    mutate(comp = factor(pretty_comp(comp),
                         levels = pretty_comp(intersect(COMP_ORDER, cs))),
           n = ifelse(grepl("lower", mdirection), -n, n))
  ggplot(d, aes(comp, n, fill = mdirection)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
    geom_text(aes(label = abs(n), vjust = ifelse(n > 0, -0.35, 1.25)),
              size = SZ_TICK / .pt, colour = "grey20") +
    scale_fill_manual(values = setNames(c(COL_DIR[["up"]], COL_DIR[["down"]]),
                                        c("miRNA higher in later stage",
                                          "miRNA lower in later stage")),
                      name = NULL) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.18, 0.18))) +
    labs(x = NULL, y = "target pairs", title = title) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank())
}
cat("FigS7B/C counts (comparison, miRNA direction, pairs):\n")
print(as.data.frame(pairs_long %>% count(comp, mdirection)), row.names = FALSE)
zecomp <- grep("^SD|^ZE", COMP_ORDER, value = TRUE)
secomp <- grep("^SE", COMP_ORDER, value = TRUE)
pB <- pair_panel(zecomp, "Zygotic series") +
  theme(legend.position = "bottom", legend.box.margin = margin(-3, 0, 0, 0))
pC <- pair_panel(secomp, "Somatic series") + theme(legend.position = "none")

## ---- D  family-to-target network for the annotated families -------------
ev <- tp %>% filter(oppDE | negCorr | isSE | isZE) %>%
  left_join(ann, by = c("miRNA", "target")) %>% filter(!is.na(family)) %>%
  mutate(tissue = case_when(isSE & isZE ~ "both series", isSE ~ "SE only",
                            isZE ~ "ZE only", TRUE ~ "correlation only"),
         edge = case_when(negCorr & oppDE ~ "both",
                          negCorr ~ "negative correlation", TRUE ~ "opposite DE"))
cat("Fig6D: families", n_distinct(ev$family), " targets", n_distinct(ev$target),
    " pairs", nrow(ev), "\n")
edges <- ev %>% transmute(from = family, to = target, edge)
nodes <- bind_rows(
  tibble(name = unique(ev$family), kind = "miRNA family", tissue = NA_character_),
  ev %>% group_by(name = target) %>%
    summarise(kind = "mRNA target",
              tissue = names(sort(table(tissue), decreasing = TRUE))[1],
              .groups = "drop"))
g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
set.seed(11)
pD <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(linetype = edge), colour = "grey70", edge_width = 0.22) +
  ## two layers: the target nodes carry the tissue colour; the family nodes have
  ## no tissue (NA) and are drawn without a colour mapping
  geom_node_point(aes(filter = kind == "mRNA target", shape = kind,
                      colour = tissue, size = kind)) +
  geom_node_point(aes(filter = kind == "miRNA family", shape = kind,
                      size = kind), colour = "grey20") +
  geom_node_text(aes(label = ifelse(kind == "miRNA family", name, NA)),
                 size = SZ_TICK / .pt, repel = TRUE, max.overlaps = 40,
                 min.segment.length = 0.1, segment.size = 0.15,
                 segment.colour = "grey60", colour = "grey10", bg.color = "white",
                 bg.r = 0.08) +
  scale_shape_manual(values = c(`miRNA family` = 16, `mRNA target` = 15),
                     name = NULL) +
  scale_size_manual(values = c(`miRNA family` = 1.9, `mRNA target` = 0.8),
                    guide = "none") +
  scale_colour_manual(values = c(`SE only` = COL_TISSUE[["Sem"]],
                                 `ZE only` = COL_TISSUE[["Zem"]],
                                 `both series` = "grey45",
                                 `correlation only` = COL_TISSUE[["SDL"]]),
                      na.value = "grey20", name = NULL, na.translate = FALSE) +
  scale_edge_linetype_manual(values = c(both = "solid",
                                        `negative correlation` = "dashed",
                                        `opposite DE` = "dotted"), name = NULL) +
  ## legend order: tissue colour, node type, edge type
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2),
         edge_linetype = guide_legend(order = 3)) +
  labs(title = NULL) +
  theme_void_ms() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"))

## ---- E  direction and evidence by family and transition -----------------
## The fold change is averaged over the distinct miRNA members of a family, not
## over pairs: averaging over pairs weights each member by its number of
## predicted targets. Cells where members change in opposite directions are
## outlined, because the mean alone hides the disagreement.
famev <- pairs_long %>%
  left_join(ann, by = c("miRNA", "target")) %>% filter(!is.na(family))
fam_memb <- famev %>% distinct(family, comp, miRNA, mlfc) %>%
  group_by(family, comp) %>%
  summarise(mixed = any(mlfc > 0) & any(mlfc < 0),
            nmemb = n(), mlfc = mean(mlfc), .groups = "drop")
fam_tr <- famev %>% group_by(family, comp) %>%
  summarise(n = n(), anycor = any(negCorr), .groups = "drop") %>%
  left_join(fam_memb, by = c("family", "comp"))
stopifnot(!any(is.na(fam_tr$mlfc)))
cat("Fig6E: cells", nrow(fam_tr), " with members disagreeing in sign",
    sum(fam_tr$mixed), "\n")
fam_tr <- fam_tr %>%
  mutate(series = factor(ifelse(grepl("^SE", comp), "Somatic series", "Zygotic series"),
                         levels = c("Zygotic series", "Somatic series")),
         comp = factor(pretty_comp(comp),
                       levels = pretty_comp(intersect(COMP_ORDER, unique(comp)))))
famord <- fam_tr %>% group_by(family) %>% summarise(n = sum(n)) %>%
  arrange(n) %>% pull(family)
lim <- max(abs(fam_tr$mlfc))
pE <- fam_tr %>% mutate(family = factor(family, levels = famord)) %>%
  ggplot(aes(comp, family, fill = mlfc)) +
  geom_tile(aes(colour = mixed, linewidth = mixed)) +
  scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "grey10"),
                      guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.5),
                         guide = "none") +
  geom_point(data = ~ subset(.x, anycor), size = 0.6, colour = "grey15") +
  scale_fill_gradient2(low = COL_DIR[["down"]], mid = "white",
                       high = COL_DIR[["up"]], midpoint = 0,
                       limits = c(-lim, lim),
                       name = expression("mean miRNA log"[2]*" fold change"),
                       guide = guide_colourbar(barwidth = unit(22, "mm"),
                                               barheight = unit(2.5, "mm"),
                                               title.position = "top")) +
  facet_grid(. ~ series, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = NULL) +
  theme_ms() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = "bottom",
        panel.spacing.x = unit(2, "mm"))

fig3 <- pA + pB + pC + pD + pE +
  plot_layout(design = "AABC\nDDEE\nDDEE", heights = c(1, 1.1, 1.1)) +
  plot_annotation(tag_levels = "A") & tag_theme()
save_fig(fig3, "FigureS7", MM_DOUBLE, 210)
