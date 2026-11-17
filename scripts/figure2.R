## Figure 2 | differential expression between adjacent developmental stages
## Sign convention: the contrasts are named later_vs_earlier, so a positive
## log2 fold change means higher in the later stage. No sign is changed.
source("scripts/fig_style.R")
suppressPackageStartupMessages({ library(ggVennDiagram) })

de <- list.files("results/miRNA_DE", pattern = "_miRresults.csv$",
                 full.names = TRUE)
names(de) <- sub("_miRresults.csv$", "", basename(de))
res <- imap_dfr(de, ~ read.csv(.x) %>% mutate(comp = .y)) %>%
  filter(miR != "Unaligned") %>%
  mutate(lfc = log2FoldChange,
         sig = !is.na(padj) & padj < 0.01 & abs(lfc) >= 0.5,
         later = ifelse(lfc > 0, "higher in later stage", "lower in later stage"))

COMP_ORDER <- c("SD2_vs_SD1", "SD3_vs_SD2", "ZE4_vs_SD3", "ZE5_vs_ZE4",
                "ZE6_vs_ZE5", "ZE7_vs_ZE6", "ZE8-9_vs_ZE7", "ZE10_vs_ZE8-9",
                "FG4_vs_SD3", "FG5_vs_FG4", "FG6_vs_FG5", "FG7_vs_FG6",
                "FG8-9_vs_FG7", "FG10_vs_FG8-9",
                paste0("SE", 2:6, "_vs_SE", 1:5))
## the names read later_vs_earlier, so the arrow is drawn from the second part
pretty_comp <- function(x) {
  p <- strsplit(as.character(x), "_vs_")
  x <- vapply(p, function(z) paste0(z[2], "\u2192", z[1]), character(1))
  x <- gsub("ZE", "Zem", x); x <- gsub("SE", "S", x); x
}
stopifnot(all(res$comp %in% COMP_ORDER))
sigd <- res %>% filter(sig)
cat("DE miRNAs:", n_distinct(sigd$miR), "\n")

## ---- A  overlap of the DE sets from the four tissues ----------------------
setof <- function(cs) unique(sigd$miR[sigd$comp %in% cs])
f2a <- list(
  SD  = setof(c("SD2_vs_SD1", "SD3_vs_SD2")),
  FG  = setof(grep("^FG", COMP_ORDER, value = TRUE)),
  Zem = setof(grep("^ZE", COMP_ORDER, value = TRUE)),
  Sem = setof(grep("^SE", COMP_ORDER, value = TRUE)))
cat("Fig2A sets:", paste(names(f2a), lengths(f2a), sep = "="), "\n")
pA <- ggVennDiagram(f2a, label = "count", label_alpha = 0,
                    label_size = SZ_TICK / .pt, set_size = SZ_ANN / .pt,
                    edge_size = 0.3) +
  scale_fill_gradient(low = "white", high = "#9ecae1", guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(title = NULL) +
  theme_void_ms()

## ---- B, C  DE counts per transition --------------------------------------
bar_panel <- function(cs, title) {
  d <- sigd %>% filter(comp %in% cs) %>%
    count(comp, later) %>%
    mutate(comp = factor(pretty_comp(comp),
                         levels = pretty_comp(intersect(COMP_ORDER, cs))),
           n = ifelse(later == "lower in later stage", -n, n))
  ggplot(d, aes(comp, n, fill = later)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
    geom_text(aes(label = abs(n), vjust = ifelse(n > 0, -0.35, 1.25)),
              size = SZ_TICK / .pt, colour = "grey20") +
    scale_fill_manual(values = c(`higher in later stage` = COL_DIR[["up"]],
                                 `lower in later stage` = COL_DIR[["down"]]),
                      name = NULL) +
    scale_y_continuous(labels = abs,
                       expand = expansion(mult = c(0.16, 0.16))) +
    labs(x = NULL, y = "DE miRNAs", title = title) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank())
}
pB <- bar_panel(grep("^SD|^ZE|^FG", COMP_ORDER, value = TRUE),
                "Zygotic series")
pC <- bar_panel(grep("^SE", COMP_ORDER, value = TRUE),
                "Somatic series")

## ---- D  fold change of every DE miRNA in every transition ----------------
keep <- sigd %>% count(comp) %>% filter(n > 0) %>% pull(comp)
hm <- res %>% filter(miR %in% unique(sigd$miR), comp %in% keep) %>%
  mutate(value = ifelse(sig, lfc, NA_real_))
wide <- hm %>% select(miR, comp, value) %>%
  pivot_wider(names_from = comp, values_from = value) %>% column_to_rownames("miR")
ord <- hclust(dist(replace(as.matrix(wide), is.na(as.matrix(wide)), 0)))$order
lim <- max(abs(hm$value), na.rm = TRUE)
pD <- hm %>%
  ## series is read from the raw comparison name, before relabelling: the label
  ## helper is not idempotent, so applying it to its own output loses the series
  mutate(series = factor(ifelse(grepl("^SE", comp),
                                "Somatic series", "Zygotic series"),
                         levels = c("Zygotic series", "Somatic series")),
         miR = factor(miR, levels = rownames(wide)[ord]),
         comp = factor(pretty_comp(comp),
                       levels = pretty_comp(intersect(COMP_ORDER, keep)))) %>%
  ggplot(aes(comp, miR, fill = value)) +
  geom_tile() +
  facet_grid(cols = vars(series), scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = COL_DIR[["down"]], mid = "white",
                       high = COL_DIR[["up"]], midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey85",
                       name = expression(log[2]~"fold change"),
                       guide = guide_colourbar(barwidth = unit(24, "mm"),
                                               barheight = unit(2.5, "mm"),
                                               title.position = "top")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = sprintf("%d DE miRNAs", n_distinct(sigd$miR)),
       title = NULL) +
  theme_ms() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = "bottom")

pB <- pB + theme(legend.position = "bottom",
                 legend.box.margin = margin(-3, 0, 0, 0))
pC <- pC + theme(legend.position = "none")
fig2 <- pA + pB + pC + pD +
  plot_layout(design = "AABB\nAACC\nDDDD\nDDDD",
              heights = c(1, 1, 1.15, 1.15)) +
  plot_annotation(tag_levels = "A") & tag_theme()
save_fig(fig2, "Figure2", MM_DOUBLE, 185)
