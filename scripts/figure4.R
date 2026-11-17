## Figure 4 | GO terms enriched among the targets of supported pairs
## Family-level enrichment is in Supplementary Table S10 and has no panel.
source("scripts/fig_style.R")

g <- read.delim("results/GO_terms_Figure4.tsv", stringsAsFactors = FALSE) %>%
  ## GO terms with BH-adjusted p < 0.05 and fewer than 5,000 annotated genes
  ## (written by scripts/GO_enrichment.R); plotted by the BH-adjusted p-value
  mutate(score = -log10(bh),
         cat = recode(cat, BP = "Biological process",
                      MF = "Molecular function", CC = "Cellular component"),
         dataset = recode(dataset, ZE = "Zygotic series", SE = "Somatic series"))
ord <- g %>% group_by(term) %>% summarise(s = max(score)) %>% arrange(s) %>% pull(term)
shared <- g %>% count(term) %>% filter(n == 2) %>% pull(term)

p4 <- g %>% mutate(term = factor(term, levels = ord)) %>%
  ggplot(aes(score, term, colour = dataset)) +
  geom_segment(aes(x = 0, xend = score, yend = term),
               position = position_dodge(width = 0.6),
               linewidth = 0.3, colour = "grey80") +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  facet_grid(rows = vars(cat), scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(`Zygotic series` = COL_SERIES[["ZE"]],
                                 `Somatic series` = COL_SERIES[["SE"]]),
                      limits = c("Zygotic series", "Somatic series"),
                      name = NULL) +
  scale_x_continuous(breaks = seq(0, 30, 10), limits = c(0, ceiling(max(g$score))),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = expression(-log[10]~"(BH-adjusted"~italic(P)*")"), y = NULL,
       title = NULL) +
  theme_ms() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.box.margin = margin(-3, 0, 0, 0),
        strip.text.y = element_text(angle = 0))
save_fig(p4, "Figure4", MM_DOUBLE, 150)
cat("terms:", n_distinct(g$term), " shared:", length(shared), "\n")
