## Figure 3 | embryo- and gametophyte-enriched miRNAs, and where they sit in
## the somatic embryo abundance distribution
source("scripts/fig_style.R")
set.seed(1)   ## ggrepel label placement is stochastic; fixed so reruns are identical
suppressPackageStartupMessages(library(ggrepel))

## the identifier is the miR column
## Classification follows the Methods: BH-adjusted p < 0.01 and |log2FC| >= 1.
## It is recomputed here from padj and log2FoldChange rather than read from the
## enriched_deseq column, so the rule in the figure cannot drift from the rule
## in the text.
d <- read.delim("results/Zem_vs_FG_enrichment.tsv", stringsAsFactors = FALSE) %>%
  mutate(id = miR,
         enriched_rule = case_when(!is.na(padj) & padj < 0.01 & log2FoldChange >= 1 ~ "Zem",
                                   !is.na(padj) & padj < 0.01 & log2FoldChange <= -1 ~ "FG",
                                   TRUE ~ "ns"))
stopifnot(identical(d$enriched_rule, d$enriched_deseq))
d <- d %>%
  mutate(cls = recode(enriched_rule, Zem = "embryo-enriched",
                      FG = "gametophyte-enriched", ns = "not enriched"),
         cls = factor(cls, levels = c("embryo-enriched", "gametophyte-enriched",
                                      "not enriched")),
         copper = family %in% c("miR397", "miR398", "miR408"))
## one class scheme for every panel, matching the Zem and FG tissue colours
## used in Figures 1 and S5; copper miRNAs are marked by outline, not colour
COL_CLS <- c(`embryo-enriched` = COL_TISSUE[["Zem"]],
             `gametophyte-enriched` = COL_TISSUE[["FG"]],
             `not enriched` = "grey70")
n_by <- table(d$cls)
cat("classes:", paste(names(n_by), as.integer(n_by), sep = "="), "\n")
cu <- d %>% filter(copper, cls == "gametophyte-enriched")
cat("gametophyte-enriched copper miRNAs:", nrow(cu), "\n")

## ---- A  enrichment against abundance rank in somatic embryos -------------
## only the miRNAs with an adjusted p-value are drawn; those removed by
## independent filtering or with all-zero counts were not tested
dA <- d %>% filter(!is.na(padj))
pA <- ggplot(dA, aes(pct_Sem, log2FoldChange)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_point(aes(colour = cls), size = 0.7, alpha = 0.85) +
  geom_point(data = cu, shape = 21, size = 2.1, stroke = 0.45,
             colour = "grey10", fill = NA) +
  geom_text_repel(data = cu, aes(label = paste0(family, " ", sub("Picea_mature_", "", id))),
                  size = SZ_TICK / .pt, colour = "grey15", min.segment.length = 0,
                  segment.size = 0.2, segment.colour = "grey55",
                  box.padding = 0.5, point.padding = 0.3, nudge_x = -0.12, nudge_y = 0.35,
                  bg.color = "white", bg.r = 0.15,
                  max.overlaps = 20) +
  scale_colour_manual(values = COL_CLS, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "abundance percentile in somatic embryos",
       y = expression("enrichment, log"[2]*" embryo / gametophyte"),
       title = NULL) +
  theme_ms() + theme(legend.position = "bottom",
                     legend.box.margin = margin(-3, 0, 0, 0))

## ---- B  where each class sits in somatic embryos ------------------------
## Sem is the result. Zem and FG are shown faded at the right because the classes
## are defined by the Zem against FG difference, so their order there is fixed by
## construction.
## The reference group is the miRNAs tested and not enriched (as in panel A).
## The miRNAs that were not tested (too rare in the paired seed libraries)
## are left out, because they are a different population: mostly rare miRNAs.
long <- d %>% filter(!is.na(padj)) %>%
  select(id, cls, Zem = pct_Zem, FG = pct_FG, Sem = pct_Sem) %>%
  pivot_longer(c(Sem, Zem, FG), names_to = "tissue", values_to = "rank") %>%
  mutate(tissue = factor(tissue, levels = c("Sem", "Zem", "FG")),
         grp = factor(ifelse(tissue == "Sem", "somatic embryos", "seed (defines classes)"),
                      levels = c("somatic embryos", "seed (defines classes)")),
         fade = tissue != "Sem",
         ckey = paste(cls, ifelse(fade, "seed", "Sem")))
COL_FADE <- c(setNames(COL_CLS, paste(names(COL_CLS), "Sem")),
              setNames(alpha(COL_CLS, 0.6), paste(names(COL_CLS), "seed")))
medB <- long %>% filter(tissue == "Sem") %>% group_by(cls) %>% summarise(m = median(rank), .groups = "drop")
mw <- wilcox.test(d$pct_Sem[d$cls == "embryo-enriched"], d$pct_Sem[d$cls == "gametophyte-enriched"])$p.value
tn <- d$cls == "not enriched" & !is.na(d$padj)
mw_fg_ns <- wilcox.test(d$pct_Sem[d$cls == "gametophyte-enriched"], d$pct_Sem[tn])$p.value
mw_em_ns <- wilcox.test(d$pct_Sem[d$cls == "embryo-enriched"], d$pct_Sem[tn])$p.value
cat("tested not enriched n =", sum(tn), "| FG vs ns p:", signif(mw_fg_ns, 2), "| embryo vs ns p:", signif(mw_em_ns, 2), "\n")
cat("Sem medians:", paste(medB$cls, round(medB$m, 3)), "| Mann-Whitney p:", signif(mw, 3), "\n")
pB <- ggplot(long, aes(tissue, rank, colour = ckey, group = interaction(tissue, cls))) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.62,
               position = position_dodge(width = 0.72)) +
  facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = COL_FADE, name = NULL, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "abundance percentile within tissue",
       title = NULL) +
  theme_ms() + theme(panel.grid.major.x = element_blank())

## ---- C  every copper-miRNA family member in the three tissues ------------
cop <- d %>% filter(copper) %>%
  select(id, family, cls, Zem = pct_Zem, FG = pct_FG, Sem = pct_Sem) %>%
  pivot_longer(c(Zem, FG, Sem), names_to = "tissue", values_to = "rank") %>%
  mutate(tissue = factor(tissue, levels = c("Zem", "FG", "Sem")),
         label = paste0(family, " ", sub("Picea_mature_", "", id)),
         flag = factor(as.character(cls), levels = levels(cls)))
pC <- ggplot(cop, aes(rank, reorder(label, rank), colour = flag)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", linewidth = 0.25,
             colour = "grey55") +
  geom_line(aes(group = label), colour = "grey80", linewidth = 0.25) +
  geom_point(aes(shape = tissue), size = 1.5) +
  scale_colour_manual(values = COL_CLS, name = NULL, drop = FALSE, guide = "none") +
  scale_shape_manual(values = c(Zem = 16, FG = 17, Sem = 15), name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "abundance percentile within tissue", y = NULL,
       title = NULL) +
  theme_ms() +
  theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey93"),
        legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"),
        legend.box.margin = margin(-3, 0, 0, 0))

fig5 <- pA + pB + pC +
  plot_layout(design = "AB\nAC", widths = c(1.25, 1), heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") & tag_theme()
save_fig(fig5, "Figure3", MM_DOUBLE, 130)
