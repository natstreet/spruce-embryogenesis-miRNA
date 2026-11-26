## Supplementary Figures S1 and S6
source("scripts/fig_style.R")

## ---- S1  length and 5' nucleotide of the replicate-checked sequences -------
## Grouped as the text describes them: sequences meeting the plant criteria,
## database-matched sequences that do not, and novel sequences that do not.
sq <- read.delim("data/replicate_checked_sequences.tsv", stringsAsFactors = FALSE)
GRP_LAB <- c(plant = "meet plant criteria", db_not_plant = "database match, criteria not met",
             novel_not_plant = "novel, criteria not met")
sq <- sq %>% mutate(group = factor(GRP_LAB[group], levels = GRP_LAB),
                    len = nchar(sequence),
                    nt5 = factor(chartr("T", "U", substr(sequence, 1, 1)), levels = c("A", "C", "G", "U")))
n_grp <- table(sq$group)
lev_n <- sprintf("%s (n = %s)", names(n_grp), format(as.integer(n_grp), big.mark = ",", trim = TRUE))
sq$group_n <- factor(lev_n[match(sq$group, names(n_grp))], levels = lev_n)
COL_GRP <- setNames(c(COL_TISSUE[["Zem"]], COL_TISSUE[["SDL"]], "grey55"), lev_n)

lenD <- sq %>% count(group_n, len) %>% group_by(group_n) %>% mutate(pct = n / sum(n)) %>% ungroup()
pA <- ggplot(lenD, aes(factor(len), pct, fill = group_n)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.75) +
  scale_fill_manual(values = COL_GRP, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "length (nt)", y = "sequences") +
  theme_ms() + theme(legend.position = "bottom", legend.direction = "vertical")
ntD <- sq %>% count(group_n, nt5) %>% group_by(group_n) %>% mutate(pct = n / sum(n)) %>% ungroup()
pB <- ggplot(ntD, aes(nt5, pct, fill = group_n)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = COL_GRP, name = NULL, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "5' nucleotide", y = "sequences") +
  theme_ms()
u5 <- ntD %>% filter(nt5 == "U")
cat("5'U by group:", paste(u5$group_n, sprintf("%.0f%%", 100 * u5$pct), collapse = "; "), "\n")
figS1 <- pA + pB + plot_layout(widths = c(2.2, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & tag_theme() & theme(legend.position = "bottom")
save_fig(figS1, "FigureS1", MM_DOUBLE, 95)

## ---- S6  the seven miRNAs shared by S2->S3 and Zem8-9->Zem10 --------------
s5 <- read.delim("results/FigureS6_shared_novel_miRNAs.tsv", stringsAsFactors = FALSE)
COMP_LAB <- c(SE3_vs_SE2 = "S2\u2192S3\ntransfer to maturation medium",
              `ZE10_vs_ZE8-9` = "Zem8-9\u2192Zem10\nto dry seed",
              SE6_vs_SE5 = "S5\u2192S6\nsomatic desiccation")
s5 <- s5 %>% mutate(comp = factor(COMP_LAB[comp], levels = COMP_LAB),
                    id = sub("Picea_mature_", "", miR),
                    id = factor(id, levels = rev(unique(id[order(log2FoldChange[comp == COMP_LAB[1]])]))))
lim <- max(abs(s5$log2FoldChange))
figS5 <- ggplot(s5, aes(comp, id, fill = log2FoldChange)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%s", log2FoldChange, ifelse(DE, "*", ""))),
            size = SZ_TICK / .pt, colour = "grey10") +
  scale_fill_gradient2(low = COL_DIR[["down"]], mid = "white", high = COL_DIR[["up"]],
                       midpoint = 0, limits = c(-lim, lim),
                       name = expression("log"[2]*" fold change")) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "novel miRNA (Picea_mature identifier)",
       title = NULL) +
  theme_ms() + theme(panel.grid = element_blank(), legend.position = "right")
save_fig(figS5, "FigureS6", MM_DOUBLE * 0.75, 75)
