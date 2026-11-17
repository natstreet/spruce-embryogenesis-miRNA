## Supplementary Figure S3 | length and 5' nucleotide of the DE miRNAs, grouped as
## DE only in the zygotic series, in both, or only in the somatic series.
## DE = adjusted p < 0.01 in at least one transition (results/miRNA_DE). Supports the
## Results sentence on 21-nt length and 5' U (Figure 2A groups).
source("scripts/fig_style.R")
de <- list.files("results/miRNA_DE", full.names = TRUE)
get_de <- function(f) { d <- read.csv(f); d$miR[!is.na(d$padj) & d$padj < 0.01 & d$miR != "Unaligned"] }
deSE <- unique(unlist(lapply(de[grepl("/SE", de)], get_de)))
deZE <- unique(unlist(lapply(de[!grepl("/SE", de)], get_de)))
fa <- readLines("data/Picea_mature.fa", warn = FALSE)
h <- grep("^>", fa); ids <- sub("^>(\\S+).*", "\\1", fa[h])
seqs <- vapply(seq_along(h), function(i) paste(fa[(h[i] + 1):(if (i < length(h)) h[i + 1] - 1 else length(fa))], collapse = ""), "")
names(seqs) <- ids
stopifnot(all(c(deSE, deZE) %in% ids))
grp <- list(`DE only in the zygotic series` = setdiff(deZE, deSE),
            `DE in both series` = intersect(deZE, deSE),
            `DE only in the somatic series` = setdiff(deSE, deZE))
lev <- sprintf("%s (n = %d)", names(grp), lengths(grp))
S <- bind_rows(lapply(seq_along(grp), function(k) tibble(group = lev[k], seq = chartr("T", "U", seqs[grp[[k]]])))) %>%
  mutate(group = factor(group, levels = lev), len = nchar(seq),
         nt5 = factor(substr(seq, 1, 1), levels = c("A", "C", "G", "U")))
COL_GRP <- setNames(c(COL_SERIES[["ZE"]], "grey55", COL_SERIES[["SE"]]), lev)
lenD <- S %>% count(group, len) %>% group_by(group) %>% mutate(pct = n / sum(n)) %>% ungroup()
ntD  <- S %>% count(group, nt5, .drop = FALSE) %>% group_by(group) %>% mutate(pct = n / sum(n)) %>% ungroup()
pA <- ggplot(lenD, aes(factor(len), pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.75) +
  scale_fill_manual(values = COL_GRP, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "length (nt)", y = "miRNAs") + theme_ms()
pB <- ggplot(ntD, aes(nt5, pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = COL_GRP, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "5' nucleotide", y = "miRNAs") + theme_ms()
figS3 <- (pA + pB + plot_layout(widths = c(1.9, 1), guides = "collect")) +
  plot_annotation(tag_levels = "A") & tag_theme() & theme(legend.position = "bottom", legend.direction = "vertical")
save_fig(figS3, "FigureS3", MM_DOUBLE, 95)
S %>% group_by(group) %>% summarise(n = n(), nt21 = round(100 * mean(len == 21), 1), U5 = round(100 * mean(nt5 == "U"), 1)) %>% print()
