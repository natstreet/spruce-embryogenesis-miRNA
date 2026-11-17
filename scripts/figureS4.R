## Supplementary Figure S4 | members of each miRNA family among the 961 retained
## miRNAs, split by the series in which each member was detected.
## Detected = at least one read in at least one library of the group.
## The somatic split uses the somatic embryo libraries S1-S6 only (f1b$Sem); the
## seedling libraries S7-S8 are excluded, as everywhere else in the paper.
## Members detected only in seedlings are shown as their own group.
source("scripts/fig_style.R")

dat <- readRDS("results/figure1_panel_data.rds")
foc <- dat$focal
ZEs <- dat$f1a$ZE          # any read in any of the 51 ZE libraries (SD, FG, Zem)
SEs <- dat$f1b$Sem         # any read in any of the 17 somatic embryo libraries S1-S6
SDLs <- dat$f1b$SDL        # seedling libraries S7-S8, used only to name the last group
stopifnot(setequal(intersect(ZEs, foc), Reduce(union, dat$f1b[c("SD", "FG", "Zem")])))

D1 <- read.delim("results/Zem_vs_FG_enrichment.tsv", stringsAsFactors = FALSE)
famv <- setNames(D1$family, D1$miR)
stopifnot(all(foc %in% D1$miR))

GRP <- c("ZE only", "SE only (S1-S6)", "both", "seedlings only")
grp <- ifelse(foc %in% ZEs & foc %in% SEs, "both",
       ifelse(foc %in% ZEs, "ZE only",
       ifelse(foc %in% SEs, "SE only (S1-S6)", "seedlings only")))
stopifnot(all(foc[grp == "seedlings only"] %in% SDLs))
cat("retained miRNAs by group:", paste(names(table(grp)), table(grp), collapse = ", "), "\n")

fd <- data.frame(miR = foc, family = famv[foc], grp = grp, stringsAsFactors = FALSE)
fd <- fd[!is.na(fd$family), ]
nf <- as.data.frame(table(fd$family), stringsAsFactors = FALSE); names(nf) <- c("family", "n")
keepf <- nf$family[nf$n >= 3]
cat("families with >= 3 members:", length(keepf), "| retained miRNAs with no family:",
    sum(is.na(famv[foc])), "\n")
pd <- fd %>% filter(family %in% keepf) %>% count(family, grp, name = "n")
ord <- nf$family[nf$family %in% keepf][order(-nf$n[nf$family %in% keepf])]
pd <- pd %>% mutate(family = factor(family, levels = rev(ord)),
                    ## stacking order: "both" starts at zero
                    grp = factor(grp, levels = c("seedlings only", "ZE only", "SE only (S1-S6)", "both")))
cat("plotted members by group:", paste(GRP, tapply(pd$n, pd$grp, sum)[GRP], collapse = ", "),
    "| total", sum(pd$n), "\n")
write.table(pd %>% arrange(desc(family), grp), "results/FigureS4_family_counts.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

COL_S4 <- setNames(c(COL_SERIES[["ZE"]], COL_SERIES[["SE"]], "grey60", COL_TISSUE[["SDL"]]), GRP)
gS4 <- ggplot(pd, aes(n, family, fill = grp)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = COL_S4, breaks = GRP, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Retained miRNAs in family", y = NULL) +
  theme_ms() +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
save_fig(gS4, "FigureS4", MM_DOUBLE, max(90, 5 + 3.2 * length(keepf)))
