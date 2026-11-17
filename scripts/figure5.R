## Figure 5 | somatic stages mapped onto the zygotic series
source("scripts/fig_style.R")

load_map <- function(f, lab) {
  read.delim(f, stringsAsFactors = FALSE, check.names = FALSE) %>%
    rename(se = 1) %>%
    pivot_longer(-se, names_to = "ze", values_to = "rho") %>%
    mutate(layer = lab)
}
## Profiles are centred within each series before correlating, so the panel asks
## which features are relatively high at a stage rather than which are abundant.
MIR_LAB <- "miRNA, 961 retained miRNAs"
MRNA_LAB <- "mRNA, 5,000 most variable genes"
m <- bind_rows(load_map("results/stage_mapping_miRNA.tsv", MIR_LAB),
               load_map("results/stage_mapping_mRNA.tsv", MRNA_LAB))
ZE_ORDER <- c("SD1", "SD2", "SD3", paste0("Zem", c(4:7, "8-9", 10)),
              paste0("FG", c(4:7, "8-9", 10)))
stopifnot(setequal(unique(m$ze), ZE_ORDER))
m <- m %>% mutate(ze = factor(ze, levels = ZE_ORDER),
                  se = factor(se, levels = rev(sort(unique(se)))),
                  layer = factor(layer, levels = c(MIR_LAB, MRNA_LAB)))
## best-matching zygotic embryo stage in each row
best <- m %>% filter(grepl("^Zem", ze)) %>% group_by(layer, se) %>%
  slice_max(rho, n = 1, with_ties = FALSE) %>% ungroup()
## overall row maximum across all zygotic groups (seed, embryo and gametophyte),
## marked separately because for some rows it is not a zygotic embryo stage
rowmax <- m %>% group_by(layer, se) %>% slice_max(rho, n = 1, with_ties = FALSE) %>% ungroup()
cat("row maximum per somatic stage:\n")
print(as.data.frame(rowmax %>% select(layer, se, ze, rho)), row.names = FALSE)
cat("best zygotic embryo match per somatic stage:\n")
print(as.data.frame(best %>% select(layer, se, ze, rho)), row.names = FALSE)

lim <- max(abs(m$rho))
fig6 <- ggplot(m, aes(ze, se, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = best, colour = COL_ALARM, linewidth = 0.6, fill = NA) +
  geom_point(data = rowmax, shape = 17, size = 1.1, colour = "grey10",
             position = position_nudge(y = 0.32)) +
  ## sprintf prints values in (-0.005, 0) as "-0.00"; show them as "0.00"
  geom_text(aes(label = sub("^-(0\\.00)$", "\\1", sprintf("%.2f", rho))), size = SZ_TICK / .pt,
            colour = ifelse(abs(m$rho) > 0.6 * lim, "white", "grey20")) +
  facet_wrap(~ layer, ncol = 1) +
  scale_fill_gradient2(low = COL_DIR[["down"]], mid = "white",
                       high = COL_DIR[["up"]], midpoint = 0,
                       limits = c(-lim, lim), name = "Spearman \u03c1",
                       guide = guide_colourbar(barwidth = unit(28, "mm"),
                                               barheight = unit(2.5, "mm"),
                                               title.position = "top")) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "somatic stage",
       title = NULL) +
  theme_ms() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.x.top = element_text(angle = 45, hjust = 0))
save_fig(fig6, "Figure5", MM_DOUBLE, 140)
