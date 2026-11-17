## Figure 1 | identification, detection and global expression of the retained miRNAs
source("scripts/fig_style.R")
suppressPackageStartupMessages({ library(ggVennDiagram); library(ggrepel) })

dat <- readRDS("results/figure1_panel_data.rds")
f1a <- dat$f1a; f1b <- dat$f1b; focal <- dat$focal; vall <- dat$vall
## known / variant / novel, from the BLAST identity of each predicted sequence
cls <- read.delim("results/miRNA_class.tsv", stringsAsFactors = FALSE)
CLS_ORDER <- c("known", "variant", "novel")
COL_CLASS <- c(known = "#c6dbef", variant = "#6baed6", novel = "#2171b5")
miR_class <- setNames(factor(cls$class, levels = CLS_ORDER), cls$miR)
stopifnot(!any(is.na(miR_class[focal])))
## the metadata calls the zygotic embryo tissue "ZE"; the manuscript calls it Zem
grp <- dat$grp %>% mutate(tissue = ifelse(tissue == "ZE", "Zem", tissue))
names(f1b)[names(f1b) == "ZE"] <- "Zem"
STAGE_ORDER <- c(paste0("SD", 1:3), paste0("Zem", c(4:7, "8-9", 10)),
                 paste0("FG", c(4:7, "8-9", 10)), paste0("S", 1:8))

## ---- A  four-set Venn over the 2,219 replicate-checked sequences -----------
## Over the sequences retained by the replicate check, not over all predictions:
## a sequence detected in a single replicate cannot be assessed in any dataset.
vennSets <- list("Plant criteria" = f1a$`plant criteria`,
                 "Database match" = f1a$`database match`,
                 "ZE"             = f1a$ZE,
                 "SE"             = f1a$SE)
pA <- ggVennDiagram(vennSets, label = "count", label_alpha = 0,
                    label_size = SZ_TICK / .pt, set_size = SZ_ANN / .pt,
                    edge_size = 0.3) +
  ## region counts only, with no fill gradient
  scale_fill_gradient(low = "white", high = "white", guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(title = NULL) +
  theme_void_ms()

## ---- B  UpSet of the 961 retained miRNAs across the five tissue types ---------
tis <- c("SD", "FG", "Zem", "Sem", "SDL")
memb <- sapply(tis, function(t) focal %in% f1b[[t]])
rownames(memb) <- focal
sig <- apply(memb, 1, function(x) paste(tis[x], collapse = "|"))
isz <- sort(table(sig), decreasing = TRUE)
isz <- isz[seq_len(min(12, length(isz)))]
setlab <- setNames(sprintf("%s (%d)", tis, lengths(f1b)[tis]), tis)
bars <- tibble(sig = names(isz), n = as.integer(isz)) %>%
  mutate(sig = factor(sig, levels = names(isz)))
## the same bars split by known / variant / novel
barsC <- tibble(miR = focal, sig = sig, class = miR_class[focal]) %>%
  filter(sig %in% levels(bars$sig)) %>%
  count(sig, class) %>%
  mutate(sig = factor(sig, levels = levels(bars$sig)),
         class = factor(class, levels = rev(CLS_ORDER)))
stopifnot(sum(barsC$n) == sum(bars$n))
mtx <- bars %>% rowwise() %>%
  do(tibble(sig = .$sig, tissue = tis,
            member = tis %in% strsplit(as.character(.$sig), "\\|")[[1]])) %>%
  ungroup() %>%
  mutate(tissue = factor(setlab[tissue], levels = rev(setlab)))
## single ggplot: intersection bars above zero, membership matrix below it
step <- max(bars$n) * 0.105
mtx <- mtx %>% mutate(row = match(as.character(tissue), rev(setlab)),
                      y = -step * row)
brk <- c(pretty(c(0, max(bars$n)), 4), -step * seq_along(setlab))
lbl <- c(format(pretty(c(0, max(bars$n)), 4), big.mark = ","), rev(setlab))
cat("retained miRNAs:", length(focal), "| in all five tissues:", sum(rowSums(memb) == length(tis)), "\n")
cat("tissue sizes:", paste(names(f1b), lengths(f1b), collapse = ", "), "\n")
print(isz)
pB <- ggplot() +
  geom_col(data = barsC, aes(sig, n, fill = class), width = 0.72) +
  scale_fill_manual(values = COL_CLASS, breaks = CLS_ORDER,
                    name = "miRNA class") +
  geom_text(data = bars, aes(sig, n, label = n), vjust = -0.4,
            size = SZ_TICK / .pt, colour = "grey20") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_line(data = subset(mtx, member), aes(sig, y, group = sig),
            colour = "grey25", linewidth = 0.4) +
  geom_point(data = mtx, aes(sig, y, colour = member), size = 1.3) +
  scale_colour_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey87"),
                      guide = "none") +
  scale_y_continuous(breaks = brk, labels = lbl,
                     expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = NULL, y = "miRNAs in intersection",
       title = NULL) +
  theme_ms() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey93"),
        ## the title starts at the plot edge because it is wider than the panel
        plot.title.position = "plot",
        plot.title = element_text(size = SZ_BASE, hjust = 0, margin = margin(t = 2, b = 2)))

## ---- C  stage-mean expression heatmap -------------------------------------
stage_lab <- grp %>%
  mutate(lab = case_when(tissue == "SD"  ~ paste0("SD", stage),
                         tissue == "FG"  ~ paste0("FG", stage),
                         tissue == "Zem" ~ paste0("Zem", stage),
                         TRUE            ~ paste0("S", stage))) %>%
  mutate(lab = factor(lab, levels = intersect(STAGE_ORDER, unique(lab))))
stopifnot(!any(is.na(stage_lab$lab)), !any(is.na(stage_lab$tissue)))
stagemean <- sapply(split(stage_lab$lib, stage_lab$lab),
                    function(ls) rowMeans(vall[, ls, drop = FALSE]))
## scaled within each miRNA, so the panel shows change across stages rather
## than abundance: the UpSet and the DE panels already cover abundance
## a few low-expressed miRNAs have the same mean in every stage; scaling leaves
## them undefined, so they are shown as no change rather than dropped
stagez <- t(scale(t(stagemean)))
flat <- apply(stagemean, 1, function(x) sd(x) == 0)
stagez[flat, ] <- 0
cat("miRNAs with no variation across stages:", sum(flat), "\n")
stopifnot(all(is.finite(stagez)))
ord <- hclust(dist(stagez))$order
hm <- as.data.frame(stagez) %>% mutate(miR = rownames(stagez)) %>%
  pivot_longer(-miR, names_to = "stage", values_to = "vst") %>%
  left_join(distinct(stage_lab, stage = lab, tissue) %>%
              mutate(stage = as.character(stage)), by = "stage") %>%
  mutate(miR = factor(miR, levels = rownames(stagez)[ord]),
         stage = factor(stage, levels = levels(stage_lab$lab)))
pC <- ggplot(hm, aes(stage, miR, fill = vst)) +
  geom_tile() +
  facet_grid(cols = vars(factor(tissue, levels = tis)), scales = "free_x",
             space = "free_x") +
  scale_fill_gradient2(low = COL_DIR[["down"]], mid = "white",
                       high = COL_DIR[["up"]], midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "z-score",
                       guide = guide_colourbar(barwidth = unit(2.5, "mm"),
                                               barheight = unit(11, "mm"),
                                               title.position = "top")) +
  scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = sprintf("%d retained miRNAs", length(focal)),
       title = NULL) +
  theme_ms() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(), legend.position = "right",
        legend.key.width = unit(2.5, "mm"), legend.key.height = unit(8, "mm"))

## ---- D-F  principal component analyses ------------------------------------
pca_panel <- function(mat, libs, title) {
  m <- mat[, libs, drop = FALSE]
  m <- m[apply(m, 1, var) > 0, , drop = FALSE]
  pc <- prcomp(t(m))
  pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  cat(title, ": ", nrow(m), " miRNAs, PC1 ", pv[1], "%, PC2 ", pv[2], "%\n", sep = "")
  d <- as.data.frame(pc$x[, 1:2]) %>%
    mutate(lib = rownames(pc$x)) %>%
    left_join(grp, by = "lib") %>%
    mutate(tissue = factor(tissue, levels = tis))
  ## one label per tissue and stage, placed at the group centroid
  cen <- d %>% group_by(tissue, stage) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  ggplot(d, aes(PC1, PC2, colour = tissue)) +
    geom_point(size = 1.1, alpha = 0.9) +
    ggrepel::geom_text_repel(data = cen, aes(label = stage), size = SZ_TICK / .pt,
                             show.legend = FALSE, min.segment.length = 0.15,
                             segment.size = 0.15, box.padding = 0.15,
                             point.padding = 0.1, max.overlaps = 40,
                             bg.color = "white", bg.r = 0.12, seed = 1) +
    ## seedlings are excluded from every expression comparison, so the scale is
    ## limited to the tissues present here; SDL still appears in panel B, which
    ## reports detection rather than expression
    scale_colour_manual(values = COL_TISSUE,
                        limits = intersect(tis, unique(as.character(grp$tissue))),
                        drop = FALSE, name = "tissue") +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
         title = title) +
    theme_ms() + theme(axis.text = element_blank(), axis.ticks = element_blank())
}
pD <- pca_panel(vall, grp$lib, "Both series")
pE <- pca_panel(vall, grp$lib[grp$series == "ZE"], "Zygotic series") +
  theme(legend.position = "none")
pF <- pca_panel(vall, grp$lib[grp$series == "SE"], "Somatic series") +
  theme(legend.position = "none")

design <- "
AABB
CCCC
DEFG
"
fig1 <- pA + pB + pC + pD + pE + pF + guide_area() +
  plot_layout(design = design, heights = c(1.05, 1.4, 0.95),
              widths = c(1, 1, 1, 0.5), guides = "collect") +
  plot_annotation(tag_levels = "A") & tag_theme()
save_fig(fig1, "Figure1", MM_DOUBLE, 210)
