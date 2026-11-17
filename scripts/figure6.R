## Figure 6: expression support for predicted targets, and the miR159-GAMYB module.
## Inputs:
##   results/target_support_chance_shuffled_targets.tsv  observed vs shuffled-target pairs
##   results/target_support_chance_random_genes.tsv     observed vs miRNA paired with random genes
##   results/miR159_GAMYB_vst_per_library.tsv           vst per library, centred per feature within series
source("scripts/fig_style.R")

COL_FOCAL <- "#D55E00"   # observed value and miR159

KEY_O <- "observed"; KEY_E <- "expected (median, 95% range)"
key_scale <- scale_colour_manual(values = setNames(c(COL_FOCAL, COL_GREY), c(KEY_O, KEY_E)), breaks = c(KEY_O, KEY_E), name = NULL,
                                 guide = guide_legend(override.aes = list(shape = c(16, 124), size = c(1.8, 3))))

## ---- A: shares for the large pair sets, against shuffled targets -------------
## opp_frac is pairs DE in opposite directions / (pairs DE in opposite directions +
## pairs DE in the same direction); a pair that is opposite in one comparison and
## same in another (37 of all predicted pairs) is counted in both
sh <- read_tsv("results/target_support_chance_shuffled_targets.tsv", show_col_types = FALSE)
setsA <- c("all predicted" = "All predicted pairs\n(26,812)",
           "four programs" = "Four-program pairs\n(18,482)",
           "conserved transcription-factor modules" = "Conserved TF modules\n(159)")
A <- sh %>% filter(subset %in% names(setsA), stat %in% c("opp_frac", "neg_frac")) %>%
  mutate(set = factor(setsA[subset], levels = rev(setsA)),
         stat = factor(c(opp_frac = "Opposite DE\n(% of opposite + same)",
                         neg_frac = "Negative correlation\n(% of significant)")[stat],
                       levels = c("Opposite DE\n(% of opposite + same)", "Negative correlation\n(% of significant)")))
pA <- ggplot(A, aes(y = set)) +
  geom_vline(xintercept = 50, linewidth = 0.25, colour = "grey70", linetype = "dashed") +
  geom_errorbarh(aes(xmin = 100 * null_lo, xmax = 100 * null_hi), height = 0.25,
                 linewidth = 0.4, colour = COL_GREY, show.legend = FALSE) +
  geom_point(aes(x = 100 * null_median, colour = KEY_E), shape = 124, size = 3) +
  geom_point(aes(x = 100 * observed, colour = KEY_O), size = 1.8) + key_scale +
  facet_wrap(~ stat, nrow = 1) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 10)) +
  labs(x = "% of pairs", y = NULL,
       title = NULL) +
  theme_ms() + tag_theme()

## ---- B: miR159-MYB counts against random genes ---------------------------------
## 28 miR159-MYB pairs, 27 of them with the three GAMYB genes; the 28th
## (Picea_mature_4734 with PA_chr06_G004097) is a MYB gene without the GAMYB label
rg <- read_tsv("results/target_support_chance_random_genes.tsv", show_col_types = FALSE)
labB <- c(opp = "DE in opposite directions", same = "DE in the same direction",
          negcor = "Negatively correlated", poscor = "Positively correlated")
B <- rg %>% filter(subset == "miR159, MYB targets", stat %in% names(labB)) %>%
  mutate(stat = factor(labB[stat], levels = rev(labB)))
pB <- ggplot(B, aes(y = stat)) +
  geom_errorbarh(aes(xmin = null_lo, xmax = null_hi), height = 0.25, linewidth = 0.4, colour = COL_GREY, show.legend = FALSE) +
  geom_point(aes(x = null_median, colour = KEY_E), shape = 124, size = 3) +
  geom_point(aes(x = observed, colour = KEY_O), size = 1.8) + key_scale +
  scale_x_continuous(limits = c(0, 13), breaks = seq(0, 12, 4)) +
  labs(x = expression("miR159 and "*italic(MYB)*" pairs (of 28)"), y = NULL,
       title = NULL) +
  theme_ms() + tag_theme()

## ---- C: miR159 and GAMYB across stages ---------------------------------------
cdat <- read_tsv("results/miR159_GAMYB_vst_per_library.tsv", show_col_types = FALSE)
ordZ <- c("SD1", "SD2", "SD3", "Zem4", "Zem5", "Zem6", "Zem7", "Zem8-9", "Zem10")
ordS <- paste0("S", 1:6)
## each GAMYB gene drawn separately, in three shades of blue with three line types
GAM <- c("PA_chr07_G000337", "PA_chr08_G002733", "PA_chr03_G008015")
FEAT <- c("miR159", GAM)
COL_C <- setNames(c(COL_FOCAL, "#08306B", "#0072B2", "#56B4E9"), FEAT)
LTY_C <- setNames(c("solid", "solid", "dashed", "dotted"), FEAT)
LAB_C <- as.expression(c(quote("miR159 (mean of 6 embryo-enriched members)"),
                         lapply(GAM, function(g) bquote(italic(GAMYB)*" "*.(g)))))
stopifnot(setequal(unique(cdat$feature), FEAT))
cdat <- cdat %>% mutate(series = factor(series, levels = c("Zygotic series", "Somatic series")),
                        stage = factor(stage, levels = c(ordZ, ordS)),
                        feature = factor(feature, levels = FEAT))
means <- cdat %>% group_by(series, stage, feature) %>% summarise(value = mean(value), .groups = "drop")
## shaded transitions, in positions within each facet's discrete axis
bands <- tibble(series = factor(c("Zygotic series", "Somatic series", "Somatic series"), levels = levels(cdat$series)),
                xmin = c(3, 2, 5), xmax = c(4, 3, 6),
                label = c("seed to\nembryo", "to ABA\nmedium", "desiccation"))
yr <- range(cdat$value)
pC <- ggplot(means, aes(stage, value, colour = feature, group = feature)) +
  geom_blank() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey92") +
  geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = yr[2] + 0.9, label = label), inherit.aes = FALSE,
            size = SZ_TICK / .pt, colour = "grey35", lineheight = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(data = cdat, size = 0.5, alpha = 0.45) +
  geom_line(aes(linetype = feature), linewidth = 0.55) +
  facet_grid(~ series, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = COL_C, labels = LAB_C, name = NULL) +
  scale_linetype_manual(values = LTY_C, labels = LAB_C, name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.6, alpha = 1)),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(limits = c(yr[1] - 0.2, yr[2] + 1.6)) +
  labs(x = NULL, y = "vst relative to series mean",
       title = NULL) +
  theme_ms() + tag_theme() +
  theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0))

pA <- pA + theme(legend.position = "bottom"); pB <- pB + theme(legend.position = "none")
fig <- (pA | pB) / free(pC) + plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")
save_fig(fig, "Figure6", MM_DOUBLE, 130)
saveRDS(list(A = A, B = B, C = cdat), "results/figure6_panel_data.rds")
