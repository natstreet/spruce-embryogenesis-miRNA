## Shared figure style, sourced by every figure script.
## Physiologia Plantarum: 84 mm single column, 174 mm double column.

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(scales)
})

invisible(Sys.setlocale("LC_COLLATE", "C"))   # sort order of labels and legend keys independent of the system locale

MM_SINGLE <- 84
MM_DOUBLE <- 174

## three font sizes, one for each role
SZ_BASE <- 8   # titles, axis labels, series identity
SZ_ANN  <- 7   # legend and annotation
SZ_TICK <- 6   # tick labels

## Okabe-Ito colours, one per tissue
COL_TISSUE <- c(SD = "#CC79A7", FG = "#D55E00", Zem = "#0072B2",
                Sem = "#009E73", SDL = "#E69F00")
## series colours: reddish purple for the zygotic series and green for the somatic
## series; vermillion means "higher in the later stage" in Figures 2, 5 and 6
COL_SERIES <- c(ZE = "#CC79A7", SE = "#009E73")
COL_DIR    <- c(up = "#D55E00", down = "#0072B2")
COL_ALARM  <- "#D55E00"   # reserved for flagged marks only
COL_GREY   <- "grey55"
COL_DIVERGE <- c("#2166AC", "#F7F7F7", "#B2182B")

theme_ms <- function(base = SZ_BASE) {
  theme_bw(base_size = base, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
      panel.border = element_rect(linewidth = 0.3, colour = "grey30"),
      axis.text = element_text(size = SZ_TICK, colour = "grey20"),
      axis.title = element_text(size = base),
      axis.ticks = element_line(linewidth = 0.25, colour = "grey30"),
      plot.title = element_text(size = base, hjust = 0, face = "plain",
                                margin = margin(b = 2)),
      legend.title = element_text(size = SZ_ANN),
      legend.text = element_text(size = SZ_ANN),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      strip.background = element_rect(fill = "grey95", linewidth = 0.3,
                                      colour = "grey30"),
      strip.text = element_text(size = SZ_ANN, margin = margin(1.5, 1.5, 1.5, 1.5)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

theme_void_ms <- function(base = SZ_BASE) {
  theme_void(base_size = base, base_family = "Helvetica") +
    theme(legend.title = element_text(size = SZ_ANN),
          legend.text = element_text(size = SZ_ANN),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = base, hjust = 0,
                                    margin = margin(b = 2)),
          plot.margin = margin(2, 2, 2, 2))
}

## panel letters: uppercase, bold, in the panel margin so they never sit on a title
tag_theme <- function() {
  theme(plot.tag = element_text(size = SZ_BASE + 2, face = "bold",
                                family = "Helvetica"),
        plot.tag.position = "topleft",
        plot.tag.location = "margin",
        plot.margin = margin(2, 2, 2, 5))
}

save_fig <- function(p, stem, width_mm, height_mm) {
  dir.create("figures", showWarnings = FALSE)
  ggsave(paste0("figures/", stem, ".pdf"), plot = p, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0("figures/", stem, ".png"), plot = p, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
  cat("wrote", stem, sprintf("(%.0f x %.0f mm)\n", width_mm, height_mm))
}
