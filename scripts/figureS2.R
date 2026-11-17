## Supplementary Figure S2 | expression of the 961 retained miRNAs in every
## library, scaled within each miRNA. Seedling libraries are excluded.
## Library tissue and stage are taken from results/figure1_panel_data.rds (dat$grp, 68
## libraries), which lists the same libraries as the vst objects after the three
## low-quality libraries were removed.
source("scripts/fig_style.R")
suppressPackageStartupMessages(library(SummarizedExperiment))

dat <- readRDS("results/figure1_panel_data.rds")
foc <- dat$focal
vze <- readRDS("data/miRNA_vst_ZE.rds")
vse <- readRDS("data/miRNA_vst_SE.rds")
gm <- function(v) as.matrix(if (inherits(v, "DESeqTransform")) assay(v) else v)
mZ <- gm(vze); mS <- gm(vse)
lowQual <- c("PabK14-03SM_SR1", "PabK14-03RO_SR2", "Pab738_2_FMG_MIR1")
mZ <- mZ[, !colnames(mZ) %in% lowQual, drop = FALSE]
mS <- mS[, !colnames(mS) %in% lowQual, drop = FALSE]
V <- cbind(mZ[foc, , drop = FALSE], mS[foc, , drop = FALSE])
## miRNAs with no reads in any of these libraries were retained on the seedling
## libraries only; with no variance they cannot be scaled, and are left out
nz <- apply(V, 1, sd) > 0
stopifnot(all(rownames(V)[!nz] %in% dat$f1b$SDL))
cat("S2: left out", sum(!nz), "miRNAs detected only in seedlings\n")
V <- V[nz, , drop = FALSE]

ann <- dat$grp; rownames(ann) <- ann$lib
stopifnot(setequal(colnames(V), ann$lib))
ann <- ann[colnames(V), ]
tiss <- ann$tissue
tiss[tiss == "ZE"] <- "Zem"
TORD <- c("SD", "FG", "Zem", "Sem", "SDL")
snum <- suppressWarnings(as.numeric(gsub("8-9", "8.5", ann$stage)))
ord <- order(match(tiss, TORD), snum)
Vo <- V[, ord]; tis_o <- tiss[ord]; stg_o <- ann$stage[ord]
Z <- t(scale(t(Vo))); Z[!is.finite(Z)] <- 0
rord <- hclust(dist(Z))$order
long <- data.frame(miR = factor(rep(rownames(Z)[rord], times = ncol(Z)), levels = rownames(Z)[rord]),
                   lib = factor(rep(colnames(Z), each = nrow(Z)), levels = colnames(Z)),
                   z = as.vector(Z[rord, ]),
                   tis = factor(rep(tis_o, each = nrow(Z)), levels = TORD), stringsAsFactors = FALSE)
labmap <- setNames(paste0(ifelse(tis_o == "Sem", "S", tis_o), stg_o), colnames(Z))   # somatic stages are S1-S6 throughout
SZ <- 8
gS2 <- ggplot(long, aes(lib, miR, fill = z)) + geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-2, 2),
                       oob = scales::squish, name = "z-score") +
  facet_grid(. ~ tis, scales = "free_x", space = "free_x") + scale_x_discrete(labels = labmap) +
  theme_bw(base_size = SZ) + theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = .5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
    strip.background = element_blank(), legend.position = "bottom", legend.key.height = unit(2, "mm"),
    legend.key.width = unit(8, "mm")) + labs(x = NULL, y = sprintf("%d retained miRNAs", nrow(Z)))
cat("S2 libraries:", ncol(Z), "| miRNAs:", nrow(Z), "| tissues:", paste(TORD[TORD %in% tis_o], collapse = ","), "\n")
save_fig(gS2, "FigureS2", MM_DOUBLE, 150)
