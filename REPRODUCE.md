# Reproducing the article

Run the commands from the repository folder after unpacking `data.tar.gz` (see `README.md`).

## 1. Analyses
Run in this order. Each script writes to `results/`.

| Script | Output |
|---|---|
| `Rscript scripts/DE_stages_voomQWB.R` | stage-to-stage DE for miRNAs and mRNAs (limma voomWithQualityWeights) |
| `Rscript scripts/correlation_support.R` | within-dataset correlation of DE miRNAs with their predicted target genes |
| `Rscript scripts/stage_mapping_mRNA.R` | mRNA stage mapping between the somatic and zygotic series (Figure 5, lower panel) |
| `Rscript scripts/S5_mapping_all_Zem_stages.R` | correlation of S5 with every zygotic embryo stage for miRNA and target gene sets (Figure S5B) |
| `Rscript scripts/TableS7_overlap.R` | overlap of DE miRNAs between somatic culture transitions and zygotic transitions (Table S7) |
| `Rscript scripts/TableS11_consistent_direction.R` | change at the matching developmental phase of the other series for miRNAs DE in only one series (Table S11) |
| `Rscript scripts/stage_mapping_contributions.R` | contributions of single miRNAs and enrichment classes to the miRNA stage mapping |
| `Rscript scripts/S6_FG89_GO.R` | GO enrichment of the genes behind the mRNA match of S6 with FG8-9 |
| `Rscript scripts/GO_enrichment.R` | GO enrichment of target genes of supported pairs by series and by miRNA family (Tables S4, S5, S10; Figure 4) |

## 2. Figures
Each script writes a PDF and a PNG to `figures/`.

| Script | Figure |
|---|---|
| `Rscript scripts/figure1.R` | Figure 1 |
| `Rscript scripts/figure2.R` | Figure 2 |
| `Rscript scripts/figure3.R` | Figure 3 |
| `Rscript scripts/figure4.R` | Figure 4 |
| `Rscript scripts/figure5.R` | Figure 5 |
| `Rscript scripts/figure6.R` | Figure 6 |
| `Rscript scripts/figureS1_S6.R` | Supplementary Figures S1 and S6 |
| `Rscript scripts/figureS2.R` | Supplementary Figure S2 |
| `Rscript scripts/figureS3.R` | Supplementary Figure S3 |
| `Rscript scripts/figureS4.R` | Supplementary Figure S4 |
| `Rscript scripts/figureS5.R` | Supplementary Figure S5 |
| `Rscript scripts/figureS7.R` | Supplementary Figure S7 |

`scripts/fig_style.R` holds the shared figure style and is sourced by every figure script.

## 3. Tables read but not written by these scripts
Some tables in `results/` are processed tables behind the article that no script in this repository
writes. Among them are the paired embryo-versus-gametophyte test (Supplementary Table S6), the miRNA
classes and families, the target-pair support flags and annotation (Supplementary Table S8), the chance
comparisons for target support (Supplementary Table S9), the embryo signature index, the miRNA stage
mapping and the Figure 1 panel data. `FILE_DESCRIPTIONS.txt` in the data archive describes each of them.

## 4. Software versions
R 4.4.3 with DESeq2 1.46.0, limma 3.62.1, edgeR 4.4.0, topGO 2.58.0 and GO.db 3.20.0
(Bioconductor 3.20); `environment.yml` gives all versions. Other Bioconductor releases can change the
differential expression and GO results slightly.
