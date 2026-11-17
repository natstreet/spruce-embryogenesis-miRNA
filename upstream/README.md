# Processing from raw reads to the deposited inputs

The scripts in this folder were run on a high-performance computing cluster (SLURM) to process the raw
sequencing reads into the inputs deposited in the SciLifeLab Data Repository item (`data/`). They are
provided as run, to document the steps described in Methods S1 of the article. They cannot be run
without the raw reads (ENA PRJEB72619), the *Picea abies* 2.0 genome and transcriptome (ENA PRJEB69221)
and the external tools listed below.

Theerarat Kochakarn wrote these scripts. `pipeline/preprocessing_for_mirDP2_fastqinput.py` is modified
from a script by Ahsan Z. Rizvi.

## Paths
Site-specific paths are replaced by variables. Set them before running a shell script:

- `PROJECT_DIR`: the project folder, with `pipeline/`, `doc/` and `data/analysis/` below it
- `REFERENCE_DIR`: genome, transcriptome, annotation and SortMeRNA reference files
- `TOOLS_DIR`: local installations of miRDeep-P2, the miRDeep-P2 pipeline scripts, TargetFinder and FASTA
- `CONTAINER_DIR`: Apptainer/Singularity images (SortMeRNA, Salmon, psRobot, TAPIR)
- `BIND_DIR`: the file system mounted into the containers

`#SBATCH -A <account>` stands for the cluster account. The R Markdown files use paths relative to the
project folder (`here()`), with `reference/` for reference files.

## Steps, in the order of Methods S1

| Step | Scripts |
|---|---|
| Read quality (FastQC) | `pipeline/runFastQC.sh`, `pipeline/submitFastQC.sh` |
| Adapter, quality and length trimming (fastp) | `pipeline/runFastp.sh`, `pipeline/submitFastp.sh`, `pipeline/nfcore_known_adapters.fa` |
| Removal of rRNA and other non-coding RNA reads (SortMeRNA) | `pipeline/runSortMeRNA.sh`, `pipeline/submitSortMeRNA.sh` |
| Merging the sequencing runs of each library | `pipeline/runMergeTechrep.sh` |
| Collapsing reads to unique sequences | `pipeline/step0.sh`, `pipeline/preprocessing_for_mirDP2_fastqinput.py` |
| Bowtie index of the genome | `pipeline/bowtieindex.sh` |
| miRNA prediction (miRDeep-P2) | `pipeline/step1.sh` |
| Combining predictions and annotation by BLAST against miRBase, PMRD and PmiREN | `pipeline/step2_moreMiR.sh`, `R/miRNAannotation.R` |
| Counting reads per predicted miRNA (Bowtie) | `pipeline/step3.sh` |
| miRNA DESeq2 objects and plant criteria | `R/two-tier_flipComp_DifferentialExpression-SE_final_omitSeedling_moreDB.Rmd`, `R/two-tier_flipComp_DifferentialExpression-ZE_finalMetadata_moreDB.Rmd` |
| Replicate check (the replicate-checked sequences) | `R/two-tier_post_miRDP2_BioQA_repCheck.Rmd` |
| Target prediction: TargetFinder, psRobot, TAPIR (FASTA and RNAhybrid engines) | `pipeline/runTargetFinder_all.sh`, `pipeline/runPsrobot.sh`, `pipeline/submitPsrobot_all.sh`, `pipeline/runTapirFasta_all.sh`, `pipeline/runTapirHybrid.sh`, `pipeline/submitTapirHybrid_all.sh`, `pipeline/submitTapirHybrid_all_fix.sh` |
| Pairs predicted by at least three of four programs, and miRNA vst matrices | `R/two-tier_flipComp_postDE-SEZE_withCorrPadj_repCheck_omitSeedling.Rmd` |
| mRNA quantification (Salmon) | `pipeline/runSalmon.sh`, `pipeline/submitSalmon.sh` |
| mRNA DESeq2 objects and vst matrices | `R/flipComp_mRNADifferentialExpression-ZE_finalMetadata.Rmd`, `R/flipComp_mRNADifferentialExpression-SE_omitSeedling.Rmd` |
| miRNA catalogue with family annotation, and gene annotation | `R/makeMasterTable.R`, `R/cleanup_miRannotations.Rmd` |

psRNATarget was run on its web server (https://www.zhaolab.org/psRNATarget/) with default settings;
its output is read by `R/two-tier_flipComp_postDE-SEZE_withCorrPadj_repCheck_omitSeedling.Rmd`.
`submitTapirHybrid_all_fix.sh` repeats the TAPIR RNAhybrid runs for the miRNA batches that did not
finish in the first submission.

Some of these scripts also test differential expression, correlate miRNAs with their targets or run GO
enrichment. The results reported in the article come from the scripts in `scripts/` (see
`REPRODUCE.md`), which start from the deposited inputs.

## Software
FastQC 0.12.1, fastp 0.23.4, SortMeRNA 4.3.7, Bowtie 1.3.1, miRDeep-P2 1.1.4 with the miRDeep-P2
pipeline scripts (https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline), BLAST+ 2.16.0, Salmon 1.10.3,
TargetFinder (https://github.com/carringtonlab/TargetFinder), psRobot 1.2 and TAPIR, with the
versions loaded as modules or containers in each script.
