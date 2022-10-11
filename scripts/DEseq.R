#install.packages("DESeq2")
#install.packages("here")
#install.packages("pheatmap")
#install.packages("apeglm")
#install.packages("ashr")
#remotes::install_github("stephenslab/mixsqp")
#install.packages("plotly")

library(DESeq2)
library(tidyverse)
library(here)
library(pheatmap)
library(apeglm)
library(knitr)
library(ashr)
library(mixsqp)
library(plotly)

theme_set(theme_bw())
getwd()
outdir <- here("results/RNA-Seq/DESeq2/")
plotdir <- here("results/RNA-Seq/DESeq2/plots")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

## Tutorial workflow if needed from http://master.bioconductor.org/packages/release/workflows/html/rnaseqGene.html
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rnaseqGene")
#########PREPARING THE DATA######################
###COMPILED FROM: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html,
###     https://biodash.github.io/tutorials/2021-01_rnaseq/09-DE.html#Getting_set_up
## Import counts files
psan_counts_file <- "/fs/project/PAS0471/linda/results/RNA-Seq/nf_dual/psan/HTSeq/pathogen_quantification_uniquely_mapped_htseq.tsv"
psan_count <- read_tsv(psan_counts_file)

## Clean sample names
psan_counts_mat <- psan_count %>% select(contains("_NumReads")) %>% as.matrix()
rownames(psan_counts_mat) <- psan_count$gene_id
colnames(psan_counts_mat) <- sub("^X", "", colnames(psan_counts_mat))
colnames(psan_counts_mat) <- sub("_S\\d+_001.*", "", colnames(psan_counts_mat))
head(psan_counts_mat)
psan_host_counts_mat <- psan_host_count %>% select(contains("_NumReads")) %>% as.matrix()

### DATASET WITH OUTLIERS REMOVED
## Remove psan samples: 24S10R3, 24S68R6, 48S68R4, 72S68R3
## Remove psoj samples: 24SRR5, 24SRR6
## Compile metadata.txt files for 80 samples
psan_60 <- psan_counts_mat[, c(1:4,6:19,29:43,45:46,55:72,74:76)]

### METADATA FILES NOW HAVE OUTLIERS REMOVED
Psan_metadata <- read.csv("/fs/project/PAS0471/linda/results/RNA-Seq/DESeq2/Psan_metadata.csv",
                          header = TRUE, na.strings = "NA")
Psan_metadata$Timepoint = as.factor(Psan_metadata$Timepoint)
Psan_metadata$Rep = as.factor(Psan_metadata$Rep)

Psan_metadata_ordered <- Psan_metadata[order(Psan_metadata$Sample), ]
Psan_metadata_ordered$Sample
psan_counts_ordered <- psan_60[, order(colnames(psan_60))]
colnames(psan_counts_ordered)
matching_names <- identical(Psan_metadata_ordered$Sample, colnames(psan_counts_ordered))
matching_names

## Creating DESeq2 objects
dseq_psan_raw <- DESeqDataSetFromMatrix(countData = psan_60,
                                        colData = Psan_metadata,
                                        design = ~ 1)

###############DE ANALYSIS - FULL DATASET##################
## Our current design has two factors: 'Timepoint' and 'Treatment'
#first start by merging into a single factor called 'group'
dseq_psan_raw$group <- factor(paste(dseq_psan_raw$Timepoint, dseq_psan_raw$Treatment,
                                    sep = "_"))
table(dseq_psan_raw$group)

## Reference = mock for each isolate
      #use Psan68 as 'true' mock for psan because all 3 isolates aligned to Psan68
dseq_psan_raw$group <- relevel(dseq_psan_raw$group, ref = "0_Psan68")
dseq_psan_raw$group

## Set analysis design
design(dseq_psan_raw) <- formula(~ group)

## Finally perform differential expression analysis
psan_dseq <- DESeq(dseq_psan_raw)

psan_res <- results(psan_dseq)
resultsNames(psan_dseq)

## Function to get nr of sig results
# psan_sig0.1_contrast <- function(contrast_psan, psan_dseq) {
#   res_sig <- results(psan_dseq,
#                      contrast = c("group", contrast_psan)) %>%
#     as.data.frame() %>%
#     dplyr::filter(padj < 0.1) %>%
#     mutate(level1 = contrast_psan[1],
#            level2 = contrast_psan[2])
#   cat(contrast_psan[1], "versus", contrast_psan[2], ":", nrow(res_sig), "significant\n")
#   return(res_sig)
# }
# psan_comps <- combn(levels(psan_dseq@colData$group), 2, simplify = FALSE)
# psan_sig0.1_all_contrasts <- do.call(rbind, lapply(psan_comps, psan_sig0.1_contrast, psan_dseq))

# 0_Psan68_vs_0_Psan10 -- everything the same
res <- results(psan_dseq, contrast = c("group", "0_Psan10", "0_Psan68"))
sum(res$padj < 0.05, na.rm = TRUE)

# 24_Psan10_vs_24_Psan65 -- 0 vs 54
res <- results(psan_dseq, contrast = c("group", "24_Psan10", "24_Psan65"))
sum(res$padj < 0.05, na.rm = TRUE)

###########DE ANALYSIS - WITH TWO FACTORS (REF = PSAN68)###########################
## Controlling for one factor:
dseq_2f_psan_raw <- dseq_psan_raw
dseq_2f_psan_raw$Treatment <- relevel(factor(dseq_2f_psan_raw$Treatment), ref = "Psan68")
dseq_2f_psan_raw$Timepoint <- relevel(factor(dseq_2f_psan_raw$Timepoint), ref = "0")
#order matters: test for the effect of the last factor while controlling for the effect of the first factor
#effect of Treatment while controlling for Timepoint:
design(dseq_2f_psan_raw) <- formula(~ Timepoint + Treatment)
psan_dseq_2fTi <- DESeq(dseq_2f_psan_raw)

resultsNames(psan_dseq_2fTi)
res <- results(psan_dseq_2fTi)
sum(res$padj < 0.05, na.rm = TRUE) 

#design(dseq_2f_psan_raw) <- formula(~ Treatment + Timepoint)
#psan_dseq_2fTr <- DESeq(dseq_2f_psan_raw)
#resultsNames(psan_dseq_2fTr)

#save a new object 
dseq_2fi_irixr_psan_raw <- dseq_2f_psan_raw
design(dseq_2fi_irixr_psan_raw) <- formula(~Timepoint + Treatment + Timepoint:Treatment)
psan_dseq_2fi_irixr <- DESeq(dseq_2fi_irixr_psan_raw)
resultsNames(psan_dseq_2fi_irixr)

## Do effects of Treatment differ among levels of Timepoint?
#interaction is last term in formula, so called by 'results()' function
res <- results(psan_dseq_2fi_irixr)
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 19
## Effect of Treatment for "control" Timepoint (0):
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Treatment","Psan10","Psan65"))
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 74
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Treatment","Psan10","Psan68"))
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 9
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Treatment","Psan65","Psan68"))
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 255

## Effect of Treatment for Timepoint24:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint24.TreatmentPsan65")))
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 54

## Effect of Treatment for Timepoint48:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan10", "Timepoint48.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 30
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 11
## Effect of Treatment for Timepoint72:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint72.TreatmentPsan10", "Timepoint72.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 46
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 30
## Effect of Timepoint for "control" Treatment (Psan68):
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 196
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 104
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 302
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 189
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","48", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 2
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 2
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "24"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 766
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 567
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 958
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 741
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 1165
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 937
## Effect of Timepoint for Treatment Psan10:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint48.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint72.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan10", "Timepoint72.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Timepoint for Treatment Psan65:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan65", "Timepoint48.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 232
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 145
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan65", "Timepoint72.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 250
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 144
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan65", "Timepoint72.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 252
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 155

###########DE ANALYSIS - WITH TWO FACTORS (REF = PSAN65)###########################
## Controlling for one factor:
#save a new object
dseq_2f_psan_raw <- dseq_psan_raw
dseq_2f_psan_raw$Treatment <- relevel(factor(dseq_2f_psan_raw$Treatment), ref = "Psan65")
dseq_2f_psan_raw$Timepoint <- relevel(factor(dseq_2f_psan_raw$Timepoint), ref = "0")
#order matters: test for the effect of the last factor while controlling for the effect of the first factor
#effect of Treatment while controlling for Timepoint:
design(dseq_2f_psan_raw) <- formula(~ Timepoint + Treatment)
psan_dseq_2fTi <- DESeq(dseq_2f_psan_raw)
resultsNames(psan_dseq_2fTi)
res <- results(psan_dseq_2fTi)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 16
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 7
#effect of Timepoint while controlling for Treatment:
design(dseq_2f_psan_raw) <- formula(~ Treatment + Timepoint)
psan_dseq_2fTr <- DESeq(dseq_2f_psan_raw)
resultsNames(psan_dseq_2fTr)
res <- results(psan_dseq_2fTr)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 4,308
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 3,569
## With an interaction term:
#save a new object 
dseq_2fi_irixr_psan_raw <- dseq_2f_psan_raw
design(dseq_2fi_irixr_psan_raw) <- formula(~Timepoint + Treatment + Timepoint:Treatment)
psan_dseq_2fi_irixr <- DESeq(dseq_2fi_irixr_psan_raw)
resultsNames(psan_dseq_2fi_irixr)
## Do effects of Treatment differ among levels of Timepoint?
#interaction is last term in formula, so called by 'results()' function
res <- results(psan_dseq_2fi_irixr)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 35
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 19
## Effect of Treatment for "control" Timepoint (0):
#same
## Effect of Treatment for Timepoint24:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint24.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 74
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 43
## Effect of Treatment for Timepoint48:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan10", "Timepoint48.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 45
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 23
## Effect of Treatment for Timepoint72:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint72.TreatmentPsan10", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 191
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 101
## Effect of Timepoint for "control" Treatment (Psan65):
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 171
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 97
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 258
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 155
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","48", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 15
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 7
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "24"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 488
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 379
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 819
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 642
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 1203
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 902
## Effect of Timepoint for Treatment Psan10:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint48.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 29
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 15
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan10", "Timepoint72.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 18
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 11
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan10", "Timepoint72.TreatmentPsan10")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 7
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 6
## Effect of Timepoint for Treatment Psan68:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan68", "Timepoint48.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 232
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 145
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan68", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 250
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 144
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan68", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 252
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 155

###########DE ANALYSIS - WITH TWO FACTORS (REF = PSAN10)###########################
## Controlling for one factor:
#save a new object
dseq_2f_psan_raw <- dseq_psan_raw
dseq_2f_psan_raw$Treatment <- relevel(factor(dseq_2f_psan_raw$Treatment), ref = "Psan10")
dseq_2f_psan_raw$Timepoint <- relevel(factor(dseq_2f_psan_raw$Timepoint), ref = "0")
#order matters: test for the effect of the last factor while controlling for the effect of the first factor
#effect of Treatment while controlling for Timepoint:
design(dseq_2f_psan_raw) <- formula(~ Timepoint + Treatment)
psan_dseq_2fTi <- DESeq(dseq_2f_psan_raw)
resultsNames(psan_dseq_2fTi)
res <- results(psan_dseq_2fTi)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 6
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 6
#effect of Timepoint while controlling for Treatment:
design(dseq_2f_psan_raw) <- formula(~ Treatment + Timepoint)
psan_dseq_2fTr <- DESeq(dseq_2f_psan_raw)
resultsNames(psan_dseq_2fTr)
res <- results(psan_dseq_2fTr)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 4,308
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 3,569
## With an interaction term:
#save a new object 
dseq_2fi_irixr_psan_raw <- dseq_2f_psan_raw
design(dseq_2fi_irixr_psan_raw) <- formula(~Timepoint + Treatment + Timepoint:Treatment)
psan_dseq_2fi_irixr <- DESeq(dseq_2fi_irixr_psan_raw)
resultsNames(psan_dseq_2fi_irixr)
## Do effects of Treatment differ among levels of Timepoint?
#interaction is last term in formula, so called by 'results()' function
res <- results(psan_dseq_2fi_irixr)
sum(res$padj < 0.1, na.rm = TRUE) 
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Treatment for "control" Timepoint (0):
#same
## Effect of Treatment for Timepoint24:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan65", "Timepoint24.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Treatment for Timepoint48:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan65", "Timepoint48.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Treatment for Timepoint72:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint72.TreatmentPsan65", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Timepoint for "control" Treatment (Psan10):
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 65
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 45
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","24", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 284
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 166
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","48", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 2
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 2
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "24"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 344
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 264
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "48"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 699
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 550
res <- results(psan_dseq_2fi_irixr,
               contrast = c("Timepoint","0", "72"))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 561
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 444
## Effect of Timepoint for Treatment Psan68:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan68", "Timepoint48.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan68", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan68", "Timepoint72.TreatmentPsan68")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 0
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 0
## Effect of Timepoint for Treatment Psan65:
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan65", "Timepoint48.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 29
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 15
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint24.TreatmentPsan65", "Timepoint72.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 18
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 11
res <- results(psan_dseq_2fi_irixr,
               contrast = list(c("Timepoint48.TreatmentPsan65", "Timepoint72.TreatmentPsan65")))
sum(res$padj < 0.1, na.rm = TRUE)
# padj < 0.1 = 7
sum(res$padj < 0.05, na.rm = TRUE) 
# padj < 0.05 = 6

################DE ANALYSIS - CONTRAST TWO CUSTOM GROUPS##################
#cannot use apeglm because using contrast statement not includes in resultsNames()
resultsNames(psan_dseq)
#1
contrast_1 <- c("0_Psan10_mock", "0_Psan65_mock")
res_c1 <- results(psan_dseq, contrast = c("group", contrast_1))
sum(res_c1$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 44
res_c1 <- results(psan_dseq, contrast = c("group", contrast_1), lfcThreshold = 1)
res_LFC_c1 <- lfcShrink(psan_dseq,
                        res = res_c1,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c1)
head(res_LFC_c1)
sum(res_LFC_c1$svalue < 0.1, na.rm = TRUE)
  #svalue < 0.1 = 0

#2
contrast_2 <- c("0_Psan10_mock", "0_Psan68_mock")
res_c2 <- results(psan_dseq, contrast = c("group", contrast_2))
sum(res_c2$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 3
res_c2 <- results(psan_dseq, contrast = c("group", contrast_2), lfcThreshold = 1)

res_LFC_c2 <- lfcShrink(psan_dseq,
                        res = res_c2,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c2)
head(res_LFC_c2)
sum(res_LFC_c2$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 2

#3
contrast_3 <- c("0_Psan65_mock", "0_Psan68_mock")
res_c3 <- results(psan_dseq, contrast = c("group", contrast_3))
sum(res_c3$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 208
res_c3 <- results(psan_dseq, contrast = c("group", contrast_3), lfcThreshold = 1)
res_LFC_c3 <- lfcShrink(psan_dseq,
                        res = res_c3,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c3)
head(res_LFC_c3)
sum(res_LFC_c3$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 228

#4
contrast_4 <- c("0_Psan10_mock", "24_Psan10")
res_c4 <- results(psan_dseq, contrast = c("group", contrast_4))
sum(res_c4$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 541
res_c4 <- results(psan_dseq, contrast = c("group", contrast_4), lfcThreshold = 1)
res_LFC_c4 <- lfcShrink(psan_dseq,
                        res = res_c4,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c4)
head(res_LFC_c4)
sum(res_LFC_c4$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 2,217

#5
contrast_5 <- c("0_Psan10_mock", "48_Psan10")
res_c5 <- results(psan_dseq, contrast = c("group", contrast_5))
sum(res_c5$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 489
res_c5 <- results(psan_dseq, contrast = c("group", contrast_5), lfcThreshold = 1)
res_LFC_c5 <- lfcShrink(psan_dseq,
                        res = res_c5,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c5)
head(res_LFC_c5)
sum(res_LFC_c5$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,224

#6
contrast_6 <- c("0_Psan10_mock", "72_Psan10")
res_c6 <- results(psan_dseq, contrast = c("group", contrast_6))
sum(res_c6$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 376
res_c6 <- results(psan_dseq, contrast = c("group", contrast_6), lfcThreshold = 1)
res_LFC_c6 <- lfcShrink(psan_dseq,
                        res = res_c6,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c6)
head(res_LFC_c6)
sum(res_LFC_c6$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 953

#7
contrast_7 <- c("0_Psan65_mock", "24_Psan65")
res_c7 <- results(psan_dseq, contrast = c("group", contrast_7))
sum(res_c7$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 345
res_c7 <- results(psan_dseq, contrast = c("group", contrast_7), lfcThreshold = 1)
res_LFC_c7 <- lfcShrink(psan_dseq,
                        res = res_c7,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c7)
head(res_LFC_c7)
sum(res_LFC_c7$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,293

#8
contrast_8 <- c("0_Psan65_mock", "48_Psan65")
res_c8 <- results(psan_dseq, contrast = c("group", contrast_8))
sum(res_c8$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 593
res_c8 <- results(psan_dseq, contrast = c("group", contrast_8), lfcThreshold = 1)
res_LFC_c8 <- lfcShrink(psan_dseq,
                        res = res_c8,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c8)
head(res_LFC_c8)
sum(res_LFC_c8$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,642

#9
contrast_9 <- c("0_Psan65_mock", "72_Psan65")
res_c9 <- results(psan_dseq, contrast = c("group", contrast_9))
sum(res_c9$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 770
res_c9 <- results(psan_dseq, contrast = c("group", contrast_9), lfcThreshold = 1)
res_LFC_c9 <- lfcShrink(psan_dseq,
                        res = res_c9,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c9)
head(res_LFC_c9)
sum(res_LFC_c9$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 2,443

#10
contrast_10 <- c("0_Psan68_mock", "24_Psan68")
res_c10 <- results(psan_dseq, contrast = c("group", contrast_10))
sum(res_c10$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 550
res_c10 <- results(psan_dseq, contrast = c("group", contrast_10), lfcThreshold = 1)
res_LFC_c10 <- lfcShrink(psan_dseq,
                        res = res_c10,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c10)
head(res_LFC_c10)
sum(res_LFC_c10$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,534

#11
contrast_11 <- c("0_Psan68_mock", "48_Psan68")
res_c11 <- results(psan_dseq, contrast = c("group", contrast_11))
sum(res_c11$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 636
res_c11 <- results(psan_dseq, contrast = c("group", contrast_11), lfcThreshold = 1)
res_LFC_c11 <- lfcShrink(psan_dseq,
                        res = res_c11,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c11)
head(res_LFC_c11)
sum(res_LFC_c11$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,525

#12
contrast_12 <- c("0_Psan68_mock", "72_Psan68")
res_c12 <- results(psan_dseq, contrast = c("group", contrast_12))
sum(res_c12$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 791
res_c12 <- results(psan_dseq, contrast = c("group", contrast_12), lfcThreshold = 1)
res_LFC_c12 <- lfcShrink(psan_dseq,
                        res = res_c12,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c12)
head(res_LFC_c12)
sum(res_LFC_c12$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,905

#13
contrast_13 <- c("24_Psan10", "24_Psan65")
res_c13 <- results(psan_dseq, contrast = c("group", contrast_13))
sum(res_c13$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c13 <- results(psan_dseq, contrast = c("group", contrast_13), lfcThreshold = 1)
res_LFC_c13 <- lfcShrink(psan_dseq,
                        res = res_c13,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c13)
head(res_LFC_c13)
sum(res_LFC_c13$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#14
contrast_14 <- c("24_Psan10", "24_Psan68")
res_c14 <- results(psan_dseq, contrast = c("group", contrast_14))
sum(res_c14$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c14 <- results(psan_dseq, contrast = c("group", contrast_14), lfcThreshold = 1)
res_LFC_c14 <- lfcShrink(psan_dseq,
                        res = res_c14,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c14)
head(res_LFC_c14)
sum(res_LFC_c14$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#15
contrast_15 <- c("24_Psan65", "24_Psan68")
res_c15 <- results(psan_dseq, contrast = c("group", contrast_15))
sum(res_c15$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 1
res_c15 <- results(psan_dseq, contrast = c("group", contrast_15), lfcThreshold = 1)
res_LFC_c15 <- lfcShrink(psan_dseq,
                        res = res_c15,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c15)
head(res_LFC_c15)
sum(res_LFC_c15$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#16
contrast_16 <- c("48_Psan10", "48_Psan65")
res_c16 <- results(psan_dseq, contrast = c("group", contrast_16))
sum(res_c16$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c16 <- results(psan_dseq, contrast = c("group", contrast_16), lfcThreshold = 1)
res_LFC_c16 <- lfcShrink(psan_dseq,
                        res = res_c16,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c16)
head(res_LFC_c16)
sum(res_LFC_c16$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#17
contrast_17 <- c("48_Psan10", "48_Psan68")
res_c17 <- results(psan_dseq, contrast = c("group", contrast_17))
sum(res_c17$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c17 <- results(psan_dseq, contrast = c("group", contrast_17), lfcThreshold = 1)
res_LFC_c17 <- lfcShrink(psan_dseq,
                        res = res_c17,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c17)
head(res_LFC_c17)
sum(res_LFC_c17$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#18
contrast_18 <- c("48_Psan65", "48_Psan68")
res_c18 <- results(psan_dseq, contrast = c("group", contrast_18))
sum(res_c18$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c18 <- results(psan_dseq, contrast = c("group", contrast_18), lfcThreshold = 1)
res_LFC_c18 <- lfcShrink(psan_dseq,
                        res = res_c18,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c18)
head(res_LFC_c18)
sum(res_LFC_c18$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#19
contrast_19 <- c("72_Psan10", "72_Psan65")
res_c19 <- results(psan_dseq, contrast = c("group", contrast_19))
sum(res_c19$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c19 <- results(psan_dseq, contrast = c("group", contrast_19), lfcThreshold = 1)
res_LFC_c19 <- lfcShrink(psan_dseq,
                        res = res_c19,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c19)
head(res_LFC_c19)
sum(res_LFC_c19$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#20
contrast_20 <- c("72_Psan10", "72_Psan68")
res_c20 <- results(psan_dseq, contrast = c("group", contrast_20))
sum(res_c20$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c20 <- results(psan_dseq, contrast = c("group", contrast_20), lfcThreshold = 1)
res_LFC_c20 <- lfcShrink(psan_dseq,
                        res = res_c20,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c20)
head(res_LFC_c20)
sum(res_LFC_c20$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#21
contrast_21 <- c("72_Psan65", "72_Psan68")
res_c21 <- results(psan_dseq, contrast = c("group", contrast_21))
sum(res_c21$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c21 <- results(psan_dseq, contrast = c("group", contrast_21), lfcThreshold = 1)
res_LFC_c21 <- lfcShrink(psan_dseq,
                        res = res_c21,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c21)
head(res_LFC_c21)
sum(res_LFC_c21$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

#22
contrast_22 <- c("24_Psan10", "48_Psan10")
res_c22 <- results(psan_dseq, contrast = c("group", contrast_22))
sum(res_c22$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 318
res_c22 <- results(psan_dseq, contrast = c("group", contrast_22), lfcThreshold = 1)
res_LFC_c22 <- lfcShrink(psan_dseq,
                         res = res_c22,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c22)
head(res_LFC_c22)
sum(res_LFC_c22$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 462

#23
contrast_23 <- c("24_Psan10", "72_Psan10")
res_c23 <- results(psan_dseq, contrast = c("group", contrast_23))
sum(res_c23$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 122
res_c23 <- results(psan_dseq, contrast = c("group", contrast_23), lfcThreshold = 1)
res_LFC_c23 <- lfcShrink(psan_dseq,
                         res = res_c23,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c23)
head(res_LFC_c23)
sum(res_LFC_c23$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 694

#24
contrast_24 <- c("48_Psan10", "72_Psan10")
res_c24 <- results(psan_dseq, contrast = c("group", contrast_24))
sum(res_c24$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 0
res_c24 <- results(psan_dseq, contrast = c("group", contrast_24), lfcThreshold = 1)
res_LFC_c24 <- lfcShrink(psan_dseq,
                         res = res_c24,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c24)
head(res_LFC_c24)
sum(res_LFC_c24$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 11

#25
contrast_25 <- c("24_Psan65", "48_Psan65")
res_c25 <- results(psan_dseq, contrast = c("group", contrast_25))
sum(res_c25$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 148
res_c25 <- results(psan_dseq, contrast = c("group", contrast_25), lfcThreshold = 1)
res_LFC_c25 <- lfcShrink(psan_dseq,
                         res = res_c25,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c25)
head(res_LFC_c25)
sum(res_LFC_c25$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 105

#26
contrast_26 <- c("24_Psan65", "72_Psan65")
res_c26 <- results(psan_dseq, contrast = c("group", contrast_26))
sum(res_c26$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 56
res_c26 <- results(psan_dseq, contrast = c("group", contrast_26), lfcThreshold = 1)
res_LFC_c26 <- lfcShrink(psan_dseq,
                         res = res_c26,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c26)
head(res_LFC_c26)
sum(res_LFC_c26$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 296

#27
contrast_27 <- c("48_Psan65", "72_Psan65")
res_c27 <- results(psan_dseq, contrast = c("group", contrast_27))
sum(res_c27$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 2
res_c27 <- results(psan_dseq, contrast = c("group", contrast_27), lfcThreshold = 1)
res_LFC_c27 <- lfcShrink(psan_dseq,
                         res = res_c27,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c27)
head(res_LFC_c27)
sum(res_LFC_c27$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 49

#28
contrast_28 <- c("24_Psan68", "48_Psan68")
res_c28 <- results(psan_dseq, contrast = c("group", contrast_28))
sum(res_c28$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 135
res_c28 <- results(psan_dseq, contrast = c("group", contrast_28), lfcThreshold = 1)
res_LFC_c28 <- lfcShrink(psan_dseq,
                         res = res_c28,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c28)
head(res_LFC_c28)
sum(res_LFC_c28$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 74

#29
contrast_29 <- c("24_Psan68", "72_Psan68")
res_c29 <- results(psan_dseq, contrast = c("group", contrast_29))
sum(res_c29$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 227
res_c29 <- results(psan_dseq, contrast = c("group", contrast_29), lfcThreshold = 1)
res_LFC_c29 <- lfcShrink(psan_dseq,
                         res = res_c29,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c29)
head(res_LFC_c29)
sum(res_LFC_c29$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 249

#30
contrast_30 <- c("48_Psan68", "72_Psan68")
res_c30 <- results(psan_dseq, contrast = c("group", contrast_30))
sum(res_c30$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 1
res_c30 <- results(psan_dseq, contrast = c("group", contrast_30), lfcThreshold = 1)
res_LFC_c30 <- lfcShrink(psan_dseq,
                         res = res_c30,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c30)
head(res_LFC_c30)
sum(res_LFC_c30$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 0

###############################
resultsNames(psoj_dseq)

#31
contrast_31 <- c("0_Psoj_mock", "24_Psoj")
res_c31 <- results(psoj_dseq, contrast = c("group", contrast_31))
sum(res_c31$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 4,270
res_c31 <- results(psoj_dseq, contrast = c("group", contrast_31), lfcThreshold = 1)
res_LFC_c31 <- lfcShrink(psoj_dseq,
                        res = res_c31,
                        type = "ashr",
                        svalue = TRUE)
summary(res_LFC_c31)
head(res_LFC_c31)
sum(res_LFC_c31$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 10,368

#32
contrast_32 <- c("0_Psoj_mock", "48_Psoj")
res_c32 <- results(psoj_dseq, contrast = c("group", contrast_32))
sum(res_c32$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 3,294
res_c32 <- results(psoj_dseq, contrast = c("group", contrast_32), lfcThreshold = 1)
res_LFC_c32 <- lfcShrink(psoj_dseq,
                         res = res_c32,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c32)
head(res_LFC_c32)
sum(res_LFC_c32$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 7,772

#33
contrast_33 <- c("0_Psoj_mock", "72_Psoj")
res_c33 <- results(psoj_dseq, contrast = c("group", contrast_33))
sum(res_c33$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 2,588
res_c33 <- results(psoj_dseq, contrast = c("group", contrast_33), lfcThreshold = 1)
res_LFC_c33 <- lfcShrink(psoj_dseq,
                         res = res_c33,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c33)
head(res_LFC_c33)
sum(res_LFC_c33$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 6,217

#34
contrast_34 <- c("24_Psoj", "48_Psoj")
res_c34 <- results(psoj_dseq, contrast = c("group", contrast_34))
sum(res_c34$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 2,482
res_c34 <- results(psoj_dseq, contrast = c("group", contrast_34), lfcThreshold = 1)
res_LFC_c34 <- lfcShrink(psoj_dseq,
                         res = res_c34,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c34)
head(res_LFC_c34)
sum(res_LFC_c34$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 7,424

#35
contrast_35 <- c("24_Psoj", "72_Psoj")
res_c35 <- results(psoj_dseq, contrast = c("group", contrast_35))
sum(res_c35$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 2,501
res_c35 <- results(psoj_dseq, contrast = c("group", contrast_35), lfcThreshold = 1)
res_LFC_c35 <- lfcShrink(psoj_dseq,
                         res = res_c35,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c35)
head(res_LFC_c35)
sum(res_LFC_c35$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 8,880

#36
contrast_36 <- c("48_Psoj", "72_Psoj")
res_c36 <- results(psoj_dseq, contrast = c("group", contrast_36))
sum(res_c36$padj < 0.1, na.rm = TRUE)
#pajd < 0.1 = 258
res_c36 <- results(psoj_dseq, contrast = c("group", contrast_36), lfcThreshold = 1)
res_LFC_c36 <- lfcShrink(psoj_dseq,
                         res = res_c36,
                         type = "ashr",
                         svalue = TRUE)
summary(res_LFC_c36)
head(res_LFC_c36)
sum(res_LFC_c36$svalue < 0.1, na.rm = TRUE)
#svalue < 0.1 = 1,443

#################VISUALLY EXPLORING THE RESULTS - MA-PLOTS################
## MA-plots provide nice overview of results for each gene:
      #y-axis: count differences in terms of LFC between two groups
      #x-axis: mean counts across both groups
            #significantly differentially expressed genes in blue
#for 24/36 contrasts: res_cX, res_LFC_cX
#d <- res_c4
#d <- res_c5
#d <- res_c6
#d <- res_c7
#d <- res_c8
#d <- res_c9
#d <- res_c10
#d <- res_c11
#d <- res_c12
#d <- res_c22
#d <- res_c23
#d <- res_c24
#d <- res_c25
#d <- res_c26
#d <- res_c27
#d <- res_c28
#d <- res_c29
#d <- res_c30
#d <- res_c31
#d <- res_c32
#d <- res_c33
#d <- res_c34
#d <- res_c35
#d <- res_c36
d <- res_LFC_c4
#d <- res_LFC_c5
#d <- res_LFC_c6
#d <- res_LFC_c7
#d <- res_LFC_c8
#d <- res_LFC_c9
#d <- res_LFC_c10
#d <- res_LFC_c11
#d <- res_LFC_c12
#d <- res_LFC_c22
#d <- res_LFC_c23
#d <- res_LFC_c24
#d <- res_LFC_c25
#d <- res_LFC_c26
#d <- res_LFC_c27
#d <- res_LFC_c28
#d <- res_LFC_c29
#d <- res_LFC_c30
#d <- res_LFC_c31
#d <- res_LFC_c32
#d <- res_LFC_c33
#d <- res_LFC_c34
#d <- res_LFC_c35
#d <- res_LFC_c36
#d <- plotMA(d, ylim = c(-5, 5), returnData = TRUE)
#ggplot(d, aes(x = mean, y = lfc, color = isDE)) +
#  geom_point(size = 0.5) +
#  scale_x_log10() +
#  scale_y_continuous(limits = c(-10, 10)) +
#  scale_color_manual(values = c("grey50", "blue")) +
#  guides(color = FALSE) +
#  labs(x = "Mean of normalized counts",
#       y = "LFC")
#make plot of Shrunken LFC 
d <- plotMA(d, ylim = c(-5, 5), returnData = TRUE)
d$gene <- (rownames(res_LFC_c1))
d <- ggplot(d, aes(x = mean, y = lfc, color = isDE, text = gene)) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c("grey50", "blue")) +
  guides(color = FALSE) +
  labs(x = "Mean of normalized counts",
       y = "Shrunken LFC")
#make plot interactive with Plotly to show gene names when hovering
ggplotly(d, tooltip = "text")

## Plot specific genes
#select 5 genes with lowest adjusted p-value
top5 <- row.names()

g10813
g19264






### re-run with p = 0.05 not p = 0.1