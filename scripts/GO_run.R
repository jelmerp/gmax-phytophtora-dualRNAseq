# SETUP ------------------------------------------------------------------------
## Load packages
library(tidyverse)
library(goseq)

## Source script with functions
source("mcic-scripts/rnaseq/rfuns/GO_fun.R") #! MAKE SURE TO UPDATE (GIT PULL OR SIMILAR): I MADE SOME CHANGES

## Input files
gene_len_file <- "results/GO/gene_lens.txt"
GO_map_file <- "results/GO/GO_map.txt"
DE_file <- "/fs/project/PAS0471/linda/results/RNA-Seq/DESeq2/pairwise_results_all/0_Psan10_vs_0_Psan65_all-res.txt"

## Output files
# ...

## Read the GO map and gene lengths
gene_len_df <- read_tsv(gene_len_file)          # Should have columns "gene_id" and "length"
GO_map <- read.delim(GO_map_file, sep = "\t")   # Should have columns "gene_id" and "go_term"
#(GO_map needs to be a regular df and not a tibble)

## Read the DE results
DE_dir <- "../../../linda/results/RNA-Seq/DESeq2/pairwise_results_all"
DE_files <- list.files(DE_dir, full.names = TRUE)
col_names <- c("gene_id", read_tsv(DE_files[1]) %>% colnames())  # These file have rownames so one column name is missing...
DE <- map_dfr(.x = DE_files, .f = read_tsv, skip = 1, col_names = col_names) %>%
  mutate(contrast = paste0(level1, "_vs_", level2))


# USE THE WRAPPER FUNCTION -----------------------------------------------------
## Run the GO analysis for 1 contrast
GO_wrap(contrast_id = DE$contrast[1],
        DE_res = DE, GO_map = GO_map, gene_lens = gene_len_df)

## Run the GO analysis for all contrasts
map_dfr(.x = unique(DE$contrast),
        .f = GO_wrap,
        DE_res = DE, GO_map = GO_map, gene_lens = gene_len_df)


# STEP-BY_STEP GO ANALYSIS FOR ONE CONTRAST ------------------------------------
## Settings

contrast_id <- "0_Psoj_vs_24_Psoj"

contrast_id <- DE$contrast[1]      # Take the first contrast by means of example
DE_direction <- "either"           # 'either' takes DE results as is, 'up'/'down' only considers LFC >0 or < 0
                                   # 'both' tests each direction separately
p_DE <- 0.05                       # Adjusted p-value threshold (default in the functions is also 0.05)
lfc_DE <- 0                        # LFC threshold (default in the functions is also 0, so no threshold)
min_in_cat <- 2                    # Min. nr. of total terms in GO category (exlude very small categories)
max_in_cat <- Inf                  # Min. nr. of total terms in GO category (exclude very large categories)
min_DE_in_cat <- 2                 # Min. nr. DE genes in GO term for a term to be significant
                                   #! I EXPLAINED THIS WRONG -- CF. `min_in_cat` AND `min_DE_in_cat`
ontologies <- c("BP", "MF", "CC")  # Ontologies to include

## Create a named vector with DE results: 0 is non-significant, 1 is significant
fDE <- DE %>%
  filter(contrast == contrast_id,
         !is.na(padj)) %>%   # Exclude genes with NA adj-p-val - not tested
  mutate(sig = ifelse(padj < p_DE & abs(log2FoldChange) > lfc_DE, 1, 0)) %>% 
  arrange(gene_id)
DE_vec <- fDE$sig
names(DE_vec) <- fDE$gene_id

## Remove rows from gene length df not in the DE_vec
fgene_lens <- gene_len_df %>% filter(gene_id %in% names(DE_vec))

## Remove elements from DE_vec not among the gene lengths
fDE_vec <- DE_vec[names(DE_vec) %in% fgene_lens$gene_id]

## Check that gene lengths and contrast vector contain the same genes in the same order
stopifnot(all(fgene_lens$gene_id == names(fDE_vec)))

## goseq probability weighting function based on gene lengths
pwf <- nullp(DEgenes = fDE_vec,
             bias.data = fgene_lens$length,
             plot.fit = FALSE)

## Run GO test with main goseq function
GO_df <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")

## Process GO results
GO_df <- GO_df %>%
  filter(numDEInCat > 0,                 # P-adjustment only for genes that were actually tested
         numInCat >= min_in_cat,
         numInCat <= max_in_cat,
         ontology %in% ontologies)
if (filter_no_descrip == TRUE) GO_df <- GO_df %>% filter(!is.na(description))
GO_df <- GO_df %>%
  # Note: we need to do the p-adjustment ourselves, goseq does not do this
  mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
         sig = ifelse(padj < 0.05 & numDEInCat >= min_DE_in_cat, 1, 0),
         contrast = contrast_id,
         DE_direction = DE_direction) %>%
  select(contrast, DE_direction,
         sig, p = over_represented_pvalue, padj,
         numDEInCat, numInCat,
         category, ontology, description = term)

cat("Contrast:", contrast_id,
    "    DE direction:", DE_direction,
    "    Nr GO cats:", nrow(GO_df),
    "    Nr DE:", sum(DE_vec),
    "    Nr sign. GO:", sum(GO_df$sig), "\n")
