## Script to prep genome metadata for GO analysis
## - We need a dataframe with gene lengths
## - We need a dataframe with gene-to-GOterm mappings (two columns only, one row per mapping)

# SETUP ------------------------------------------------------------------------
## Load packages
library(ape)
library(tidyverse)

## Input files
gff_file <- "/fs/project/PAS0471/linda/results/genome_annotation/braker2/P68_no_mask/braker/augustus.hints.gff3"
GO_annot_file <- "/fs/project/PAS0471/jelmer/assist/2021-12_linda/results/entap/p65/final_results/final_annotations_lvl0_enrich_geneid_go.tsv"

## Output files
outdir <- "results/GO"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
gene_len_file <- file.path(outdir, "gene_lens.txt")
GO_map_file <- file.path(outdir, "GO_map.txt")

## Read the GFF file
gff <- read.gff(gff_file)


# GENE LENGTHS -----------------------------------------------------------------
## Create a df with just gene IDs and gene lengths
gene_len_df <- gff %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub(".*ID=([^;]+).*", "\\1", attributes),
         length = end - start + 1) %>%
  select(gene_id, length) %>%
  arrange(gene_id)
message("Nr of genes in gene length df: ", nrow(gene_len_df))

## Write the df to file
write_tsv(gene_len_df, gene_len_file)


# GO MAP -----------------------------------------------------------------------
GO_map <- read_tsv(GO_annot_file) %>%      # Should have columns "gene_id" and "go_term"
  arrange(gene_id) %>%
  mutate(gene_id = sub(".t\\d+", "", gene_id)) %>%   # Turn transcript IDs into gene IDs
  distinct() %>%     # Get rid of duplicates (different transcripts for same gene annotated with same GO term)
  as.data.frame()
message("Nr of genes in gene-to-term mappings in the GO map: ", nrow(GO_map))
message("Nr of unique genes in the GO map: ", length(unique((GO_map$gene_id))))
message("Nr of unique GO terms in the GO map: ", length(unique((GO_map$go_term))))

## Write the df to file
write_tsv(GO_map, GO_map_file)
