library(clusterProfiler)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)

# Filter significant genes
significant_genes <- read_csv("results/significant_genes.csv")

# Remove version numbers from NM codes
nm_codes <- significant_genes$...1
#nm_codes <- sub("\\..*", "", nm_codes)

# Map NM codes to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = nm_codes,
  keytype = "REFSEQ",
  column = "ENTREZID"
)

# Add Entrez IDs to the significant genes dataframe
significant_genes$EntrezID <- entrez_ids

# Perform KEGG pathway enrichment
kegg_enrichment <- enrichKEGG(
  gene = significant_genes$EntrezID,
  organism = "hsa", # Human
)

# View the results
head(kegg_enrichment)

# Visualize KEGG enrichment results
dotplot(kegg_enrichment, showCategory = 20)  # Show top 20 enriched pathways
