# Volcano plot for the differential expression analysis results
library(ggplot2)
library(ggrepel)
library(dbplyr)


# Load all genes
all_genes <- read_csv("results/deseq2_results.csv")
top_genes <- all_genes %>%
  arrange(padj) %>%
  head(20)

ggplot(all_genes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = abs(log2FoldChange)), alpha = 0.6, size = 2) +  # Color by fold change
  scale_color_gradient(low = "blue", high = "red", name = "Absolute Log2 Fold Change") +  # Gradient colors
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_genes,  
    aes(label = ...1),  
    box.padding = 0.5, 
    point.padding = 0.2,
    color = "black",  
    size = 2, 
    max.overlaps = 19
  ) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = "Volcano Plot of Differential Expression",
    color = "Absolute Log2 Fold Change"
  )
