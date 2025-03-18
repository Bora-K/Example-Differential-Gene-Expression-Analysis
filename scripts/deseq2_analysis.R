# Load the necessary libraries
library(tidyverse)
library(GEOquery)
library(dplyr)
library(DESeq2)
library(purrr)
library(stringr)


# Download the GEO dataset
gse <- getGEO(GEO = "GSE223954", GSEMatrix = TRUE)

# Untar the raw counts
untar("data/GSE223954_RAW.tar", exdir = "data")

# Get the filenames
files <- list.files("data", pattern = "\\.gz$", full.names = TRUE)

# Function to read a single file and rename the "Count" column to the GSM ID
# Function to read and rename "Count" to GSM ID (using base R renaming)
read_counts <- function(file) {
  gsm_id <- str_extract(basename(file), "^GSM\\d+")  # Extract GSM ID
  
  # Read the file and rename "Count" using base R (no dplyr::rename)
  df <- read.delim(file, check.names = FALSE) %>%
    dplyr::select(Id, Symbol, Count)
  
  # Rename "Count" column to GSM ID using base R
  colnames(df)[colnames(df) == "Count"] <- gsm_id
  return(df)
}

# Read and merge all files
count_matrix <- files %>%
  purrr::map(read_counts) %>%
  purrr::reduce(full_join, by = c("Id", "Symbol"))

# Set row names to gene/transcript IDs
rownames(count_matrix) <- count_matrix$Id
count_matrix <- count_matrix[, !(colnames(count_matrix) %in% c("Id", "Symbol"))]

# Verify the cleaned matrix
head(count_matrix)

# Get the condition name
name_condition <- "characteristics_ch1.3"

# Get metadata
metadata <- pData(phenoData(gse[[1]])) %>%
  rownames_to_column("GSM_ID") %>%
  dplyr::select(GSM_ID, condition = name_condition)

metadata$condition <- str_remove(metadata$condition, "condition: ")
metadata$condition <- as.factor(metadata$condition)

# Check alignment
all(colnames(count_matrix) %in% metadata$GSM_ID)  # Should return TRUE

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

# Filter and run DESeq2 (as before)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

# See results
res <- results(dds, name = "condition_shFGFR1_vs_Control", alpha = 0.05)

significant_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(significant_genes[order(significant_genes$padj), ], 10)  # Top 10 genes

# Save the results
write.csv(res, "results/deseq2_results.csv", row.names = TRUE)
write.csv(significant_genes, "results/significant_genes.csv", row.names = FALSE)
