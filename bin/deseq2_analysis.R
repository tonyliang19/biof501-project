#!/usr/bin/env Rscript

# Use this script to run deseq analysis on feature counts data

# Load libraries
library(DESeq2)
# first parse cli, this gives as string for every arg
# So args is just character vector
args <- commandArgs(trailingOnly = TRUE)

# Args are in this order:
# Rds, metadata.csv
rds_path <- args[1] # Full featureCounts object, only using counts there for now
metadata_path <- args[2] # contains experiment design data
# Read inputs
fc <- readRDS(rds_path)
# Extract the raw counts directly from fc 
raw_counts <- fc$counts
# Handle metadata here
metadata_df <- read.csv(metadata_path)
# Assign sample identifier as rownames
sample_identifier <- "sample_name"
rownames(metadata_df) <- metadata_df[[sample_identifier]]

# Then make sure the condition of metadata is factor
metadata_df$condition <- as.factor(metadata_df$condition)

# Resort our counts based on the metadata and make sure order match
raw_counts <- raw_counts[, rownames(metadata_df)]
if (!all(rownames(metadata_df) == colnames(raw_counts))) stop("Order mismatched between samples")

# =============================
# Run DESEQ here
# =============================

# First construct DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata_df,
  design = ~ condition
  )

dds <- DESeq(dds, betaPrior=FALSE)
# Add the contrast 
contrast <- c("condition", "control", "treatment")
# Then extract results and shrink its log fold changes
res <- results(dds, contrast = contrast)
res <- lfcShrink(dds, contrast = contrast, res=res, type = 'normal')

# Get the relevant columns for plotting later and write them to file
output_df <- data.frame(
  log2FoldChange  = res$log2FoldChange,
  pvalue          = res$pvalue, 
  ensembl_id      = rownames(res)
)

# Remove the ones that do not have log2FC
output_df <- output_df[!is.na(output_df$log2FoldChange), ]

# Modify this filename
outfile <- "deseq2_result.csv"

write.csv(output_df, file = outfile, row.names=FALSE)

message("\nSaved result to ", outfile)