#!/usr/bin/env Rscript

# Use this script to map ensembl gene id to gene symbols

# Load libraries
library(org.Hs.eg.db)

# Parse cli
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
# Load inputs
deseq_df <- read.csv(data_path)
# Get the annotations and drop those that do not have a corresponding 
annots <- select(org.Hs.eg.db, keys=deseq_df$ensembl_id, columns="SYMBOL", keytype="ENSEMBL_ID")
merged <- merge(deseq_df, annots, by.x="ensembl_id", by.y="ENSEMBL_ID")
mapped_df <- merged[!is.na(merged$SYMBOL), ]

message("\nDropped ", nrow(mapped_df), " genes that did not have gene symbol from total ", nrow(merged), " genes")


# lastly write out to csv
outfile <- "mapped_id_deseq_result.csv"
write.csv(mapped_df, file=outfile, row.names=FALSE)




