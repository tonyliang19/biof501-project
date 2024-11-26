#!/usr/bin/env Rscript

library(Rsubread)

# first parse cli, this gives as string for every arg
# So args is just character vector
args <- commandArgs(trailingOnly = TRUE)
# BAM could be more than 1 files, these are separated by space
bam_list <- args[1]
# The underlying R version is old, so native pipe |> do not support
bam <- unlist(strsplit(bam_list, split = " "))
gtf <- args[2]

# Generates count matrix
fc <- featureCounts(files = bam, annot.ext = gtf, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,
GTF.featureType = "gene")

# Also rename those column names of counts to remove the _sorted.bam suffix
pattern <- "_sorted.bam"
colnames(fc$counts) <- gsub(pattern, "", colnames(fc$counts))

# Lastly write it to file
file <- "feature_counts.rds"
saveRDS(fc, file = file)

# And print the first 6 rows of the counts data to see matching and targets
cat("\nTargets are: ", fc$targets)

# And return the head of counts object in log
cat("\n")
print(head(fc$counts))


