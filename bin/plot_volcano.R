#!/usr/bin/env Rscript

library(EnhancedVolcano)
library(ggplot2)

# first parse cli, this gives as string for every arg
# So args is just character vector
args <- commandArgs(trailingOnly = TRUE)
deseq_result_path <- args[1] # This is deseq result with mapped gene symbol
# The underlying R version is old, so native pipe |> do not support

# Load the data first
plot_df <- read.csv(deseq_result_path)
# Plotting parameters (Could be abstracted to nextflow config)
pCutOff <- 1e-05
FCcutoff <- 0.5
pointSize <- 3.0
labSize <- 6.0
# Then plot it and save to memory
# x and y refers to actual columns in the plot df

volcano_plot <- EnhancedVolcano(plot_df, lab=final_df$SYMBOL,
                x= "log2FoldChange",
                y = "pvalue",
                pCutoff = pCutOff,
                FCcutoff = FCcutoff,
                pointSize = pointSize ,
                labSize = labSize)

# Lastly save it using ggsave
filename <- "volcano_plot.png"
ggsave(filename=filename, plot=volcano_plot)
