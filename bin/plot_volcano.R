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
rownames(plot_df) <- plot_df$SYMBOL

# Plotting parameters (Could be abstracted to nextflow config)
# Given this runs on biologically not significant testing dataset
# So adjust pCutOff accordingly
# Default starting values
pCutOff <- 0.05
FCcutoff <- 0.5
# Initial filtering of significant points
significant <- abs(plot_df$log2FoldChange) > FCcutoff & plot_df$pvalue < pCutOff
num_significant <- sum(significant)

# Adjust thresholds if there isnt too many significant points
if (num_significant < 5) {
  cat("Few significant points found. Relax thresholds.\n")
  pCutOff <- 0.2  # Relax p-value threshold
  FCcutoff <- FCcutoff - 0.4 # Relax fold-change threshold
} else if (num_significant > 50) {
  cat("Many significant points found. Tighten thresholds.\n")
  pCutOff <- 0.01  # Tighten p-value threshold
  FCcutoff <- FCcutoff + 0.4  # Tighten fold-change threshold
}

# Graphic parameters
pointSize <- 3.0
labSize <- 6.0
# Then plot it and save to memory
# x and y refers to actual columns in the plot df

volcano_plot <- EnhancedVolcano(
  toptable=plot_df, 
  lab=rownames(plot_df),
  x= "log2FoldChange",
  y = "pvalue",
  pCutoff = pCutOff,
  FCcutoff = FCcutoff,
  pointSize = pointSize,
  labSize = labSize
)

# Lastly save it using ggsave
filename <- "volcano_plot.png"

message("\nSaving volcanot plot to ", filename)
ggsave(filename=filename, plot=volcano_plot)
