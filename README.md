# Biof 501 Special Topics in Bioinformatics Project

Differential Gene Expression Analysis

Author: Tony Liang

## Project Outline


**Question**: Identify genes that are differential expressed between different conditions like disease vs control

**Data input**: Fastq file of RNA-seq data, single or paired allowed

**Output**: Summary report that contains visualization of results and patterns observed from differential expression analysis

### Workflow stages/steps

1. Read data from fastq format
2. Quality check and filtering
    -  Apply initial QC
    - Trim reads
3. Read alignment or quantifcation, with or without reference genome
4. Gene expression quantification to generate count matrix
5. Normalization of count data
6. Perform differential expression analysis
7. Visualization and summary report

