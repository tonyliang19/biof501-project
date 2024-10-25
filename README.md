# Biof 501 Special Topics in Bioinformatics Project

Differential Gene Expression Analysis

Author: Tony Liang

## Project Outline


**Question**: Identify genes from 10 genes of interest that are differential expressed between different conditions like disease vs control

**Data input**: Fastq file of RNA-seq data, single or paired allowed

**Output**: Summary report that contains visualization of results and patterns observed from differential expression analysis

### Workflow stages/steps

Take  data from fastq format that has disease and control and perform the following:
1. Quality check and filtering
    -  Apply initial QC using **fastqc**
    - Trim reads **cutadapt** or ***trimmmomatic**
3. Read alignment reference genome (containing 10 genes of interest only) **STAR** or **hisat2**
   - Gene expression quantification to generate count matrix
6. Perform differential expression analysis in R using **DESeq2**
7. Visualization and summary report in R using **pheatmap**

