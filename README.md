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


## Sample data

Using the [GSE272659 study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272659) and looking its 2 samples only:

1. [GSM8408768 NT-ctrl-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408768), which corresponds to run [SRR29891678](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891678&display=metadata)
2. [GSM8408771 NT-Doxo-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408771), which corresponds to run [SRR29891675](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891675&display=metadata)

Construct a metadata for your samples like the following:

```csv
sample_name,condition
SRR29891678,control
SRR29891675,treatment
```

The `sample_name` here is the run name.

TODO: this sample data section need to be rewritten, since the metadata part can be inside the samplesheet
