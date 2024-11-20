# Biof 501 Special Topics in Bioinformatics Project

Differential Gene Expression Analysis

Author: Tony Liang

---

**Table of Contents**:

1. [Project Outline](#project-outline)
1. [Repository Structure](#directory-contents)
1. [Setup](#setup)
1. [Preparing Input](#preparing-input)
1. [Running Instructions](#running-instruction)

## Repository Structure

The following ilustrates the tree structure of this repository, those that ends with `/` symbol are directories, otherwise it is treated as a regular file. 

This structure is adapted from **nf-core** standard worfklow from the [rnaseq](https://github.com/nf-core/rnaseq/tree/master) repository.

```bash
./
├── Makefile
├── README.md
├── bin/
│   └── featureCounts.R*
├── conf/
│   ├── base.config
│   └── test.config
├── data/
│   ├── SRR29891674_sample_1.fastq.gz
│   ├── SRR29891674_sample_2.fastq.gz
│   ├── SRR29891675_sample_1.fastq.gz
│   ├── SRR29891675_sample_2.fastq.gz
│   ├── SRR29891677_sample_1.fastq.gz
│   ├── SRR29891677_sample_2.fastq.gz
│   ├── SRR29891678_sample_1.fastq.gz
│   ├── SRR29891678_sample_2.fastq.gz
│   └── samplesheet.csv
├── main.nf
├── modules/
│   ├── local/
│   │   ├── download_reference/
│   │   │   └── main.nf
│   │   ├── fastqc/
│   │   │   └── main.nf
│   │   ├── feature_counts/
│   │   │   └── main.nf
│   │   ├── hisat2/
│   │   │   ├── hisat2_align/
│   │   │   └── hisat2_build/
│   │   ├── samtools/
│   │   │   ├── samtools_sort/
│   │   │   └── samtools_to_bam/
│   │   └── trimgalore/
│   │       └── main.nf
│   └── nf-core/
│       └── main.nf
├── nextflow.config
└── run_nxf.sh*
```

### Directories

- `bin`: This directory contains **executable** binary scripts used in the processes, provided if you have the relevant binary installed like `Rscript` or `python3`. For more information please see [here](bin/README.md).
- `modules`: This directory contains nextflow processes definitions following these conventions:
    - `local/<tool_name>/main.nf`: `<tool_name>` would be an actual binary tool that could be treated as smallest unit software to use in the workflow.
    - `nf-core/main.nf`: Code adapted from nf-core
- `data`:
- `conf`:




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

Using the [GSE272659 study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272659) and looking its 4 samples only:

1. [GSM8408768 NT-Ctrl-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408768), which corresponds to run [SRR29891678](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891678&display=metadata)
2. [GSM8408769 NT-Ctrl-2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408769), which corresponds to run [SRR29891677](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891677&display=metadata)
3. [GSM8408771 NT-Doxo-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408771), which corresponds to run [SRR29891675](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891675&display=metadata)
4. [GSM8408772 NT-Doxo-2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408772), which corresponds to run [SRR29891674](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891674&display=metadata)



## Setup

There are some software dependencies in order to run this pipeline, corresponding softwares and versions are tested in the following:

```bash
# <sofware>-<version_number>
bash-5.1.16
jave-11.0.22
Docker--20.10.23
git-2.34.1
nextflow-24.04.4
```

Other tools used in the pipeline are from public repository docker images retrieved mostly from the [biocontainers organization](https://quay.io/organization/biocontainers)

For a list of software images used in this pipeline see the configuration file [here](nextflow.config) under the `process` scope.

### Preparing Input

> [!NOTE]
> This instruction is yet not completed, still under construction

First, prepare a samplesheet with your input data that looks as follows:

**data/samplesheet.csv**:

```csv
sample_name,fastq1,fastq2
control_rep1,data/SRR29891678_sample_1.fastq.gz,data/SRR29891678_sample_2.fastq.gz
control_rep2,data/SRR29891677_sample_1.fastq.gz,data/SRR29891677_sample_2.fastq.gz
treatment_rep1,data/SRR29891675_sample_1.fastq.gz,data/SRR29891675_sample_2.fastq.gz
treatment_rep2,data/SRR29891674_sample_1.fastq.gz,data/SRR29891674_sample_2.fastq.gz
```

The `sample_name` here ilustrates condition and replicate number together


Each row represents a pair of fastq files (paired end) corresponding to some condition + replicate.

TODO: this sample data section need to be rewritten, since the metadata part can be inside the samplesheet

### Running Instruction


First, clone the this github repository and change the working directory to the clone repo using:

```bash
# Assuming you use one of bash/zsh or other unix systems
# NO Powershell
git clone https://github.com/tonyliang19/biof501-project.git
cd biof501-project/
```

Then the pipeline could be run as the following using a test data contained already in this repository:

```bash
# The outdir could also be replace to some other directory of your preference
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --outdir results \
    -profile docker,test
```

Overall the running command should follow this structure:

```bash
nextflow run main.nf \
    --samplesheet <SOME_SAMPLESHEET_CSV> \
    --outdir <OUTDIR> \
    -profile docker,test
```

where `<SOME_SAMPLESHEET_CSV>` is the csv data that follows format in [preparing-input section](#preparing-input) and `OUTDIR` being the directory you want the output files to store.

> [!NOTE]
> The pipeline should take some time run for the very first time, because of the containerized images it have to pull from internet, and downloads of some reference data
>
> Also the profile there should be NO SPACE, i.e.:
>
> -profile docker,test works
>
> -profile docker, test DO NOT WORK


## Reference

Nextflow related

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

    > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

    > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

Pipeline tools

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [featureCounts](https://pubmed.ncbi.nlm.nih.gov/24227677/)

  > Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PubMed PMID: 24227677.
