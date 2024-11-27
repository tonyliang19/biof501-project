# Data

## Initial directory Structure

This directory stores a raw sample data used in the pipeline for testing purposes. 

The initial file structure is the following:

```bash
data/
├── README.md
├── SRR29891674_sample_1.fastq.gz
├── SRR29891674_sample_2.fastq.gz
├── SRR29891675_sample_1.fastq.gz
├── SRR29891675_sample_2.fastq.gz
├── SRR29891677_sample_1.fastq.gz
├── SRR29891677_sample_2.fastq.gz
├── SRR29891678_sample_1.fastq.gz
├── SRR29891678_sample_2.fastq.gz
└── samplesheet.csv
```

The `samplesheet.csv` is the required input of the pipeline, it takes the following format:

```csv
sample_name,condition,rep,fastq1,fastq2
<sample_name>,<condition>,<repNumber>,<data/path_to_sample_name_fastq1.gz>,<data/path_to_sample_name_fastq2.gz>
<sample_name>,<condition>,<repNumber>,<data/path_to_sample_name_fastq1.gz>,<data/path_to_sample_name_fastq2.gz>
```

The `<sample_name>` is unique identifier tell which sample is which one, `<condition>` is experiment metadata, telling if sample belongs to certain experiment condition like control or treatment, `<repNumber>` is number of biological replicate, `<data/path_to_sample_name_fastq1.gz>` and `<data/ path_to_sample_name_fastq2.gz>` are corresponding paths to a paired-end sequencing data in `fastq.gz` format.


## Sample data

The `SRR*.fastq.gz` are sampled fastq files from their respective `SRR` run identifier. Due to its large size, we sampled 10,000 reads randomly from each sample from the [GSE272659 study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272659) and looking these 4 samples only:

1. [GSM8408768 NT-Ctrl-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408768), which corresponds to run [SRR29891678](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891678&display=metadata)
2. [GSM8408769 NT-Ctrl-2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408769), which corresponds to run [SRR29891677](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891677&display=metadata)
3. [GSM8408771 NT-Doxo-1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408771), which corresponds to run [SRR29891675](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891675&display=metadata)
4. [GSM8408772 NT-Doxo-2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8408772), which corresponds to run [SRR29891674](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR29891674&display=metadata)

The [SRA-toolkit-3.1.1](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) was used to download raw fastq files of the sequencing runs, then [Seqtk-1.4-r132-dirty](https://github.com/lh3/seqtk) is used to take random 10000 reads from each sample of the study.

To reproduce this, you could follow these steps assuming you have the relevant tools installed:

1. Download accession run file and its dependencies using `prefetch-3.1.1` from the toolkit:

```bash
# accession_number would be like of SRR29891678
# This would create a directory <accession_number> 
prefetch <accession_number>
```

2. Extract raw fastq files from the downloaded files from step 1

```bash
# This would extract fastq files out from the <accession_number> directory
fasterq-dump <accession_number>
```

3. Sample 10000 read pairs from the extracted fastq files

```bash
# This requires seqtk in your PATH variable
# The -s100 option is to use seed 100 to reproduce the sampling result
# The 10000 option is number of read pairs to take sample
# The <accession_number>/<accession_number> refers to fastq file in its run directory from prefetch + fasterq-dump
# This commands runs on both forward fastq and reverse fastq i.e. _1 and _2 
seqtk sample -s100 <accession_number>/<accession_number>_1.fastq 10000 > <accession_number>/<accession_number>_sample_1.fastq
seqtk sample -s100 <accession_number>/<accession_number>_2.fastq 10000 > <accession_number>/<accession_number>_sample_2.fastq
```

These instructions are based on official guide of [sra-toolkit](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) and [seqtk](https://github.com/lh3/seqtk/blob/master/README.md) respectively.

## Intermediate files from pipeline

This pipeline focus on differential gene expression analysis on given fastq files. Hence, it requires some other genome and genome annotation files, including sample metadata information to setup experiment design for downstream analyses. 

The metadata is also included in the `samplesheet.csv`, and its separated automatically as a different file during pipeline computation as `metadata.csv` with the following format:

```csv
sample_name,condition,rep
SRR29891678,control,rep1
SRR29891677,control,rep2
SRR29891675,treatment,rep1
SRR29891674,treatment,rep2
```

Then, for all other required files to complete the pipeline, those are downloaded during computation from the internet due to its nature of large size.

The following is the final directory structure of this data directory after execution of the pipeline:

```bash
data/
├── Homo_sapiens.GRCh38.113.gtf.gz
├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
├── Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
├── README.md
├── SRR29891674_sample_1.fastq.gz
├── SRR29891674_sample_2.fastq.gz
├── ...
├── SRR29891678_sample_1.fastq.gz
├── SRR29891678_sample_2.fastq.gz
├── metadata.csv
└── samplesheet.csv
```
