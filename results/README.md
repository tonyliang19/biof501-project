# Results

This directory store results of the pipeline, using default parameters for input:

```bash
nextflow run main.nf --samplesheet data/samplesheet.csv --outdir results -profile docker,test
```

Then its expected directory structure would be:

```bash
results/
├── README.md
├── deseq2
│   └── deseq2_result.csv
├── enhanced_volcano
│   └── volcano_plot.png
├── fastqc
│   ├── SRR29891674
│   │   ├── SRR29891674_sample_1_fastqc.html
│   │   ├── SRR29891674_sample_1_fastqc.zip
│   │   ├── SRR29891674_sample_2_fastqc.html
│   │   └── SRR29891674_sample_2_fastqc.zip
│   ├── SRR29891675
│   │   ├── SRR29891675_sample_1_fastqc.html
│   │   ├── SRR29891675_sample_1_fastqc.zip
│   │   ├── SRR29891675_sample_2_fastqc.html
│   │   └── SRR29891675_sample_2_fastqc.zip
│   ├── SRR29891677
│   │   ├── SRR29891677_sample_1_fastqc.html
│   │   ├── SRR29891677_sample_1_fastqc.zip
│   │   ├── SRR29891677_sample_2_fastqc.html
│   │   └── SRR29891677_sample_2_fastqc.zip
│   └── SRR29891678
│       ├── SRR29891678_sample_1_fastqc.html
│       ├── SRR29891678_sample_1_fastqc.zip
│       ├── SRR29891678_sample_2_fastqc.html
│       └── SRR29891678_sample_2_fastqc.zip
├── feature_counts
│   ├── feature_counts.log
│   └── feature_counts.rds
├── hisat2_align
│   ├── SRR29891674
│   │   ├── SRR29891674.sam
│   │   └── hisat2-SRR29891674.log
│   ├── SRR29891675
│   │   ├── SRR29891675.sam
│   │   └── hisat2-SRR29891675.log
│   ├── SRR29891677
│   │   ├── SRR29891677.sam
│   │   └── hisat2-SRR29891677.log
│   └── SRR29891678
│       ├── SRR29891678.sam
│       └── hisat2-SRR29891678.log
├── hisat2_build
│   └── hisat2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.1.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.2.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.3.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.4.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.5.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.6.ht2
│       ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa.7.ht2
│       └── Homo_sapiens.GRCh38.dna.chromosome.1.fa.8.ht2
├── map_ensembl_id
│   └── mapped_id_deseq_result.csv
├── pipeline_info
│   ├── dge_analysis_versions.yml
│   ├── execution_report_2024-12-03_17-26-38.html
│   ├── execution_timeline_2024-12-03_17-26-38.html
│   ├── execution_trace_2024-12-03_17-26-38.txt
│   └── pipeline_dag_2024-12-03_17-26-38.html
├── samtools_sort
│   ├── SRR29891674_sorted.bam
│   ├── SRR29891675_sorted.bam
│   ├── SRR29891677_sorted.bam
│   └── SRR29891678_sorted.bam
├── samtools_to_bam
│   ├── SRR29891674.bam
│   ├── SRR29891675.bam
│   ├── SRR29891677.bam
│   └── SRR29891678.bam
└── trimgalore
    ├── SRR29891674
    │   ├── SRR29891674_sample_1_val_1.fq.gz
    │   └── SRR29891674_sample_2_val_2.fq.gz
    ├── SRR29891675
    │   ├── SRR29891675_sample_1_val_1.fq.gz
    │   └── SRR29891675_sample_2_val_2.fq.gz
    ├── SRR29891677
    │   ├── SRR29891677_sample_1_val_1.fq.gz
    │   └── SRR29891677_sample_2_val_2.fq.gz
    └── SRR29891678
        ├── SRR29891678_sample_1_val_1.fq.gz
        └── SRR29891678_sample_2_val_2.fq.gz
```

Due to size limit of files tracked in git, this result directory only contain small output files on github, and this README.