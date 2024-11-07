/*
    This is the main entrance of the workflow for differential
    gene expression analysis
*/

// Including modules
include { FASTQC }                  from "./modules/local/fastqc"
include { HISAT2_BUILD }            from "./modules/local/hisat2/hisat2_build"
include { HISAT2_ALIGN }            from "./modules/local/hisat2/hisat2_align"
include { DOWNLOAD_REFERENCE }      from "./modules/local/download_reference"
include { TRIMGALORE }              from "./modules/local/trimgalore"
include { SAMTOOLS_TO_BAM }         from "./modules/local/samtools/samtools_to_bam"
include { SAMTOOLS_SORT }           from "./modules/local/samtools/samtools_sort"
include { FEATURE_COUNTS }          from "./modules/local/feature_counts"
include { softwareVersionsToYAML }  from "./modules/nf-core/main.nf"


workflow {
    println "Hello World"
    /*
        Should take in a samplesheet.csv where each entry is geo accession code
        then have process to download fastq file from the accesion code

        Say format is:

        accession_code,some_other_column,...
        12345678,some_other_value,...
    */

    // Initialize version file to store
    ch_versions = Channel.empty()

    // This channel reads the samplesheet, and construct path of the records
        Channel
        .fromPath( params.samplesheet )
        .splitCsv( header: true )
        // This gives [ sample_name, [ name, condition ], [ reads ] ]
        .map { row ->
            // Splits sample name from condition_rep
            def split_name = row.sample_name.split("_")
            // Get run name from the fastq file
            def run_name = file(row.fastq1.split("_1.fastq.gz")[0]).getSimpleName()
            [run_name , split_name[0], split_name[1], [ file(row.fastq1), file(row.fastq2) ] ]
        }
        .set { raw_record }
    
    // Use only these relevant values for process running
    raw_record
        .map { run_name, condition, rep, reads -> [ run_name, reads ] }
        .set { record }

    // Then extract the metadata to save it as csv for later use
    // Need extra row as header of columns
    header_row = Channel.fromList(['run_name,condition,rep'])
    // Then concat it with the relevant metadata
    header_row
        .concat(  
            record.map { it -> it.join(",") }
        )
        .collectFile(name: 'metadata.csv', storeDir: 'data', newLine: true, sort: false)
        .set { metadata }

    // Download a genome if not provided from params
    if ( file("data/${params.genome.tokenize('/')[-1]}").exists() &&
        file("data/${params.genome_annotation.tokenize('/')[-1]}").exists() ) {
        log.info "Found local files for genome and annotation file"
        // When both the genome and its annotation exists in local then used them directly
        genome = file("data/${params.genome.tokenize('/')[-1]}")
        gtf = file("data/${params.genome_annotation.tokenize('/')[-1]}")
    } else {
        // Otherwise download to disk
        ch_refs = Channel.fromList([params.genome, params.genome_annotation])
        //
        // PROCESS: Download relevant references like fa, gtf to data directory
        //
        DOWNLOAD_REFERENCE ( ch_refs , file("data"))
        DOWNLOAD_REFERENCE.out.reference
                        .branch { it ->
                            genome: it.toString().contains('fa.gz')
                            gtf: it.toString().contains('gtf.gz')
                        }
                        .set { refs }
        genome = refs.genome
        gtf = refs.gtf
    }
    //
    // PROCESS: Execute initial quality control on fastq data
    //
    FASTQC ( record )
    //
    // PROCESS: Trim on low quality reads
    //
    TRIMGALORE ( record )
    //
    // PROCESS: Build genome index for align later
    HISAT2_BUILD ( genome )
    //
    // PROCESS: Align the trimmed reads to reference fasta (genome)
    //
    // [ sample_name, [read1, read2] ], [ fa_name, hisat2 ]
    HISAT2_ALIGN ( TRIMGALORE.out.reads, HISAT2_BUILD.out.index )
    //
    // PROCESS: Convert the the aligned sam files to bam
    //
    SAMTOOLS_TO_BAM ( HISAT2_ALIGN.out.sam )
    // 
    // PROCESS: Sort these bam files
    //
    SAMTOOLS_SORT ( SAMTOOLS_TO_BAM.out.bam )
    // 
    // PROCESS: Collect all the bam files generated and pass in 1 as one arg, and count
    //
    FEATURE_COUNTS ( SAMTOOLS_SORT.out.bam.collect(), gtf )

    // =============================================================================
    // Collect versions from modules
    ch_versions = ch_versions
        .mix( HISAT2_BUILD.out.versions )
        .mix( HISAT2_ALIGN.out.versions )
        .mix( SAMTOOLS_TO_BAM.out.versions )
        .mix( SAMTOOLS_SORT.out.versions )
        .mix( FEATURE_COUNTS.out.versions )
//       .mix ( FASTQC.out.versions )
//        .mix ( TRIMGALORE.out.versions )

    // Lastly collect all software versions and to YAML
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'dge_analysis_versions.yml',
            sort: false,
            newLine: true
            )

}
