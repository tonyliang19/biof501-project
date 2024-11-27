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
include { DESEQ2 }                  from "./modules/local/deseq2"
include { MAP_ENSEMBL_ID }          from "./modules/local/map_ensembl_id"
include { ENHANCED_VOLCANO }        from "./modules/local/enhanced_volcano"
include { softwareVersionsToYAML; 
        print_start_params; 
        print_finish_summary }      from "./modules/nf-core/main.nf"


workflow {
    print_start_params()
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
        // This gives [ sample_name, condition, rep, [ read1, read2 ] ]
        .map { row ->
            [ row.sample_name , row.condition, row.rep, [ file(row.fastq1), file(row.fastq2) ] ]
        }
        .multiMap { sample_name, condition, rep, reads ->
            record: [ sample_name , reads ] 
            meta: [ sample_name, condition, rep ] 
        }
        .set { ch_input }
    
    // Then extract the metadata to save it as csv for later use
    // Need extra row as header of columns
    header_row = Channel.fromList(['sample_name,condition,rep'])
    // Then concat it with the relevant metadata
    header_row
        .concat(  
            ch_input.meta
                    // For the metadata info store it for later use at downstream analysis
                    .map { it -> it.join(",") }
        )
        .collectFile(name: 'metadata.csv', storeDir: 'data', newLine: true, sort: false)
        .set { metadata }

    metadata_ch = metadata

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
    FASTQC ( ch_input.record )
    //
    // PROCESS: Trim on low quality reads
    //
    TRIMGALORE ( ch_input.record )
    //
    // PROCESS: Build genome index for align later
    HISAT2_BUILD ( genome )
    //
    // PROCESS: Align the trimmed reads to reference fasta (genome)
    //
    // [ sample_name, [read1, read2] ], [ hisat2 ], fa_name
    HISAT2_ALIGN ( TRIMGALORE.out.reads, HISAT2_BUILD.out.index, genome.baseName )
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
    //
    // PROCESS: Perform downstream analysis using deseq2
    DESEQ2 ( FEATURE_COUNTS.out.feature_count, metadata )
    // 
    // PROCESS: Map the emsembl ids of genes to its symbol
    MAP_ENSEMBL_ID ( DESEQ2.out.deseq2_result )
    //
    // PROCESS: Make volcano plot of log fold changes of genes
    //
    ENHANCED_VOLCANO ( MAP_ENSEMBL_ID.out.mapped_id_path )




    // =============================================================================
    // Collect versions from modules
    ch_versions = ch_versions
        .mix ( FASTQC.out.versions )
        .mix ( TRIMGALORE.out.versions )
        .mix ( HISAT2_BUILD.out.versions )
        .mix ( HISAT2_ALIGN.out.versions )
        .mix ( SAMTOOLS_TO_BAM.out.versions )
        .mix ( SAMTOOLS_SORT.out.versions )
        .mix ( FEATURE_COUNTS.out.versions )
        .mix ( DESEQ2.out.versions )
        .mix ( MAP_ENSEMBL_ID.out.versions )
        .mix ( ENHANCED_VOLCANO.out.versions )

    // Lastly collect all software versions and to YAML
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'dge_analysis_versions.yml',
            sort: false,
            newLine: true
            )
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    println "Successfully completed"
    /*
    // This bit cannot be run interactively????, only try when sending as pipeline 
    jsonStr = JsonOutput.toJson(params)
    file("${params.outdir}/params.json").text = JsonOutput.prettyPrint(jsonStr)
    */
    print_finish_summary()
}

workflow.onError = {
    println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}

// Use this fun to log the input parameters and other metadata information
// of the pipeline upon starting



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/  


