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
        // This gives [ sample_name, [ reads ] ]
        .map { row ->
            [ row.sample_name, [ file(row.fastq1), file(row.fastq2) ] ]
        }
        .set { record }

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
        //ch_refs.map { it -> it.tokenize['/'][-1] }.view()
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

    // Execute initial quality control on fastq data
    //FASTQC ( record )

    // Then run trimming adapters from qced fastq
    // TODO: This step might be failing now?
    // TRIMGALORE ( record )
    HISAT2_BUILD ( genome )

    HISAT2_ALIGN ( record, HISAT2_BUILD.out.index )
    
//     // TRINITY ( TRIMGALORE.out.reads )
//     // Collect versions from modules
    ch_versions = ch_versions
        .mix( HISAT2_BUILD.out.versions )
        .mix( HISAT2_ALIGN.out.versions)
//       .mix ( FASTQC.out.versions )
//        .mix ( TRIMGALORE.out.versions )

    // Lastly collect all software versions and to YAML
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'dge_analysis_versions.yml',
            sort: true,
            newLine: true
            )

}
