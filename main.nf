/*
    This is the main entrance of the workflow for differential
    gene expression analysis
*/

// Including modules
include { FASTQC }                  from "./modules/local/fastqc"
include { TRIMGALORE }              from "./modules/local/trimgalore"
include { TRINITY }                 from "./modules/local/trinity"
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

    // record.view()
    // Execute initial quality control on fastq data
    FASTQC ( record )

    // Then run trimming adapters from qced fastq
    // TODO: This step might be failing now?
    TRIMGALORE ( record )

    TRINITY ( TRIMGALORE.out.reads )
    // Collect versions from modules
    ch_versions
        .mix ( FASTQC.out.versions )
        .mix ( TRIMGALORE.out.versions )

    // Lastly collect all software versions and to YAML
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info", 
            name: 'dge_analysis_versions.yml', 
            sort: true, 
            newLine: true
            )

}