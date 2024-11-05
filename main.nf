/*
    This is the main entrance of the workflow for differential
    gene expression analysis
*/

// Including modules
include { FASTQC }                  from "./modules/local/fastqc"
include { HISAT2_BUILD }            from "./modules/local/hisat2/hisat2_build"
include { HISAT2_ALIGN }            from "./modules/local/hisat2/hisat2_align"
include { DOWNLOAD_GENOME }         from "./modules/local/download_genome"
//include { TRIMGALORE }              from "./modules/local/trimgalore"
//include { TRINITY }                 from "./modules/local/trinity"
include { softwareVersionsToYAML }  from "./modules/nf-core/main.nf"

def getGenomeName(genome_url) {
    // Parse the url containing genome and get its name
    def name = genome_url.tokenize('/')[-1]
    return name.replaceAll(/\.(fa|fasta|fa\.gz|fasta\.gz)$/, '') // Strips extensions like ".fa", ".fasta", ".fa.gz", or ".fasta.gz"
    
}

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
    // if (params.genome.startsWith('http') || params.genome.startsWith('ftp')) {
    //     def genome_file = file("data/${getGenomeName(params.genome)}.fa.gz")

    //     if (genome_file.exists()) {
    //         println "Found genome, skipping download"
    //         genome = genome_file
    //         return genome
    //     } else {
    //         println "Starting download of genome"
    //         genome = DOWNLOAD_GENOME(params.genome, file("data/"))
    //         return genome
    //     }
    // }
    // record.view()
    // Execute initial quality control on fastq data
    //FASTQC ( record )

    // Then run trimming adapters from qced fastq
    // TODO: This step might be failing now?
    // TRIMGALORE ( record )
    HISAT2_BUILD ( genome )

    HISAT2_ALIGN ( HISAT2_BUILD.out.index )
    
    // TRINITY ( TRIMGALORE.out.reads )
    // Collect versions from modules
    ch_versions
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
