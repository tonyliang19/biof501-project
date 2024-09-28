process FASTQC {
    debug false
    // Use docker container for now
    container "biocontainers/fastqc:v0.11.9_cv8"
    tag "${sample_name}"
    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}/${sample_name}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
	)

    // Input, output, script block
    // This process should be running initial quality control on fastq files
    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    tuple val(sample_name), path(reads)
    output:
    tuple val(sample_name), path("*.html"), emit: html
    tuple val(sample_name), path("*.zip") , emit: zip
    path("versions.yml"),                   emit: versions

    script:
    // Compute first, then collect the version of binary ran
    """
    fastqc -f fastq -q ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}