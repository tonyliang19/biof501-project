process TRINITY {
    debug false
    // Use docker container for now
    container 'biocontainers/trinity:2.11.0--h5ef6573_1'
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
    //tuple val(sample_name), path("*.html"), emit: html
    //tuple val(sample_name), path("*.zip") , emit: zip
    path("trinity_out_dir/Trinity.fasta"), emit: fasta, optional: true
    path("versions.yml"),                   emit: versions

    script:
    // Compute first, then collect the version of binary ran
    """
    Trinity --seqType fq \
        --left ${reads[0]} \
        --right ${reads[1]} \
        --max_memory 2G \
        --CPU 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Trinity: \$( Trinity --version )
    END_VERSIONS
    """
}