process TRIMGALORE {
    debug true
    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'

    tag "${sample_name}"
    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}/${sample_name}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
	)

    input:
    tuple val(sample_name), path(reads)
    output:
    // See here https://github.com/FelixKrueger/TrimGalore/issues/60
    tuple val(sample_name), path("*val*.fq.gz"),    emit: reads
    path("versions.yml"),                           emit: versions

    script:
    """
    trim_galore \\
        --paired \\
        --gzip \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """


}