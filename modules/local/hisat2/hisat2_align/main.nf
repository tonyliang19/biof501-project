process HISAT2_ALIGN {
    tag "Aligning genome using hisat2 at sample ${sample_name}"
    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}/",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
	)

    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    tuple val(sample_name), path(reads)
    tuple val(fa_name), path(index)
    output:
    path("${sample_name}.sam"), emit: sam
    path("*.log"),              emit: log
    path('versions.yml'),       emit: versions
    script:
    // Compute first, then collect the version of binary ran
    """
    hisat2 -x ${index}/${fa_name} \
        -1 ${reads[0]} -2 ${reads[1]} \
        -S ${sample_name}.sam \
        --summary-file hisat2-${sample_name}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS 
    """
}