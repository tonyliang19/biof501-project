process HISAT2_BUILD {
    tag "Indexing genome using hisat2"
    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}/",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
	)

    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    path(fa)
    output:
    tuple val(fa.baseName), path("hisat2"), emit: index
    path('versions.yml'),                   emit: versions
    script:
    // Compute first, then collect the version of binary ran

    """
    gunzip ${fa}
    mkdir hisat2
    hisat2-build -p ${task.cpus} \
        ${fa.baseName} \
        hisat2/${fa.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2-build: \$(hisat2-build --version | grep -o 'version [^ ]*' | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS 
    """
}