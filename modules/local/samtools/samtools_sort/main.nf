process SAMTOOLS_SORT {
    debug true
    tag "Sorting ${bam}"
    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )


    input:
    tuple val(sample_name), path(bam)
    output:
    path("${sample_name}_sorted.bam"),  emit: bam
    path("versions.yml"),               emit: versions

    script:
    """
    samtools sort -o ${sample_name}_sorted.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}