process SAMTOOLS_TO_BAM {
    debug true
    tag "Converting ${sam} to bam"
    container 'biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
	)


    input:
    tuple val(sample_name), path(sam)
    output:
    path("${sample_name}.bam"), emit: bam
    path("versions.yml"), emit: versions

    script:
    """
    samtools view -Sb -o ${sample_name}.bam \
        --threads ${task.cpus}  \
        ${sam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}