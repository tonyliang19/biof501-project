process SAMTOOLS {
    debug true
    container 'biocontainers/samtools:1.21--h50ea8bc_0'
    input:
    path(x)
    output:
    path(x)
    path("versions.yml"), emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}