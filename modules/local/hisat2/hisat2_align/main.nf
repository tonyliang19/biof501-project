process HISAT2_ALIGN {
    debug true
    tag "Aligning genome using hisat2"
    container "biocontainers/hisat2:2.0.4--py35_1"

    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    tuple val(sample_name), path(reads)
    output:
    path('versions.yml'), emit: versions
    script:
    // Compute first, then collect the version of binary ran
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS 
    """
}