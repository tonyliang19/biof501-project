process HISAT2_BUILD {
    debug true
    tag "Indexing genome using hisat2"
    container "biocontainers/hisat2:2.0.4--py35_1"

    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    path(fa)
    output:
    path("hisat2") ,        emit: index
    path('versions.yml'),   emit: versions
    script:
    // Compute first, then collect the version of binary ran
    def 

    """
    mkdir hisat2
    hisat2-build -p ${task.cpus} \
        ${fa} \
        hisat2/${fa.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2-build: \$(hisat2-build --version | grep -o 'version [^ ]*' | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS 
    """
}