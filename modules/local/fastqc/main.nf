process FASTQC {
    debug true
    // Use docker container for now
    container "biocontainers/fastqc:v0.11.9_cv8"
    tag "${row.sample_name}"
    // This process should be running initial quality control on fastq files
    input:
    // This is take in as a map, able to retrieve element from map_name.key_name
    val(row)
    output:
    path("versions.yml"), emit: versions
    stdout

    script:
    // Compute first, then collect the version of binary ran
    """
    echo "hello world"
    echo "This is sample: ${row.sample_name}"
    echo "This is read1: ${row.fastq1}"
    echo "This is read2: ${row.fastq2}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}