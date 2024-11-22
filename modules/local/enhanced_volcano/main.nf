process ENHANCED_VOLCANO {
    debug true

    container "quay.io/biocontainers/bioconductor-enhancedvolcano:1.20.0--r43hdfd78af_0"

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )

    input:
    path(mapped_id_path)

    output:
    path("*.png"),                  emit: volcano_plot
    path("*.log"), optional: true,  emit: log
    // path("versions.yml"),           emit: versions

    script:
    """
    plot_volcano.R ${mapped_id_path} > ${task.process.tokenize(':')[-1]}.log
    """
}