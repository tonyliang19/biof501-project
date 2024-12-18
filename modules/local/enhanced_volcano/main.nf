process ENHANCED_VOLCANO {
    debug true
    tag "Plotting volcano of log2fc"
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
    path("versions.yml"),           emit: versions

    script:
    """
    plot_volcano.R ${mapped_id_path}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version) | sed -n 's/^R version \\([0-9.]*\\).*/\\1/p')
        EnhancedVolcano: \$(Rscript -e "cat(as.character(packageVersion('EnhancedVolcano')))")
        ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}