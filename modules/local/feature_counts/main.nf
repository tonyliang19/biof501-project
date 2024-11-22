process FEATURE_COUNTS {
    tag "Running feature counts"
    container 'quay.io/biocontainers/bioconductor-rsubread:2.4.0--r40h037d062_0'

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )


    input:
    path(bam)
    path(gtf)

    output:
    path("*.rds"),          emit: feature_count
    path("*.log"),          emit: log
    path("versions.yml"),   emit: versions

    script:
    """
    gunzip -f ${gtf}
    featureCounts.R "${bam}" ${gtf.baseName} > feature_counts.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version) | sed -n 's/^R version \\([0-9.]*\\).*/\\1/p')
        Rsubread: \$(Rscript -e "cat(as.character(packageVersion('Rsubread')))")
    END_VERSIONS
    """
}