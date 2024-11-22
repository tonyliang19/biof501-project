process DESEQ2 {
    debug true

    container "quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_2"

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )

    input:
    path(fc_rds_path)
    path(metadata_path)

    output:
    path("*.csv"),                  emit: deseq2_result
    path("*.log"), optional: true,  emit: log
    // path("versions.yml"),           emit: versions

    script:
    """
    deseq2_analysis.R ${fc_rds_path} ${metadata_path} > deseq2_analysis.log
    """
}