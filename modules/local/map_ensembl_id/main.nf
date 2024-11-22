process MAP_ENSEMBL_ID {
    debug true
    tag "Converting emsembl gene id to gene symbol"
    container "quay.io/biocontainers/bioconductor-org.hs.eg.db:3.18.0--r43hdfd78af_0"

    publishDir (
		path: "${params.outdir}/${task.process.tokenize(':').join('/').toLowerCase()}",
		mode: "${params.publish_dir_mode}",
        // https://nextflow.slack.com/archives/C02T98A23U7/p1648120122138739
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    )

    input:
    path(deseq2_result_path)

    output:
    path("*.csv"),                  emit: mapped_id_path
    path("versions.yml"),           emit: versions

    script:
    """
    map_ensembl_id.R ${deseq2_result_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version) | sed -n 's/^R version \\([0-9.]*\\).*/\\1/p')
        org.Hs.eg.db: \$(Rscript -e "cat(as.character(packageVersion('org.Hs.eg.db')))")
    END_VERSIONS
    """
}