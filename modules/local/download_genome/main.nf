process DOWNLOAD_GENOME {
    // Given this downloads large file from internet, so allows some retry
    errorStrategy 'retry', maxRetries: 2
    tag "Downloading genome"
    
    publishDir (
		path: "${download_dir}",
		mode: "${params.publish_dir_mode}"
	)

    input:
    val(link)
    path(download_dir)

    output:
    path("hg38.fa.gz"), emit: genome
    
    script:
    // Check https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips
    """
    wget ${link} -O hg38.fa.gz
    """
}