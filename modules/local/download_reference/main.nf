process DOWNLOAD_REFERENCE {
    // Given this downloads large file from internet, so allows some retry
    errorStrategy 'retry', maxRetries: 2
    // Use split here rather than tokenize since its string
    tag "Downloading reference ${link.split('/')[-1]}"
    
    publishDir (
		path: "${download_dir}",
		mode: "${params.publish_dir_mode}"
	)

    input:
    val(link)
    path(download_dir)

    output:
    path("*.gz"), emit: reference
    
    script:
    // Check https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips
    // Check https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/
    def filename = link.split('/')[-1]
    """
    wget ${link} -O ${filename}
    """
}