process FEATURE_COUNTS {
  debug true
  container 'biocontainers/bioconductor-rsubread:2.4.0--r40h037d062_0'
  input:
  path(bam)
  path(gtf)
  output:
  path("*.rds")

  script:
  //featureCounts.R --bam "${bam.join()}" --annot ${gtf}
  """
  gunzip -f ${gtf}
  featureCounts.R "${bam}" ${gtf.baseName}
  """
}