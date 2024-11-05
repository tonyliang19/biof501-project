process FEATURE_COUNTS {
  input:
  path(bam)
  path(gtf)
  output:
  path(x)

  script:
  """
  featureCounts.R --bam ${bam} --annot ${gtf}
  """
}