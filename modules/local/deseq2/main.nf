process DESEQ2 {
  debug true
  input:
  path(fc_rds_path)
  path(metadata_path)

  output:
  // path("*.csv"),                  emit: deseq2_result
  path("*.log"), optional: true,  emit: log
  // path("versions.yml"),           emit: versions

  script:
  //   deseq2_analysis.R ${fc_rds_path} ${metadata_path} > deseq2_analysis.log
  """
  echo "This is fc_rds_path: ${fc_rds_path}"
  echo "This is metdata path: ${metadata_path}"
  """
}