#!/bin/bash
# =============================================================
# Use this script now to run nextflow with relevant parameters
# And allowing for running with resume
# Author: Tony Liang
# =============================================================

# Configuration options
NXF_SRC_MAIN=main.nf # Main script to run
BATCH=FALSE # By default logs are printed into the screen, if true then to log
RESUME=FALSE # Initially run without resume
# Add unique timestamp to the log file
timestamp=$(date +"%Y%m%d_%H%M%S")
LOG_FILE=${timestamp}-DGE_pipeline.log
# Input options for the pipeline
SAMPLESHEET=data/samplesheet.csv
OUTDIR=results # Output directory for publishing files
# These profiles should not contain spaces, and order matters
# docker,test might not be equivalent to test,docker
PROFILE="docker"
# Add a parser from cli arg i.e.
# bash run_nxf.sh -r allows for resuming workflow
# bahs run_nxf.sh -b allows for redireting output of the log to DGE_pipeline.log
while getopts ":rb" option;  # include both 'r' and 'b' options
do 
  case $option in
  r)
    echo -e "Resuming pipeline\n"
    RESUME=TRUE
    ;;
  b)
    echo -e "Running as batch job mode, log will be output to ${LOG_FILE}\n"
    BATCH=TRUE
    ;;
  *)
    echo "Invalid option: -$OPTARG"
    exit 1
    ;;
  esac
done

# Then parse the run cmd
#export NXF_OFFLINE='true'
NEXTFLOW_CMD="nextflow run ${NXF_SRC_MAIN} --outdir ${OUTDIR} --samplesheet ${SAMPLESHEET} -profile ${PROFILE}"

if [ "$RESUME" == "TRUE" ]; then
  NEXTFLOW_CMD+=" -resume"
fi

if [ "$BATCH" == "TRUE" ]; then
  NEXTFLOW_CMD+=" -ansi-log false > ${LOG_FILE} 2>&1"  # Redirect output and errors to the log file
fi

# Then execute the run cmd
eval ${NEXTFLOW_CMD}
# if [ "${RESUME}" == "TRUE" ];
# then
#   # Only when provide -r in cli, then run the pipeline with resume
#   nextflow run ${NXF_SRC_MAIN} \
#   --outdir ${OUTDIR} \
#   -resume
# else
#   # Otherwise run without resume
#   nextflow run ${NXF_SRC_MAIN} \
#   --outdir ${OUTDIR}
# fi
