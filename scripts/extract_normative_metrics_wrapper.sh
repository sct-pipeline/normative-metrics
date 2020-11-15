#!/bin/bash

# This is wrapper for extract_normative_metrics.py to allow parallel run using SCT function sct_run_batch
#
# USAGE:
# - parallel mode across multiple subjects (using SCT sct_run_batch function and extract_normative_metrics_wrapper.sh wrapper):
#	    sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results -continue-on-error 1 -script scripts/extract_normative_metrics_wrapper.sh
# - same as above, but with yml file containing labels to process(passed by -script-args option):
#	    sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results -continue-on-error 1 -script-args "-yml-file scripts/labels_to_process.yml" -script scripts/extract_normative_metrics_wrapper.sh
#
# Author: Jan Valosek

# Variables $1 and $PATH_RESULTS are passed from sct_run_batch script
SUBJECT=$1
PATH_DATA=$(echo $PATH_RESULTS | sed 's|\(.*\)/.*|\1|')		# delete "/results" from the end of PATH_RESULTS variable
# Input yml file containing labels to process
INPUT_YML=$3

# Run without yml file
if [[ ${INPUT_YML} == "" ]];then
	  # following command expect that this script is called from normative-metrics repo
	  python scripts/extract_normative_metrics.py -sub ${SUBJECT} -path-data ${PATH_DATA}
# Run with passed yml file
else
    # following command expect that this script is called from normative-metrics repo
    python scripts/extract_normative_metrics.py -sub ${SUBJECT} -path-data ${PATH_DATA} -yml-file ${INPUT_YML}
fi

