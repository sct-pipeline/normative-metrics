#!/bin/bash

# This is wrapper for extract_normative_metrics.py to allow parallel run using SCT function sct_run_batch
#
# USAGE (using SCT sct_run_batch function):
# 	    sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results -continue-on-error 1 -script scripts/extract_normative_metrics_wrapper.sh
#

# Variables $1 and $PATH_RESULTS are passed from sct_run_batch script
SUBJECT=$1
PATH_DATA=$(echo $PATH_RESULTS | sed 's|\(.*\)/.*|\1|')		# delete "/results" from the end of PATH_RESULTS variable

python /usr/local/normative-metrics/scripts/extract_normative_metrics.py -sub ${SUBJECT} -path-data ${PATH_DATA}