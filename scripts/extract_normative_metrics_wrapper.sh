#!/bin/bash

# This is wrapper for extract_normative_metrics.py to allow parallel run using SCT function sct_run_batch

# Variables $1 and $PATH_RESULTS are passed from sct_run_batch script
SUBJECT=$1
PATH_DATA=$(echo $PATH_RESULTS | sed 's|\(.*\)/.*|\1|')		# delete "/results" from the end of PATH_RESULTS variable

python extract_normative_metrics.py -sub ${SUBJECT} -path-data ${PATH_DATA}