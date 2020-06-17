#!/bin/bash
#
# Check the presence of input files.
#
# Usage:
#   ./check_input_files.sh <SUBJECT>

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

if [[ $# -ne 1 ]] && [[ $1 != "-h" ]];then

  PATH_LOG="`pwd`"

  # Retrieve input params
  SUBJECT=$1

  # Create BIDS architecture
  PATH_IN="`pwd`/${SUBJECT}"

  # Check if passed subject's directory exists
  if [[ ! -d $PATN_IN ]];then echo "Path $PATN_IN is wrong. Exiting.";exit;fi

  # Verify presence of output files and write log file if error
  FILES_TO_CHECK=(
    "$PATH_IN/anat/${SUBJECT}_acq-T1w_MTS.nii.gz"
    "$PATH_IN/anat/${SUBJECT}_acq-MTon_MTS.nii.gz"
    "$PATH_IN/anat/${SUBJECT}_acq-MToff_MTS.nii.gz"
    "$PATH_IN/anat/${SUBJECT}_T2w.nii.gz"
    "$PATH_IN/anat/${SUBJECT}_T2star.nii.gz"
    "$PATH_IN/anat/${SUBJECT}_T1w.nii.gz"
    "$PATH_IN/dwi/${SUBJECT}_dwi.nii.gz"
  )

  for file in ${FILES_TO_CHECK[@]}; do
    if [ ! -e $file ]; then
      echo "${file} does not exist"
      echo "${file} does not exist" >> $PATH_LOG/_error_check_input_files.log
    fi
  done

else
  echo -n "Usage:\n\t./check_input_files.sh <SUBJECT>"
fi