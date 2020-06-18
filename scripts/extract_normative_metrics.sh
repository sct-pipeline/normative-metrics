#!/bin/bash

# -------------------------------------------------------
# Extract qMRI metrics (FA, MD, AD, RD, MTR, MTsat) from
# individual ROI perlevel between C2 and C5
#
# Usage:
#     ./extract_normative_metrics.sh <SUBJECT>
# -------------------------------------------------------

# Immediately exit if error
set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Retrieve input params
SUBJECT=$1

# Go to subject results folder
cd ${PATH_RESULTS}/data/${SUBJECT}/

ROI_TO_EXTRACT="50,51,52,53,54,55"
# 50 - spinal cord
# 51 - white matter
# 52 - gray matter
# 53 - dorsal columns
# 54 - lateral columns/funiculi
# 55 - ventral columns/funiculi

# -------------------------------------------------------
# DWI
# -------------------------------------------------------

# Go to dwi folder where DTI metrics (FA, MD, AD and RD) are located
cd dwi

# Compute FA, MD, AD and RD perlevel between C2 and C5

DWI_METRICS_TO_PROCESS=(
"FA"
"MD"
"AD"
"RD"
)

for metric in ${DWI_METRICS_TO_PROCESS[@]};do
  sct_extract_metric -i dti_${metric}.nii.gz -f label/atlas -l ${ROI_TO_EXTRACT} -vert 2:5 -perlevel 1 -o ${PATH_RESULTS}/DWI_${metric}_perlevel.csv -append 1
done

# -------------------------------------------------------
# Magnetic transfer
# -------------------------------------------------------
# Go to anat folder MTR and MTsat metrics are located
cd ../anat

# Compute MTR and MTsat perlevel between C2 and C5
MT_METRICS_TO_PROCESS=(
"mtr"
"mtsat"
)

for metric in ${MT_METRICS_TO_PROCESS[@]};do
  sct_extract_metric -i ${metric}.nii.gz -f label_axT1w/atlas -l ${ROI_TO_EXTRACT} -vert 2:5 -vertfile label_axT1w/template/PAM50_levels.nii.gz -perlevel 1 -o ${PATH_RESULTS}/$(echo $metric | tr '[:lower:]' '[:upper:]')_perlevel.csv -append 1
done
