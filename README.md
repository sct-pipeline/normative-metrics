# Normative qMRI metrics

Set of scripts for computing quantitative MRI (qMRI) metrics from diffusion weighted (DWI) and magnetic transfer (MT) spinal cord MRI data.

## Motivation

To compute normative qMRI metrics from large cohort of multi-center multi-vendor healthy subjects for various spinal cord ROI based on PAM50 atlas and analyze their variability per individual vertebral levels and also per individual vendors.

## Data 

Multi-center multi-vendor [data](https://spine-generic.readthedocs.io/en/latest/index.html) of healthy subjects acquired and analyzed within [spine generic protocol](https://github.com/sct-pipeline/spine-generic) project.

## Prerequisites
[Spinal Cord Toolbox (SCT) v4.3](https://github.com/neuropoly/spinalcordtoolbox)

[Spine generic v1.1](https://github.com/sct-pipeline/spine-generic)

## Pipeline
- [Download](https://openneuro.org/datasets/ds001919/versions/1.0.8/download) multi-center multi-subject data from OpenNeuro platform, e.g.:

`aws s3 sync --no-sign-request s3://openneuro.org/ds001919 ds001919-download/`

- Rename downloaded data:
 
`mv ds001919-download ~/spineGeneric-multi-subject`

- Create a folder where results will be generated:

`mkdir ~/spineGeneric-multi-subject_results`

- Analyze data using [process_data.sh](https://github.com/sct-pipeline/spine-generic/blob/master/processing/process_data.sh) script:

`sct_run_batch -jobs -1 -path-data ~/data/spineGeneric-multi-subject/ -path-output ~/spineGeneric-multi-subject_results/ processing/process_data.sh`

- Compute qMRI metrics from various ROI per individual vertebral levels:

`...`
