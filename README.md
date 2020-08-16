# Normative qMRI metrics

Set of scripts for computing quantitative MRI (qMRI) metrics from diffusion weighted (DWI) and magnetic transfer (MT) spinal cord MRI data.

## Motivation

To compute normative qMRI metrics from large cohort of multi-center multi-vendor healthy subjects for various spinal cord ROI based on PAM50 atlas and analyze their variability per individual vertebral levels and also per individual vendors.

## Data 

Multi-center multi-vendor [data](https://spine-generic.readthedocs.io/en/latest/index.html) of healthy subjects [acquired](https://osf.io/tt4z9/) and [analyzed](https://spine-generic.readthedocs.io/en/latest/documentation.html#getting-started) within _spine generic protocol_ project.

## Prerequisites
[Spinal Cord Toolbox (SCT) v4.3](https://github.com/neuropoly/spinalcordtoolbox)

[Spine generic v2.1](https://github.com/sct-pipeline/spine-generic)

## Pipeline
- [Download](https://github.com/spine-generic/data-multi-subject#download-zip-package-recommended) multi-center multi-subject data from GitHub webpage using ``curl`` command:

`curl -o spinegeneric.zip -L https://github.com/spine-generic/data-multi-subject/archive/master.zip`

- Unzip downloaded data:

`unzip spinegeneric.zip`

- Create a folder where results will be generated:

`mkdir ~/data-multi-subject-master_results`

- Analyze multi-subject dataset in parallel mode using [process_data.sh](https://github.com/spine-generic/spine-generic/blob/master/process_data.sh) script:

`sct_run_batch -jobs -1 -path-data ~/data-multi-subject-master -path-output ~/data-multi-subject-master_results/ -continue-on-error 1 -script <PATH_TO_SPINE-GENERIC>/spine-generic/process_data.sh`

- Compute qMRI metrics from various ROI per individual vertebral levels using `extract_normative_metrics.sh` script and save them as *csv files:

`sct_run_batch -jobs -1 -path-data ~/spineGeneric-multi-subject/ -path-output ~/spineGeneric-multi-subject_results/ -include <SUBJECT> -script scripts/extract_normative_metrics.sh`

- Generate figures

`python generate_figures.py -path-results ~/data-multi-subject-master_results/results/`
