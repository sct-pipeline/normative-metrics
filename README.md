# Normative qMRI metrics

Set of scripts for computing quantitative MRI (qMRI) metrics from diffusion weighted (DWI) and magnetic transfer (MT) spinal cord MRI data.

## Motivation

To compute normative qMRI metrics from large cohort of multi-center multi-vendor healthy subjects for various spinal cord ROI based on PAM50 atlas and analyze their variability per individual vertebral levels and also per individual vendors.

## Data 

Multi-center multi-vendor [data](https://spine-generic.readthedocs.io/en/latest/index.html) of healthy subjects [acquired](https://osf.io/tt4z9/) and [analyzed](https://spine-generic.readthedocs.io/en/latest/documentation.html#getting-started) within _spine generic protocol_ project.

## Prerequisites
[Spinal Cord Toolbox (SCT) v4.3](https://github.com/neuropoly/spinalcordtoolbox)

[Spine generic v2.1](https://github.com/sct-pipeline/spine-generic)

## Requirements

pandas

pyyaml

matplotlib

numpy

(All these packages are included in SCT virtual env (`venv_sct`), activate it:

```
source ${SCT_DIR}/python/etc/profile.d/conda.sh
conda activate venv_sct
```

## Pipeline
- [Download](https://github.com/spine-generic/data-multi-subject#download) multi-center multi-subject data from GitHub webpage using ``git annex``:

```
git clone https://github.com/spine-generic/data-multi-subject && cd data-multi-subject && git annex init && git annex get
```

- Create a folder where results will be generated:

```
mkdir ~/data-multi-subject_results
```

- Analyze multi-subject dataset in parallel mode using [process_data.sh](https://github.com/spine-generic/spine-generic/blob/master/process_data.sh) script:

```
sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results/ -continue-on-error 1 -script <PATH_TO_SPINE-GENERIC>/process_data.sh
```

- Compute qMRI metrics from various ROI per individual vertebral levels across all subjects using `extract_normative_metrics_wrapper.sh` script :

```
sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results/ -continue-on-error 1 -script scripts/extract_normative_metrics_wrapper.sh
```

(You can run the script only for specific subjects using `-include <SUBJECT>` flag)

(Individual \*perlevel.csv files will be stored in `/results/perlevel` folder)

- Generate figures

```
python generate_figures.py -path-results ~/data-multi-subject-master_results/results/
```
