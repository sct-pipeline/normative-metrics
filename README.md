# Normative qMRI metrics

Set of scripts for computing of quantitative MRI (qMRI) metrics from diffusion-weighted (DWI/dMRI) and magnetization transfer (MT) spinal cord data.

## Motivation

To compute the normative qMRI metrics from a large cohort of healthy subjects scanned across multiple sites for various spinal cord ROI based on PAM50 atlas and analyze their variability per individual vertebral levels (C2-C5) and also per individual vendors (Siemens, Philips, GE).

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
cd ~
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

(Analysis will automatically use manually-corrected labels and segmentation located in `derivatives/labels/`)

- Compute qMRI metrics from various ROI per individual vertebral levels across all subjects using `extract_normative_metrics_wrapper.sh` script :

```
cd ~
git clone https://github.com/sct-pipeline/normative-metrics.git
cd normative-metrics
sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results/ -continue-on-error 1 -script scripts/extract_normative_metrics_wrapper.sh
```

(You can run the script only for specific subjects using `-include <SUBJECT>` flag)

(Individual \*perlevel.csv files will be stored in `/results/perlevel` folder)

- Generate figures and compute statistics

```
python generate_figures.py -path-results ~/data-multi-subject-master_results/results/perlevel -config ~/data-multi-subject_results/results/exclude.yml -participants-file  ~/data-multi-subject_results/results/participants.tsv
```
where:

`-path-results` - path to directory with *csv perlevel files (created in the previous step)

`-participants-file` - path to .tsv file with participants characteristics (age, sex, ...)

`-config` - path to .yml file containing subjects to exclude
