#!/usr/bin/env python

# -------------------------------------------------------
# Extract qMRI metrics (FA, MD, AD, RD, MTR, MTsat) from
# individual ROI perlevel between C2 and C5
#
# Extracted metrics will be saved as *perlevel.csv files
# in results/perlevel directory
#
# USAGE:
# - parallel mode across multiple subjects (using SCT sct_run_batch function and extract_normative_metrics_wrapper.sh wrapper):
#	    sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results -continue-on-error 1 -script scripts/extract_normative_metrics_wrapper.sh
#
# (you can run the script only on some subjects, using -include flag, see sct_run_batch -h)
#
# - single subject mode with default labels/ROIs:
#       extract_normative_metrics.py -path-data ~/data-multi-subject_results -sub sub-amu01
# - single subject mode with labels/ROIs defined by yml file:
#       extract_normative_metrics.py -path-data ~/data-multi-subject_results -sub sub-amu01 -config labels_to_process.yml
#
# Authors: Jan Valosek, Julien Cohen-Adad
# -------------------------------------------------------

import os
import argparse
import yaml
import re

# dict with qMRI metrics to process - metric: dir_location
metrics_to_process = {
    'FA': 'dwi',
    'MD': 'dwi',
    'AD': 'dwi',
    'RD': 'dwi',
    'MTR': 'anat',
    'MTsat': 'anat',
}

# if no input yml file with labels/ROI to process is passed, then, following labels/ROIs will be used
default_labels_to_process = {
    '50': 'spinal cord',
    '51': 'white matter',
    '52': 'gray matter',
    '53': 'dorsal columns',
    '54': 'lateral funiculi',
    '55': 'ventral funiculi',
}

def get_parameters():
    parser = argparse.ArgumentParser(
        description="Extract qMRI metrics (FA, MD, AD, RD, MTR, MTsat) from individual ROI perlevel between"
                    "C2 and C5 vertebral levels.",
        add_help=True,
        prog=os.path.basename(__file__))

    parser.add_argument(
        '-sub',
        required=True,
        metavar='<sub_name>',
        help="Subject name.")
    parser.add_argument(
        '-path-data',
        required=True,
        metavar='<data_path>',
        help="Path to directory with processed data."
    )
    parser.add_argument(
        '-config',
        required=False,
        metavar='<fname>',
        help="Config yaml file listing ROIs/labels to process.")
    #TODO - add example of yaml file

    args = parser.parse_args()
    return args

def main():

    # fetch input arguments
    args = get_parameters()

    # create dict with labels/ROIs to process based on input yml file
    if args.config is not None:
        # check if input yml config file exists
        if os.path.isfile(args.config):
            fname_yml = args.config
        else:
            raise FileNotFoundError("Input yml file '{}' was not found.".format(args.config))

        # fetch input yml file as dict
        with open(fname_yml, 'r') as stream:
            try:
                dict_exclude_subj = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
    else:
        # if no input yml file with labels/ROIs to process is passed, then, use the default labels/ROIs
        print("No input yml file with labels/ROIs to processed was passed. Continuing with default labels/ROI:")
        labels_to_process = default_labels_to_process
        print("(ID, name)")
        for key, value in labels_to_process.items():
            print('{}: {}'.format(key, value))

    if args.sub is not None and args.path_data is not None:
        # put together -path_data and -sub and create subject_path; check if this dir exists
        subject_path = os.path.join(args.path_data, 'data_processed', args.sub)
        if not os.path.isdir(subject_path):
            raise FileNotFoundError("Path to subject dir {} is invalid.".format(subject_path))
        # create dir perlevel in results folder where *perlevel.csv files will be saved
        results_path = os.path.join(args.path_data, 'results', 'perlevel')
        if not os.path.isdir(results_path):
            os.mkdir(results_path)

    print("Processing subject: {}".format(subject_path))

    # string containing labels/ROIs IDs to process (required in this format by sct_extract_metric fucntion)
    # create str from dict keys
    labels_to_process_str = ' '.join([str(label) for label in labels_to_process.keys()])
    # replace spaces by commas (,)
    labels_to_process_str = re.sub(" ", ",", labels_to_process_str.strip())


    # Extract qMRI metrics perlevel
    for metric, folder in metrics_to_process.items():       # loop across metrics
        # dwi metrics (FA, MD, AD, RD)
        if folder == 'dwi':
            # nifti metric file
            metric_fname = os.path.join(subject_path, folder, 'dti_' + metric + 'nii.gz')
            # path to atlas
            atlas_path = os.path.join(subject_path, folder, 'label', 'atlas')
            # name of output csv file where perlevel values will be saved
            csv_fname = os.path.join(results_path, 'DWI_' + metric + '_perlevel.csv')
            # create command
            command = 'sct_extract_metric -i ' + metric_fname + ' -f ' + atlas_path + ' -l ' + labels_to_process_str \
                      + ' -vert 2:5 ' + ' -perlevel 1 ' + ' -o ' + csv_fname + ' -append 1'

        # mt metrics (MTR, MTsat)
        if folder == 'anat':
            # nifti metric file
            metric_fname = os.path.join(subject_path, folder, metric.lower() + 'nii.gz')
            # path to atlas
            atlas_path = os.path.join(subject_path, folder, 'label_axT1w', 'atlas')
            # path to vertfile
            vertfile_fname = os.path.join(subject_path, folder, 'label_axT1w', 'template', 'PAM50_levels.nii.gz')
            # name of output csv file where perlevel values will be saved
            csv_fname = os.path.join(results_path, metric + '_perlevel.csv')
            # create command
            command = 'sct_extract_metric -i ' + metric_fname + ' -f ' + atlas_path + ' -l ' + labels_to_process_str \
                      + ' -vert 2:5 ' + ' -vertfile ' + vertfile_fname + ' -perlevel 1 ' + ' -o ' + csv_fname \
                      + ' -append 1'

        # run shell command
        if os.path.isfile(metric_fname):
            os.system(command)
        else:
            raise FileNotFoundError("Metric file {} does not exist.".format(metric_fname))


if __name__ == "__main__":
    main()
