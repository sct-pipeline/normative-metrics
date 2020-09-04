#!/usr/bin/env python

# -------------------------------------------------------
# Extract qMRI metrics (FA, MD, AD, RD, MTR, MTsat) from
# individual ROI perlevel between C2 and C5
#
# USAGE:
# - parallel mode across multiple subjects (using SCT function sct_run_batch and extract_normative_metrics.sh wrapper):
#	    sct_run_batch -jobs -1 -path-data ~/data-multi-subject/ -path-output ~/data-multi-subject_results -continue-on-error 1 -script scripts/extract_normative_metrics.sh
#
# (you can run the script only on some subjects, using -include flag)
#
# - single subject mode:
#       extract_normative_metrics.py -path-data ~/data-multi-subject_results -sub sub-amu01
# -------------------------------------------------------

import os
import argparse
import yaml

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

    args = get_parameters()

    # create dict with subjects to exclude if input yml config file is passed
    if args.config is not None:
        # check if input yml config file exists
        if os.path.isfile(args.config):
            fname_yml = args.config
        else:
            raise FileNotFoundError("Input yml config file '{}' was not found.".format(args.config))

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
        for key, value in labels_to_process.items():
            print('{}: {}'.format(key, value))

    # put together -path_data and -sub and create subject_path; check if this dir exists
    if args.sub is not None and args.path_data is not None:
        subject_path = os.path.join(args.path_data, 'data_processed', args.sub)
        if not os.path.isdir(subject_path):
            raise FileNotFoundError("Path to subject dir {} is invalid.".format(subject_path))

    print("Subject: {}".format(subject_path))

    for key in labels_to_process.keys():
        print(key)

if __name__ == "__main__":
    main()
