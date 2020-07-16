#!/usr/bin/env python
#
# Generate figures for individual qMRI metrics
#
# Usage:
#   python generate_figures.py -path-results ~/spineGeneric-multi-subject_results/results/
#
# Optional option:
#   -path-results - directory with *.csv files
#
# Inspired by - https://github.com/sct-pipeline/spine-generic/blob/master/processing/generate_figure.py
# Authors: Julien Cohen-Adad, Jan Valosek

import os
import argparse
import tqdm
import sys
import glob
import csv
import pandas as pd

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import logging

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)

# List subject to remove, associated with contrast
# TODO - update this
SUBJECTS_TO_REMOVE = [
    #{'subject': 'sub-oxfordFmrib04', 'metric': 'csa_t1'},
]

# List of sites to exclude based on the metric
# TODO - update this
SITES_TO_EXCLUDE = {
    'mtr': ['stanford',  # Used different TR.
    ]
    #         'sapienza']  # TODO: check what's going on with this site
    }


# color to assign to each MRI model for the figure
vendor_to_color = {
    'GE': 'black',
    'Philips': 'dodgerblue',
    'Siemens': 'limegreen',
    }

# marker change for individual vendors for scatter plots
# https://matplotlib.org/3.2.2/api/markers_api.html
vendor_to_marker = {
    'GE': 'o',          # circle
    'Philips': 'v',     # triangle down
    'Siemens': 's',     # square
}

# fetch qMRI metrics based on csv file
file_to_metric = {
    'DWI_FA_perlevel.csv': 'dti_fa',
    'DWI_MD_perlevel.csv': 'dti_md',
    'DWI_AD_perlevel.csv': 'dti_ad',
    'DWI_RD_perlevel.csv': 'dti_rd',
    'MTR_perlevel.csv': 'mtr',
    'MTsat_perlevel.csv': 'mtsat',
    }

# fetch metric field
metric_to_field = {
    'dti_fa': 'WA()',
    'dti_md': 'WA()',
    'dti_ad': 'WA()',
    'dti_rd': 'WA()',
    'mtr': 'WA()',
    'mtsat': 'WA()',
    }

# fetch metric field
metric_to_label = {
    'dti_fa': 'Fractional anisotropy',
    'dti_md': 'Mean diffusivity [$mm^2.s^-1$]',
    'dti_ad': 'Axial diffusivity [$mm^2.s^-1$]',
    'dti_rd': 'Radial diffusivity [$mm^2.s^-1$]',
    'mtr': 'Magnetization transfer ratio [%]',
    'mtsat': 'Magnetization transfer saturation [a.u.]',
    }

# scaling factor (for display)
scaling_factor = {
    'dti_fa': 1,
    'dti_md': 1000,
    'dti_ad': 1000,
    'dti_rd': 1000,
    'mtr': 1,
    'mtsat': 1,
    }

# region-of-interest
roi_to_label =  {
    'spinal cord': 'Spinal cord',
    'white matter': 'White matter',
    'gray matter': 'Gray matter',
    'dorsal columns': 'Dorsal columns',
    'lateral funiculi': 'Lateral funiculi/columns',
    'ventral funiculi': 'Ventral funiculi/columns',
}

levels_to_label = {
    '2': 'C2',
    '3': 'C3',
    '4': 'C4',
    '5': 'C5',
}

# FIGURE PARAMETERS
FONTSIZE = 15
TICKSIZE = 10
LABELSIZE = 15

# TODO - modify this function to save metrics for individual ROI and levels
def aggregate_per_site(dict_results, metric):
    """
    Aggregate metrics per site. This function assumes that the file participants.tsv is present in folder ./data/
    :param dict_results:
    :param metric: Metric type
    :return:
    """
    # Build Panda DF of participants based on participants.tsv file
    if os.path.isfile('data/participants.tsv'):
        participants = pd.read_csv(os.path.join('data/participants.tsv'), sep="\t")
    else:
        raise FileNotFoundError("File \"participants.tsv\" was not found.")

    summary_per_vendor(participants)

    # Fetch specific field for the selected metric
    metric_field = metric_to_field[metric]
    # Build a dictionary that aggregates values per site
    results_agg = {}
    # Loop across lines and fill dict of aggregated results
    subjects_removed = []
    for i in tqdm.tqdm(range(len(dict_results)), unit='iter', unit_scale=False, desc="Loop across subjects",
                       ascii=False,
                       ncols=80):
        filename = dict_results[i]['Filename']
        logger.debug('Filename: ' + filename)
        # Fetch metadata for the site
        # dataset_description = read_dataset_description(filename, path_data)
        # cluster values per site
        subject = fetch_subject(filename)
        # check if subject needs to be discarded
        if not remove_subject(subject, metric):
            # Fetch index of row corresponding to subject
            rowIndex = participants[participants['participant_id'] == subject].index
            # Add column "val" with metric value
            participants.loc[rowIndex, 'val'] = dict_results[i][metric_field]
            site = participants['institution_id'][rowIndex].array[0]
            if not site in results_agg.keys():
                # if this is a new site, initialize sub-dict
                results_agg[site] = {}
                results_agg[site][
                    'site'] = site  # need to duplicate in order to be able to sort using vendor AND site with Pandas
                results_agg[site]['vendor'] = participants['manufacturer'][rowIndex].array[0]
                results_agg[site]['model'] = participants['manufacturers_model_name'][rowIndex].array[0]
                # initialize empty sub-dict for metric values with values as list (will be list of
                # metrics for individual subjects within site
                results_agg[site]['val'] = defaultdict(list)
            # get label (ROI name)
            label = dict_results[i]['Label']
            # get val for site (ignore None)
            val = dict_results[i][metric_field]
            # get vertlevel for site
            vertlevel = dict_results[i]['VertLevel']
            if not val == 'None':
                # append data into sub-dict {'vertlevel' 'label': 'metric value'} (key is tuple, values are list)
                results_agg[site]['val'][(vertlevel, label)].append(float(val))
        else:
            subjects_removed.append(subject)
    logger.info("Subjects removed: {}".format(subjects_removed))
    return results_agg

# def compute_stats(results_dict):
#     # loop across individual sites
#     for key in results_dict:
#         # initialize mean value for individual sites
#         results_dict[key]['mean'] = ()
#
#     return results_dict

def summary_per_vendor(participants):
    """
    Some pervendor computation. Only for debug, not used yet.
    :return:
    """
    # compute number of subjects pervendor
    print('Number of subjects per vendor:')
    for vendor in vendor_to_color.keys():
        print('{}: {}'.format(vendor, sum(participants['manufacturer'] == vendor)))

    # compute number of sites pervendor
    print('Number of sites per vendor:')
    for vendor in vendor_to_color.keys():
        print('{}: {}'.format(vendor, pd.unique(participants.loc[participants['manufacturer'] == vendor,
                                                                 'institution_id']).size))


def fetch_subject(filename):
    """
    Get subject from filename
    :param filename:
    :return: subject
    """
    path, file = os.path.split(filename)
    subject = path.split(os.sep)[-2]
    return subject

def remove_subject(subject, metric):
    """
    Check if subject should be removed
    :param subject:
    :param metric:
    :return: Bool
    """
    for subject_to_remove in SUBJECTS_TO_REMOVE:
        if subject_to_remove['subject'] == subject and subject_to_remove['metric'] == metric:
            return True
    return False

def get_parameters():
    parser = argparse.ArgumentParser(
        description="Generate figures. This script needs to be run within the folder with *.csv "
                    "files or this folder can be passed usin -path-results flag.",
        add_help=True,
        prog=os.path.basename(__file__))

    parser.add_argument(
        '-path-results',
        required=False,
        metavar='<dir_path>',
        help="Path to directory with results.")

    args = parser.parse_args()
    return args

def main():

    args = get_parameters()

    # change directory to where are .csv files are located
    if args.path_results is not None:
        if os.path.isdir(args.path_results):
            # Go to results directory defined by user
            os.chdir(args.path_results)
        else:
            raise FileNotFoundError("Directory '{}' was not found.".format(args.path_results))
    else:
        # Stay in current directory (assume it is results directory)
        os.chdir(os.getcwd())

    # fetch perlevel .csv files
    csv_files = glob.glob('*perlevel.csv')

    if not csv_files:
        raise RuntimeError("No *.csv files were found in current directory. You can specify directory with *.csv files "
                           " by -path-results flag.")

    # loop across individual .csv files (i.e., across individual qMRI metrics)
    for csv_file in csv_files:
        print(csv_file)

        # open .csv file and create dict
        logger.info('\nProcessing: ' + csv_file)
        dict_results = []
        with open(csv_file, newline='') as f_csv:
            reader = csv.DictReader(f_csv)
            for row in reader:
                dict_results.append(row)

        # fetch metric name from .csv file
        _, csv_file_small = os.path.split(csv_file)
        metric = file_to_metric[csv_file_small]

        # fetch metric values per site
        results_dict = aggregate_per_site(dict_results, metric)

        #results_dict = compute_stats(results_dict)

        # make it a pandas structure (easier for manipulations)
        df = pd.DataFrame.from_dict(results_dict, orient='index')

        # get individual sites
        site_sorted = df.sort_values(by=['vendor', 'model', 'site']).index.values

        # ------------------------------------------------------------------
        # generate figure - level evolution per ROI
        # ------------------------------------------------------------------
        fig = plt.subplots(figsize=(14, 7))

        # loop across sites
        for site in site_sorted:
            # initialize dict for each loop - {label: mean_metric}
            mean_dict = dict()
            # loop across roi/labels
            for index, label in enumerate(roi_to_label.keys()):
                # create individual subplots
                ax = plt.subplot(2,3,index+1)
                mean_value = list()
                # loop across levels
                for level in levels_to_label.keys():
                    # append mean metric value across individual subjects within site perlevel (C2,C3,C4,C5)
                    mean_value.append(np.mean(df['val'][site][level,label]))
                    # fill dict - {label: mean_metric}, e.g. - {'spinal cord': [0.6, 0.6, 0.5, 0.3], 'white matter': [0.6, 0.7, 0.5, 0.4], ...}
                    mean_dict[label] = mean_value
                # get vendor for current site
                vendor = df.vendor[site]
                # plot mean values perlevel (C2, C3, C4, C5)
                plt.plot([float(key) for key in levels_to_label],
                            mean_dict[label],
                            marker=vendor_to_marker[vendor],    # change marker symbol based on vendor
                            markersize=8,                       # size of marker symbol
                            alpha=0.5,                          # transparency
                            label=site)
                # rename xticks to C2, C3, C4, C5
                plt.xticks([float(key) for key in levels_to_label], levels_to_label.values())
                plt.ylabel(metric_to_label[metric],fontsize=FONTSIZE)
                plt.grid(axis='y', linestyle="--")
                plt.title(roi_to_label[label])
        # place legend next to last subplot
        plt.legend(bbox_to_anchor=(1.5, 1), fontsize=FONTSIZE-5)
        # Move subplots closer to each other
        plt.subplots_adjust(wspace=-0.5)
        plt.tight_layout()
        plt.show()

        print('done')

if __name__ == "__main__":
    main()