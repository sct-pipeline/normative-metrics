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
    'dti_md': 'Mean diffusivity [$mm^2.s^{-1}$]',
    'dti_ad': 'Axial diffusivity [$mm^2.s^{-1}$]',
    'dti_rd': 'Radial diffusivity [$mm^2.s^{-1}$]',
    'mtr': 'Magnetization transfer ratio [%]',
    'mtsat': 'Magnetization transfer saturation [a.u.]',
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
    :return: results_agg: nested dict with metric values per site
    """
    # Build Panda DF of participants based on participants.tsv file
    if os.path.isfile('data/participants.tsv'):
        participants = pd.read_csv(os.path.join('data/participants.tsv'), sep="\t")
    else:
        raise FileNotFoundError("File \"participants.tsv\" was not found.")

    # Fetch specific field for the selected metric
    metric_field = metric_to_field[metric]
    # Build a dictionary that aggregates values per site
    results_agg = {}
    # Loop across lines and fill dict of aggregated results
    subjects_removed = []
    for i in tqdm.tqdm(range(len(dict_results)), unit='iter', unit_scale=False, desc="Loop across lines.",
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
                # initialize empty sub-dict for metric values with values as list
                # metrics for individual subjects within site
                results_agg[site]['val'] = defaultdict(list)
                # initialize empty sub-dict for metric mean values
                # metrics mean for within site
                results_agg[site]['mean'] = defaultdict(int)
            # get label (ROI name)
            label = dict_results[i]['Label']
            # get val for site (ignore None)
            val = dict_results[i][metric_field]
            # get vertlevel for site
            vertlevel = dict_results[i]['VertLevel']
            if not val == 'None':
                # append data into sub-dict {'vertlevel' 'label': 'metric value'} (key is tuple, values are list)
                results_agg[site]['val'][(vertlevel, label)].append(float(val))
            # compute mean perlevel per ROI/label
            results_agg[site]['mean'][(vertlevel, label)] = np.mean(results_agg[site]['val'][(vertlevel, label)])

        else:
            subjects_removed.append(subject)
    logger.info("Subjects removed: {}".format(subjects_removed))

    return results_agg

def summary_per_vendor(df):
    """
    Compute number of used (so, after exclusion) subjects and sites per vendor.
    :return:
    """
    # compute number of subjects pervendor
    print('Number of subjects per vendor:')

    # loop across vendors
    for vendor in vendor_to_color.keys():
        num_of_sub = 0
        # loop across sites for given vendor
        for site in df[df['vendor'] == vendor]['site']:
            # get number of used subjects for given site and add it to num_of_sub variable
            num_of_sub = num_of_sub + len(df[df['vendor'] == vendor]['val'][site]['5', 'spinal cord'])

        print('{}: {}'.format(vendor, num_of_sub))

    # compute number of sites pervendor
    print('Number of sites per vendor:')
    # loop across vendors
    for vendor in vendor_to_color.keys():
        print('{}: {}'.format(vendor, sum(df['vendor'] == vendor)))


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

        # make it a pandas structure (easier for manipulations)
        df = pd.DataFrame.from_dict(results_dict, orient='index')

        # get individual sites
        site_sorted = df.sort_values(by=['vendor', 'model', 'site']).index.values

        # compute per vendor summary
        summary_per_vendor(df)

        # ------------------------------------------------------------------
        # generate figure - level evolution per ROI for individual sites
        # ------------------------------------------------------------------
        fig, _ = plt.subplots(figsize=(14, 7))
        # add master title for whole figure
        fig.suptitle(metric_to_label[metric], fontsize=FONTSIZE)

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
                    mean_value.append(df['mean'][site][level,label])
                    # fill dict - {label: mean_metric}, e.g. - {'spinal cord': [0.6, 0.6, 0.5, 0.3], 'white matter': [0.6, 0.7, 0.5, 0.4], ...}
                    mean_dict[label] = mean_value
                # get vendor for current site
                vendor = df.vendor[site]
                # plot mean values perlevel (C2, C3, C4, C5)
                x = [float(key) for key in levels_to_label]    # individual levels - 2,3,4,5
                y = mean_dict[label]                            # mean values per levels
                plt.plot(x,
                         y,
                         marker=vendor_to_marker[vendor],    # change marker symbol based on vendor
                         markersize=8,                       # size of marker symbol
                         alpha=0.5,                          # transparency
                         label=site)
                # rename xticks to C2, C3, C4, C5
                plt.xticks(x, levels_to_label.values())
                # add ylabel for every subplot
                #plt.ylabel(metric_to_label[metric],fontsize=FONTSIZE)
                # add grid
                plt.grid(axis='y', linestyle="--")
                # add title to individual subpolots (i.e., ROI/label)
                plt.title(roi_to_label[label])

        # place legend next to last subplot
        plt.legend(bbox_to_anchor=(1.1, 1), fontsize=FONTSIZE-5)
        # Move subplots closer to each other
        plt.subplots_adjust(wspace=-0.5)
        plt.tight_layout()
        # tight layout of whole figure and shift master title up
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)
        # save figure
        fname_fig = os.path.join(args.path_results, metric + '.png')
        plt.savefig(fname_fig, dpi=200)
        logger.info('Created: ' + fname_fig)

        #plt.show()

        print('done')

if __name__ == "__main__":
    main()