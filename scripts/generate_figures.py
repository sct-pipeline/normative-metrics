#!/usr/bin/env python
#
# Generate figures for individual qMRI metrics
#
# USAGE (run this script in the directory with *csv files, i.e. "results" directory):
#   python generate_figures.py
# OR specify directory with csv files using -path-results flag:
#   python generate_figures.py -path-results ~/spineGeneric-multi-subject_results/results/
#
# Optional option:
#   -path-results   - directory with *.csv files
#   -config         - input yml config file with subjects to remove (e.g., due to back data quality)
#
# Inspired by - https://github.com/sct-pipeline/spine-generic/blob/master/processing/generate_figure.py
# Authors: Jan Valosek, Julien Cohen-Adad, Alexandru Foias

import os
import argparse
import tqdm
import sys
import glob
import csv
import pandas as pd
import yaml

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import logging

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)


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
    'mtsat': 'Magnetization transfer saturation [%]',
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
def aggregate_per_site(dict_results, metric, dict_exclude_subj):
    """
    Aggregate metrics per site. This function assumes that the file participants.tsv is present in -path-results folder
    :param dict_results:
    :param metric: Metric type
    :return: results_agg: nested dict with metric values per site
    """

    participants = load_participants()

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
        if not remove_subject(subject, metric, dict_exclude_subj):
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
            # because we iterate across individual lines, same subject can be included more than one time in
            # subjects_removed list -> use conversion to set in next command
            subjects_removed.append(subject)
    logger.info("Subjects removed: {}".format(set(subjects_removed)))

    return results_agg

# TODO - impelement remove_subject feature into this function
def aggregate_age_and_sex_per_vendor():
    """
    Aggregate age and sex per individual vendors
    :return: df_age: pandas df with age statistics (min, max, mean, std) per vendor
    :return: df_sex: pandas df with sex (M, F) per vendor
    """

    participants = load_participants()

    # Replace '-' by NaN (some subjects do not have listed age and sex)
    participants = participants.replace('-', np.nan)
    # Convert age from str to float to fit with NaN
    participants['age'] = participants['age'].astype(float)
    # Aggregate age per vendors
    df_age = participants.groupby('manufacturer').age.agg(['min', 'max', 'mean', np.std])
    df_age.columns = ['min age', 'max age', 'mean age', 'std age']

    dict_sex = defaultdict(tuple)
    # Loop across subjects grouped by vendors (i.e, 3 groups - GE, Philips, Siemens)
    for vendor, value in participants.groupby('manufacturer'):
        # Insert number of males and females as a tuple into dict pervendors
        dict_sex[vendor] = (sum(value['sex'].values == 'M'), sum(value['sex'].values == 'F'))

    # convert dict to pandas df
    df_sex = pd.DataFrame.from_dict(dict_sex, orient='index', columns=['M', 'F'])

    # Another option without df.groupby with exclusion of subjects
    # # Iterate across lines (individual subjects)
    # dict_age = defaultdict(list)
    # for _, sub in participants.iterrows():
    #     subject = sub['participant_id']
    #     # loop across vendors
    #     for vendor in vendor_to_color.keys():
    #         if participants.loc[participants['participant_id'] == subject]['manufacturer'].values == vendor:
    #             dict_age[vendor].append(participants.loc[participants['participant_id'] == subject]['age'].values)

    return df_age, df_sex

def summary_per_vendor(df):
    """
    Compute mean values pervendor for individual ROI/labels perlevel and number of used (so, after exclusion) subjects
    and sites per vendor.
    :param df: pandas df with metric values per site
    :return: df_vendor: pandas df with mean metric values per vendor
    :return: df_summary_vendor: pandas df with number of subjects and sites per vendor
    """

    # compute mean per vendor (GE, Philips, Siemens) for individual ROI/labels perlevel
    # initialize dict
    dict_vendor = {}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        # if this is a new vendor, initialize sub-dict
        dict_vendor[vendor] = {}
        # loop across levels and ROI (e.g., '5', 'spinal cord'; '4', 'spinal cord'; ...)
        for label in list(df['mean']['amu'].keys()):
            # if this is a new label, initialize sub-dict
            dict_vendor[vendor][label] = {}
            mean_values = list()
            # loop across sites for given vendor
            for site in df[df['vendor'] == vendor]['site']:
                # collect mean values from all sites for given vendor to one list
                mean_values.append(float(df['mean'][site][label]))
            # computer mean from all sites' mean values for given vendor and insert it into dict
            dict_vendor[vendor][label] = np.mean(mean_values)

    df_vendor = pd.DataFrame.from_dict(dict_vendor, orient='index')

    # compute number of subjects pervendor
    dict_sub_per_vendor = dict()        # e.g.: {'GE': 7, 'Philips': 2, 'Siemens': 9}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        num_of_sub = 0
        # loop across sites for given vendor
        for site in df[df['vendor'] == vendor]['site']:
            # get number of used subjects for given site and add it to num_of_sub variable
            num_of_sub = num_of_sub + len(df[df['vendor'] == vendor]['val'][site]['5', 'spinal cord'])
        dict_sub_per_vendor[vendor] = num_of_sub

    # compute number of sites pervendor
    dict_sites_per_vendor = dict()      # e.g.: {'GE': 2, 'Philips': 1, 'Siemens': 4}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        dict_sites_per_vendor[vendor] = sum(df['vendor'] == vendor)

    # convert and concatenate dict_sub_per_vendor and dict_sites_per_vendor to pandas DF
    df_summary_vendor = pd.concat([pd.DataFrame.from_dict(dict_sub_per_vendor, orient='index'),
                                   pd.DataFrame.from_dict(dict_sites_per_vendor, orient='index')], axis=1, sort=False)
    df_summary_vendor.columns = ['num. of sub. per vendor', 'num. of sites per vendor']

    return df_vendor, df_summary_vendor

def generate_level_evolution_persite(df, metric, path_output):
    """
    Generate figure for each metric - level evolution (C2, C3, C4, C5) per ROI for individual sites
    :param df: pandas df with input data
    :param metric: currently processed qMRI metric (e.g. FA, MD, MTsat,...)
    :param path_output: path where figures will be generated (manually entered -path-results path or current dir)
    :return:
    """
    # get individual sites
    site_sorted = df.sort_values(by=['vendor', 'model', 'site']).index.values

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
            ax = plt.subplot(2, 3, index + 1)
            mean_value = list()
            # loop across levels
            for level in levels_to_label.keys():
                # append mean metric value across individual subjects within site perlevel (C2,C3,C4,C5)
                mean_value.append(df['mean'][site][level, label])
                # fill dict - {label: mean_metric}, e.g. - {'spinal cord': [0.6, 0.6, 0.5, 0.3], 'white matter': [0.6, 0.7, 0.5, 0.4], ...}
                mean_dict[label] = mean_value
            # get vendor for current site
            vendor = df.vendor[site]
            # plot mean values perlevel (C2, C3, C4, C5)
            x = [float(key) for key in levels_to_label]  # individual levels - 2,3,4,5
            y = mean_dict[label]  # mean values per levels
            plt.plot(x,
                     y,
                     marker=vendor_to_marker[vendor],  # change marker symbol based on vendor
                     markersize=8,  # size of marker symbol
                     alpha=0.5,  # transparency
                     label=site)
            # rename xticks to C2, C3, C4, C5
            plt.xticks(x, levels_to_label.values())
            # add ylabel for every subplot
            # plt.ylabel(metric_to_label[metric],fontsize=FONTSIZE)
            # add grid
            plt.grid(axis='y', linestyle="--")
            # add title to individual subpolots (i.e., ROI/label)
            plt.title(roi_to_label[label])

    # place legend next to last subplot
    plt.legend(bbox_to_anchor=(1.1, 1), fontsize=FONTSIZE - 5)
    # Move subplots closer to each other
    plt.subplots_adjust(wspace=-0.5)
    plt.tight_layout()
    # tight layout of whole figure and shift master title up
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    # save figure
    fname_fig = os.path.join(path_output, metric + '_per_sites.png')
    plt.savefig(fname_fig, dpi=200)
    logger.info('Created: ' + fname_fig)

    # plt.show()

def load_participants():
    # Build Panda DF of participants based on participants.tsv file
    if os.path.isfile('participants.tsv'):
        participants = pd.read_csv(os.path.join('participants.tsv'), sep="\t")
    else:
        raise FileNotFoundError("File \"participants.tsv\" was not found in {} folder.".format(os.getcwd()))

    return participants

def fetch_subject(filename):
    """
    Get subject from filename
    :param filename:
    :return: subject
    """
    path, file = os.path.split(filename)
    subject = path.split(os.sep)[-2]
    return subject

def remove_subject(subject, metric, dict_exclude_subj):
    """
    Check if subject should be removed
    :param subject:
    :param metric:
    :param dict_exclude_subj: dictonary with subjects to remove from analysis (due to bad data quality, etc.)
    :return: Bool
    """
    if metric in dict_exclude_subj.keys():
        if subject in dict_exclude_subj[metric]:
            return True
    return False

def get_parameters():
    parser = argparse.ArgumentParser(
        description="Generate figures. This script needs to be run within the folder with *.csv "
                    "files or this folder can be passed using -path-results flag.",
        add_help=True,
        prog=os.path.basename(__file__))

    parser.add_argument(
        '-path-results',
        required=False,
        metavar='<dir_path>',
        help="Path to directory with results (*.csv files).")
    parser.add_argument(
        '-config',
        required=False,
        metavar='<fname>',
        help="Config yaml file listing images that you want to remove from processing.")

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
        # initialize empty dict if no config yml file is passed
        dict_exclude_subj = dict()

    # change directory to where are .csv files are located (and where figures will be generated)
    if args.path_results is not None:
        if os.path.isdir(args.path_results):
            path_output = args.path_results
            # Go to results directory defined by user
            os.chdir(path_output)
        else:
            raise FileNotFoundError("Directory '{}' was not found.".format(args.path_results))
    else:
        # Stay in current directory (assuming it is results directory)
        path_output = os.getcwd()
        print("-path-results flag has not been set. Assuming current directory as a directory with *csv files.")
        os.chdir(path_output)

    # fetch perlevel .csv files
    csv_files = glob.glob('*perlevel.csv')

    if not csv_files:
        raise RuntimeError("No *.csv files were found in the current directory. You can specify directory with *.csv "
                           "files by -path-results flag.")

    # loop across individual .csv files (i.e., across individual qMRI metrics)
    for csv_file in csv_files:

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
        results_dict = aggregate_per_site(dict_results, metric, dict_exclude_subj)
        # make it a pandas structure (easier for manipulations)
        df = pd.DataFrame.from_dict(results_dict, orient='index')

        # ------------------------------------------------------------------
        # compute per vendor summary
        # ------------------------------------------------------------------
        df_vendor, df_summary_vendor = summary_per_vendor(df)
        # and save it as .csv file
        fname_csv_per_vendor = os.path.join(os.getcwd(), metric) + '_per_vendors.csv'
        df_vendor.to_csv(fname_csv_per_vendor)
        logger.info('Created: ' + fname_csv_per_vendor)

        print(df_vendor.head())

        # ------------------------------------------------------------------
        # compute age and sex per vendor
        # ------------------------------------------------------------------
        df_age, df_sex = aggregate_age_and_sex_per_vendor()

        # Concatenate number of sites and subjects with sex and age pervendor
        df_summary_vendor = pd.concat([df_summary_vendor, df_age, df_sex], sort=False, axis=1)

        print(df_summary_vendor)

        # ------------------------------------------------------------------
        # generate figure - level evolution per ROI for individual sites
        # ------------------------------------------------------------------

        generate_level_evolution_persite(df, metric, path_output)

        print('done')

if __name__ == "__main__":
    main()