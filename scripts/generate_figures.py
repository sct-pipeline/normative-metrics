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
# Authors: Jan Valosek, Julien Cohen-Adad

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

FNAME_LOG = 'log_stats.txt'

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

# method used for qMRI metric extraction (e.g., 'WA()' or 'MAP()'), see help of sct_extract_metric function
extraction_method = 'MAP()'

# fetch metric field
metric_to_field = {
    'dti_fa': extraction_method,
    'dti_md': extraction_method,
    'dti_ad': extraction_method,
    'dti_rd': extraction_method,
    'mtr': extraction_method,
    'mtsat': extraction_method,
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
roi_to_label = {
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

# scaling factor (for display)
scaling_factor = {
    'dti_fa': 1,
    'dti_md': 1000,     # scale to mm^2s^-1
    'dti_ad': 1000,     # scale to mm^2s^-1
    'dti_rd': 1000,     # scale to mm^2s^-1
    'mtr': 1,
    'mtsat': 1,
    }

# FIGURE PARAMETERS
FONTSIZE = 15
TICKSIZE = 10
LABELSIZE = 15

# TODO - modify this function to save metrics for individual ROI and levels
def aggregate_per_site(dict_results, metric, dict_exclude_subj, path_participants):
    """
    Aggregate metrics per site.
    :param dict_results:
    :param metric: Metric type
    :param path_participants: path to participants.tsv file
    :return: results_agg: nested dict with metric values per site
    """

    participants_df = load_participants_file(path_participants)

    # Fetch specific field for the selected metric
    metric_field = metric_to_field[metric]
    # Build a dictionary that aggregates values per site
    results_agg = {}
    # Loop across lines and fill dict of aggregated results
    subjects_removed = []
    for i in tqdm.tqdm(range(len(dict_results)), unit='iter', unit_scale=False, desc="Loop across lines",
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
            rowIndex = participants_df[participants_df['participant_id'] == subject].index
            # Add column "val" with metric value
            participants_df.loc[rowIndex, 'val'] = dict_results[i][metric_field]
            site = participants_df['institution_id'][rowIndex].array[0]
            if not site in results_agg.keys():
                # if this is a new site, initialize sub-dict
                results_agg[site] = {}
                # need to duplicate site in order to be able to sort using vendor AND site with Pandas
                results_agg[site]['site'] = site
                # save vendor (e.g., Siemens)
                results_agg[site]['vendor'] = participants_df['manufacturer'][rowIndex].array[0]
                # save model (e.g., Verio)
                results_agg[site]['model'] = participants_df['manufacturers_model_name'][rowIndex].array[0]
                # initialize empty sub-dict for metric values with values as a list
                # metrics for individual subjects within site
                results_agg[site]['val'] = defaultdict(list)
                # initialize empty sub-dict for metric mean values
                # metrics mean for within site
                results_agg[site]['mean'] = defaultdict(int)
            # get label (ROI name) (e.g., spinal cord, white matter, ...)
            label = dict_results[i]['Label']
            # if label is in roi_to_label dict, process this label, i.e., skip labels which we do not want to process
            if label in roi_to_label.keys():
                # get val for site
                val = dict_results[i][metric_field]
                if not val == 'None':
                    val = float(dict_results[i][metric_field]) * scaling_factor[metric]     # scale metric
                    # get vertlevel for site
                    vertlevel = dict_results[i]['VertLevel']        # e.g., 5
                    # append data into sub-dict  - {'vertlevel' 'label': 'metric values'} (key is tuple, value is list)
                    results_agg[site]['val'][(vertlevel, label)].append(float(val))
                # compute mean perlevel per ROI/label - {'vertlevel' 'label': 'mean value'} (key is tuple, value is float)
                results_agg[site]['mean'][(vertlevel, label)] = np.mean(results_agg[site]['val'][(vertlevel, label)])
            else:
                logger.info('Skipping {}.'.format(label))

        else:
            # because we iterate across individual lines, same subject can be included more than one time in
            # subjects_removed list -> use conversion to set in next command
            subjects_removed.append(subject)
    logger.info("Subjects removed: {}".format(set(subjects_removed)))

    return results_agg

# TODO - impelement remove_subject feature into this function
def aggregate_age_and_sex_per_vendor(path_participants):
    """
    Aggregate age and sex per individual vendors
    :param path_participants: path to participants.tsv file
    :return: df_age: pandas df with age statistics (min, max, mean, std) per vendor
    :return: df_sex: pandas df with sex (M, F) per vendor
    """

    participants_df = load_participants_file(path_participants)

    # Replace '-' by NaN (some subjects do not have listed age and sex)
    participants_df = participants_df.replace('-', np.nan)
    # Convert age from str to float to fit with NaN
    participants_df['age'] = participants_df['age'].astype(float)
    # Aggregate age per vendors
    df_age = participants_df.groupby('manufacturer').age.agg(['min', 'max', 'mean', np.std])
    df_age.columns = ['min age', 'max age', 'mean age', 'std age']

    dict_sex = defaultdict(tuple)
    # Loop across subjects grouped by vendors (i.e, 3 groups - GE, Philips, Siemens)
    for vendor, value in participants_df.groupby('manufacturer'):
        # Insert number of males and females as a tuple into dict pervendors
        dict_sex[vendor] = (sum(value['sex'].values == 'M'), sum(value['sex'].values == 'F'))

    # convert dict to pandas df
    df_sex = pd.DataFrame.from_dict(dict_sex, orient='index', columns=['M', 'F'])

    # Another option without df.groupby with exclusion of subjects
    # # Iterate across lines (individual subjects)
    # dict_age = defaultdict(list)
    # for _, sub in participants_df.iterrows():
    #     subject = sub['participant_id']
    #     # loop across vendors
    #     for vendor in vendor_to_color.keys():
    #         if participants_df.loc[participants_df['participant_id'] == subject]['manufacturer'].values == vendor:
    #             dict_age[vendor].append(participants_df.loc[participants_df['participant_id'] == subject]['age'].values)

    return df_age, df_sex

def summary_per_vendor(df):
    """
    Compute mean and std values pervendor for individual ROI/labels perlevel and number of used (so, after exclusion)
    subjects and sites per vendor.
    :param df: pandas df with metric values per site
    :return: df_vendor_mean: pandas df with mean metric values per vendor
    :return: df_vendor_std: pandas df with std metric values per vendor
    :return: df_summary_vendor: pandas df with number of subjects and sites per vendor
    """

    # compute mean and std per vendor (GE, Philips, Siemens) for individual ROI/labels perlevel
    # initialize dict
    dict_vendor_mean = {}
    dict_vendor_std = {}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        # if this is a new vendor, initialize sub-dict
        dict_vendor_mean[vendor] = {}
        dict_vendor_std[vendor] = {}
        # loop across levels and ROI (e.g., '5', 'spinal cord'; '4', 'spinal cord'; ...)
        for label in list(df['mean']['amu'].keys()):        # <class 'tuple'>: (e.g., '5', 'ventral funiculi')
            # if this is a new label, initialize sub-dict
            dict_vendor_mean[vendor][label] = {}
            dict_vendor_std[vendor][label] = {}
            mean_values = list()
            # loop across sites for given vendor
            for site in df[df['vendor'] == vendor]['site']:
                # collect mean values from all sites for given vendor to one list
                mean_values.append(float(df['mean'][site][label]))
            # compute mean from all sites' mean values for given vendor and insert it into dict
            dict_vendor_mean[vendor][label] = np.mean(mean_values)
            # compute std from all sites' mean values for given vendor and insert it into dict
            dict_vendor_std[vendor][label] = np.std(mean_values)

    # convert dict to pandas df for easy csv save
    df_vendor_mean = pd.DataFrame.from_dict(dict_vendor_mean, orient='index')
    df_vendor_std = pd.DataFrame.from_dict(dict_vendor_std, orient='index')


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

    return df_vendor_mean, df_vendor_std, df_summary_vendor

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

def generate_level_evolution_pervendor(df_vendor, df_summary_vendor, metric, path_output):
    """
    Generate figure for each metric - level evolution (C2, C3, C4, C5) per ROI for individual vendors
    :param df_vendor: pandas df with aggregated data pervendor
    :param df_summary_vendor: pandas df with number of subjects and sites per vendor
    :param metric: currently processed qMRI metric (e.g. FA, MD, MTsat,...)
    :param path_output: path where figures will be generated (manually entered -path-results path or current dir)
    :return:
    """

    fig, _ = plt.subplots(figsize=(14, 7))
    # add master title for whole figure
    fig.suptitle(metric_to_label[metric], fontsize=FONTSIZE)

    # loop across vendors
    for vendor, row in df_vendor.iterrows():
        # loop across roi/labels
        for index, label in enumerate(roi_to_label.keys()):
            # create individual subplots
            ax = plt.subplot(2, 3, index + 1)
            # loop across levels
            y = list()
            e = list()
            for level in levels_to_label.keys():
                # get mean value for given label (e.g, C2, C3, etc) and given label/roi (e.g., spinal cord etc.)
                y.append(row[level,label][0])
                # get std value for given label (e.g, C2, C3, etc) and given label/roi (e.g., spinal cord etc.)
                e.append(row[level, label][1])

            # plot mean values pervendor for each level (C2, C3, C4, C5)
            x = [float(key) for key in levels_to_label]  # individual levels - 2,3,4,5
            plt.errorbar(x,
                         y,
                         e,
                         marker=vendor_to_marker[vendor],  # change marker symbol based on vendor
                         markersize=8,  # size of marker symbol
                         capsize=5,     # the length of the error bar caps
                         alpha=0.5,  # transparency
                         label=vendor)
            # rename xticks to C2, C3, C4, C5
            plt.xticks(x, levels_to_label.values())
            # add grid
            plt.grid(axis='y', linestyle="--")
            # add title to individual subpolots (i.e., ROI/label)
            plt.title(roi_to_label[label])

    # place legend next to last subplot
    leg = plt.legend(bbox_to_anchor=(1.1, 1), fontsize=FONTSIZE - 5)
    # insert number of subjects and number of sites per vendor into legend
    # loop across vendors
    for num in range(0,len(leg.get_texts())):
        leg.get_texts()[num].set_text('{}: {} subjects, {} sites'.
                                    format(df_summary_vendor.index.values[num], df_summary_vendor.iloc[num, 0],
                                           df_summary_vendor.iloc[num, 1]))

    # Move subplots closer to each other
    plt.subplots_adjust(wspace=-0.5)
    plt.tight_layout()
    # tight layout of whole figure and shift master title up
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    # save figure
    fname_fig = os.path.join(path_output, metric + '_per_vendor.png')
    plt.savefig(fname_fig, dpi=200)
    logger.info('Created: ' + fname_fig)


def load_participants_file(path_participants):
    """
    Build Pandas DF of participants with their characteristics (age, sex, ...) based on participants.tsv file
    :param path_participants: path to participants.tsv file
    :return participants_df: pandas df
    """
    if os.path.isfile(path_participants):
        participants_df = pd.read_csv(path_participants, sep="\t")
    else:
        raise FileNotFoundError("File \"participants.tsv\" was not found in {} folder. You can specify path to this "
                                "file by -participants-file flag.".format(os.getcwd()))

    return participants_df

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
        help='Path to directory with results (*.csv files).')
    parser.add_argument(
        '-participants-file',
        required=False,
        metavar='<fname>',
        help='Path to .tsv file with participants characteristics.')
    parser.add_argument(
        '-config',
        required=False,
        metavar='<fname>',
        help='Config yaml file listing images that you want to remove from processing.')

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

    # fetch where is participants.tsv file located
    if args.participants_file is not None:
        if os.path.isfile(args.participants_file):
            path_participants = args.participants_file
        else:
            raise FileNotFoundError("Participants file '{}' was not found.".format(args.participants_file))
    else:
        # if not passsed, assuming it is located in same dir as a *csv files
        path_participants = os.path.join(os.getcwd(), 'participants.tsv')

    # fetch perlevel .csv files
    csv_files = glob.glob('*perlevel.csv')

    if not csv_files:
        raise RuntimeError("No *.csv files were found in the current directory. You can specify directory with *.csv "
                           "files by -path-results flag.")

    # Dump log file there
    if os.path.exists(FNAME_LOG):
        os.remove(FNAME_LOG)
    fh = logging.FileHandler(os.path.join(os.path.abspath(os.curdir), FNAME_LOG))
    logging.root.addHandler(fh)

    # loop across individual .csv files (i.e., across individual qMRI metrics)
    for csv_file in csv_files:

        logger.info('\nProcessing: ' + csv_file)
        # open .csv file and create dict
        dict_results = []
        with open(csv_file, newline='') as f_csv:
            reader = csv.DictReader(f_csv)
            for row in reader:
                dict_results.append(row)

        # fetch metric name from .csv file
        _, csv_file_small = os.path.split(csv_file)
        metric = file_to_metric[csv_file_small]

        # fetch metric values per site
        results_dict = aggregate_per_site(dict_results, metric, dict_exclude_subj, path_participants)
        # make it a pandas structure (easier for manipulations)
        df = pd.DataFrame.from_dict(results_dict, orient='index')

        # ------------------------------------------------------------------
        # compute per vendor summary
        # ------------------------------------------------------------------
        df_vendor_mean, df_vendor_std, df_summary_vendor = summary_per_vendor(df)
        # and save mean as a .csv file
        fname_csv_per_vendor = os.path.join(os.getcwd(), metric) + '_per_vendors.csv'
        df_vendor_mean.to_csv(fname_csv_per_vendor)
        logger.info('Created: ' + fname_csv_per_vendor)

        logger.info(df_vendor_mean.head())

        # ------------------------------------------------------------------
        # compute age and sex per vendor
        # ------------------------------------------------------------------
        df_age, df_sex = aggregate_age_and_sex_per_vendor(path_participants)

        # Concatenate number of sites and subjects with sex and age pervendor
        df_summary_vendor = pd.concat([df_summary_vendor, df_age, df_sex], sort=False, axis=1)

        logger.info(df_summary_vendor.to_string(index=False))

        # ------------------------------------------------------------------
        # generate per VENDORs figure - level evolution per ROI for individual vendors
        # ------------------------------------------------------------------

        # concatenate mean and std values
        df_vendor = pd.concat([df_vendor_mean, df_vendor_std], axis=1)
        generate_level_evolution_pervendor(df_vendor, df_summary_vendor, metric, path_output)

        # ------------------------------------------------------------------
        # generate per SITEs figure - level evolution per ROI for individual sites
        # ------------------------------------------------------------------

        generate_level_evolution_persite(df, metric, path_output)

        print('Finished.')

if __name__ == "__main__":
    main()
