#!/usr/bin/env python

#
# Generate figures and compute statistics for individual qMRI metrics (FA, MD, AD, RD, MTR, MTsat) for different
# labels/ROIs (SC, WM, GM, ...) along individual cervical levels (C2, C3, C4, C5) pervendor and persite
#
# USAGE:
#
#   python generate_figures.py
#   -path-results ~/spineGeneric-multi-subject_results/results/perlevel
#   -config ~/spineGeneric-multi-subject_results/results/exclude.yml
#   -participants-file ~/spineGeneric-multi-subject_results/results/participants.tsv
#
# Input arguments:
#   -path-results           - directory with perlevel *.csv files (computed by extract_normative_metrics.py script)
#   -config                 - input yml config file with subjects to exclude (e.g., due to bad data quality, noise, ...)
#   -participants-file      - input .tsv file with participants characteristics (sex, age, ...)
#
# Inspired by - https://github.com/sct-pipeline/spine-generic/blob/master/processing/generate_figure.py
# Authors: Jan Valosek, Julien Cohen-Adad
#

# TODO - combine this script (and probably also whole repo) with spine-generic repo

import os
import argparse
import tqdm
import sys
import glob
import csv
import pandas as pd
import yaml

import numpy as np
import matplotlib.pyplot as plt
import logging

from collections import defaultdict
from scipy.stats import f_oneway
from matplotlib.lines import Line2D


# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)

FNAME_LOG = 'log_stats.txt'
log_line = '========================================================'

# color to assign to each MRI model for the figure
vendor_to_color = {
    'GE': 'orange',
    'Philips': 'dodgerblue',
    'Siemens': 'limegreen',
    }

# change marker for individual vendors for scatter plots
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

# fetch metric field (column name in .csv file)
metric_to_field = {
    'dti_fa': extraction_method,
    'dti_md': extraction_method,
    'dti_ad': extraction_method,
    'dti_rd': extraction_method,
    'mtr': extraction_method,
    'mtsat': extraction_method,
    }

# abbreviation for log
metric_to_abbreviation = {
    'dti_fa': 'FA',
    'dti_md': 'MD',
    'dti_ad': 'AD',
    'dti_rd': 'RD',
    'mtr': 'MTR',
    'mtsat': 'MTsat',
    }

# master titles for figures with subplots
metric_to_title = {
    'dti_fa': 'Fractional anisotropy (FA)',
    'dti_md': 'Mean diffusivity (MD)',
    'dti_ad': 'Axial diffusivity (AD)',
    'dti_rd': 'Radial diffusivity (RD)',
    'mtr': 'Magnetization transfer ratio (MTR)',
    'mtsat': 'Magnetization transfer saturation (MTsat)',
    }

# ylabel for subplots
metric_to_label = {
    'dti_fa': 'FA',
    'dti_md': 'MD [$× 10^{-3} mm^{2}/s$]',
    'dti_ad': 'AD [$× 10^{-3} mm^{2}/s$]',
    'dti_rd': 'RD [$× 10^{-3} mm^{2}/s$]',
    'mtr': 'MTR [%]',
    'mtsat': 'MTsat [%]',
    }

# titles for individual subplots
# TODO - add also another ROIs/labels to the dict below
roi_to_label = {
    'spinal cord': 'Spinal cord',
    'white matter': 'White matter',
    'gray matter': 'Gray matter',
    'dorsal columns': 'Dorsal columns',
    'lateral funiculi': 'Lateral funiculi/columns',
    'ventral funiculi': 'Ventral funiculi/columns',
}

# xticks in subplots
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
    Aggregate qMRI metrics per site. Loop across lines in dict_results (loaded rows from input .csv file) and build
    nested dict with aggregated metrics persite
    :param dict_results: dictionary with loaded rows from input .csv file for given qMRI metric
    :param metric: currently processed metric (e.g., dti_md)
    :param dict_exclude_subj: dictionary with subjects to exclude from analysis (due to bad data quality, etc.)
    :param path_participants: path to participants.tsv file
    :return: results_agg: nested dict with metric values persite
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
        subject = fetch_subject(filename)       # e.g., sub-beijingGE01
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
                # initialize empty sub-list for processed subjects for given site
                if not 'processed_subjects' in results_agg[site].keys():
                    results_agg[site]['processed_subjects'] = list()
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
                # get vertlevel for site
                vertlevel = dict_results[i]['VertLevel']  # e.g., 5
                if not val == 'None':
                    # append currently processed subject into list to have information which subjects were analysed
                    if not subject in results_agg[site]['processed_subjects']:
                        results_agg[site]['processed_subjects'].append(subject)
                    val = float(dict_results[i][metric_field]) * scaling_factor[metric]     # scale metric
                    # append data into sub-dict  - {'vertlevel' 'label': 'metric values'} (key is tuple, value is list)
                    results_agg[site]['val'][(vertlevel, label)].append(float(val))
                # if value is missing, report it to log
                else:
                    logger.info('Value for {} at level {} in {} is None or missing.'.format(subject, vertlevel, label))
                # compute mean perlevel per ROI/label - {'vertlevel' 'label': 'mean value'} (key is tuple, value is float)
                results_agg[site]['mean'][(vertlevel, label)] = np.mean(results_agg[site]['val'][(vertlevel, label)])
            else:
                logger.info('Skipping {}.'.format(label))

        else:
            # because we iterate across individual lines, same subject can be included more than one time in
            # subjects_removed list -> use conversion to set in next command
            subjects_removed.append(subject)

    logger.info('\n{} subjects were removed: {}\n'.format(len(set(subjects_removed)),set(subjects_removed)))

    return results_agg


def aggregate_age_and_sex_per_vendor(path_participants, subjects_processed):
    """
    Aggregate age and sex per individual vendors for successfully processed subjects (i.e., without excluded subjects)
    :param path_participants: path to participants.tsv file
    :param subjects_processed: list of all successfully processed subjects
    :return: df_age: pandas df with age statistics (min, max, mean, std) per vendor
    :return: df_sex: pandas df with sex (M, F) per vendor
    """

    participants_df = load_participants_file(path_participants)

    # Replace '-' by NaN (some subjects do not have listed age and sex)
    participants_df = participants_df.replace('-', np.nan)
    # Convert age from str to float to fit with NaN
    participants_df['age'] = participants_df['age'].astype(float)

    # Let only successfully processed participants
    processed_participants_df = participants_df[participants_df['participant_id'].isin(subjects_processed)]

    # Aggregate age per vendors
    df_age = processed_participants_df.groupby('manufacturer').age.agg(['min', 'max', 'mean', np.std])
    df_age.columns = ['min age', 'max age', 'mean age', 'std age']
    df_age = df_age.round(1)    # round to one decimal

    dict_sex = defaultdict(tuple)
    # Loop across subjects grouped by vendors (i.e, 3 groups - GE, Philips, Siemens)
    for vendor, value in processed_participants_df.groupby('manufacturer'):
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

    # ANOVA between vendors in age
    stat, pvalue = f_oneway(processed_participants_df[processed_participants_df['manufacturer'] == 'Siemens'].age,
                            processed_participants_df[processed_participants_df['manufacturer'] == 'Philips'].age,
                            processed_participants_df[processed_participants_df['manufacturer'] == 'GE'].age,
                            )

    logger.info('\nANOVA between vendors in age: p{}\n'.format(format_pvalue(pvalue)))

    return df_age, df_sex


def summary_per_vendor(df, metric):
    """
    Compute mean, sd and cov values pervendor (Siemens, GE, Philips) for individual labels/ROIs (SC, WM, GM, ...)
    perlevel (C2, C3, C4, C5) and number of used (so, after exclusion) subjects and sites per vendor.
    :param df: pandas df with metric values per site
    :param metric: currently processed metric (e.g., dti_fa)
    :return: df_vendor_mean: pandas df with mean metric values per vendor
    :return: df_vendor_sd: pandas df with sd metric values per vendor
    :return: df_vendor_mean_and_sd: pandas df with mean and sd as a single string (great for Excel import)
    :return: df_vendor_cov: pandas df with COV values per vendor
    :return: df_summary_vendor: pandas df with number of subjects and sites per vendor
    """

    # compute mean and sd values per vendor (GE, Philips, Siemens) for individual ROI/labels perlevel
    # initialize dict
    dict_vendor_mean = {}
    dict_vendor_sd = {}
    dict_vendor_cov = {}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        # if this is a new vendor, initialize sub-dict
        dict_vendor_mean[vendor] = {}
        dict_vendor_sd[vendor] = {}
        dict_vendor_cov[vendor] = {}
        # loop across levels and ROI (e.g., '5', 'spinal cord'; '4', 'spinal cord'; ...)
        for label in list(df['mean']['amu'].keys()):        # <class 'tuple'>: (e.g., '5', 'ventral funiculi')
            # if this is a new label, initialize sub-dict
            dict_vendor_mean[vendor][label] = {}
            dict_vendor_sd[vendor][label] = {}
            dict_vendor_cov[vendor][label] = {}
            mean_values = list()
            # loop across sites for given vendor
            for site in df[df['vendor'] == vendor]['site']:
                # collect mean values from all sites for given vendor to one list
                mean_values.append(float(df['mean'][site][label]))
            # compute mean from all sites' mean values for given vendor and insert it into dict
            dict_vendor_mean[vendor][label] = round(np.mean(mean_values),2)
            # compute std from all sites' mean values for given vendor and insert it into dict
            dict_vendor_sd[vendor][label] = round(np.std(mean_values),2)
            # compute COV in %
            dict_vendor_cov[vendor][label] = round((np.std(mean_values)/np.mean(mean_values))*100,2)

    # Create dict with mean and SD as a single string (e.g, '0.73+-0.07') for each level and ROI
    # Great for saving as a .csv file and importing as a table into Excel
    # This dict is then converted to pandas DF and saved as csv table
    dict_vendor_mean_and_sd = {}    # initialize dict
    for vendor in dict_vendor_mean.keys():      # loop across vendors
        dict_vendor_mean_and_sd[vendor] = {}        # initialize sub-dict for each vendor
        for key, value in dict_vendor_mean[vendor].items():     # loop across items (e.g., '5', 'spinal cord')
            # combine mean and SD values into single strings (e.g., '0.73+-0.07')
            dict_vendor_mean_and_sd[vendor][key] = (str(dict_vendor_mean[vendor][key]) +        # mean value
                                                    '\u00B1' +                                   # +- sign
                                                    str(dict_vendor_sd[vendor][key]))            # sd value


    # convert dict to pandas df for easy csv save
    df_vendor_mean = pd.DataFrame.from_dict(dict_vendor_mean, orient='index')
    df_vendor_sd = pd.DataFrame.from_dict(dict_vendor_sd, orient='index')
    df_vendor_cov = pd.DataFrame.from_dict(dict_vendor_cov, orient='index')
    df_vendor_mean_and_sd = pd.DataFrame.from_dict(dict_vendor_mean_and_sd, orient='index')

    # find out which ROI/label and vertebral level has the largest COV and print it to the log
    logger.info('')    # print empty line to log
    cov_max = df_vendor_cov.max(axis=1)             # get max across all ROI and levels
    cov_max_col = df_vendor_cov.idxmax(axis=1)      # get column name (e.g., '3, spinal cord')
    # loop across vendors
    for vendor in vendor_to_color.keys():
        logger.info('{} showed the largest COV for {} in {} at {} level ({}%).'.format(
                    metric_to_abbreviation[metric],  # metric (e.g., FA)
                    vendor,                  # vendor (e.g., Siemens)
                    cov_max_col[vendor][1],  # label (e.g., spinal cord)
                    cov_max_col[vendor][0],  # vertebral level (e.g., 3)
                    cov_max[vendor],         # COV value in %
                    ))

    # compute number of subjects pervendor
    dict_sub_per_vendor = dict()        # e.g.: {'GE': 7, 'Philips': 2, 'Siemens': 9}
    # loop across vendors
    for vendor in vendor_to_color.keys():
        num_of_sub = 0
        # loop across sites for given vendor
        for site in df[df['vendor'] == vendor]['site']:
            # get number of used subjects for given site and add it to num_of_sub variable
            num_of_sub = num_of_sub + len(df[df['vendor'] == vendor]['processed_subjects'][site])

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

    return df_vendor_mean, df_vendor_sd, df_vendor_mean_and_sd, df_vendor_cov, df_summary_vendor


def generate_level_evolution_persite(df, df_summary_vendor, metric, path_output):
    """
    Generate figure for each metric - level evolution (C2, C3, C4, C5) per ROI for individual sites
    :param df: pandas df with input data
    :param: df_summary_vendor: pandas df with number of subjects and sites per vendor
    :param metric: currently processed qMRI metric (e.g. dti_md)
    :param path_output: path where figures will be generated (manually entered -path-results path or current dir)
    :return:
    """
    # get individual sites
    site_sorted = df.sort_values(by=['vendor', 'model', 'site']).index.values

    fig, _ = plt.subplots(figsize=(17, 10))
    # add master title for whole figure
    fig.suptitle('{} for totally {} subjects across {} sites.'.format(metric_to_title[metric],
                                                                      df_summary_vendor['M'].sum() +
                                                                      df_summary_vendor['F'].sum(),
                                                                      len(site_sorted)), fontsize=FONTSIZE)

    # loop across sites
    for site in site_sorted:
        # initialize dict for each loop - {label: mean_metric}
        mean_dict = dict()
        # loop across roi/labels
        for index, label in enumerate(roi_to_label.keys()):
            # create individual subplots
            # TODO - modify command below to be able to work with different number of ROIs
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

            # # show legend only one time (in right up corner)
            # if index == 2:
            #     # place legend next to last subplot
            #     leg = plt.legend(bbox_to_anchor=(1.1, 1.03), fontsize=FONTSIZE - 5)
            #     # insert number of subjects for individual sites into legend
            #     # loop across sites
            #     for num in range(0,len(leg.get_texts())):
            #         leg.get_texts()[num].set_text('{}: {} subjects'.format(site_sorted[num],        # site
            #                                                                len(df[df['site'] == site_sorted[num]]['processed_subjects'].values[0])))    # number of subjects for given site

    # TODO - find how to show in legend each site only once
    leg = fig.legend()
    leg.texts = leg.get_texts()[:len(site_sorted)]
    # loop across sites
    for num in range(0,len(site_sorted)):
        leg.texts[num].set_text('{}: {} subjects'.format(site_sorted[num],        # site
                                                               len(df[df['site'] == site_sorted[num]]['processed_subjects'].values[0])))


    # Move subplots closer to each other
    plt.subplots_adjust(wspace=-0.5)
    plt.tight_layout()
    # tight layout of whole figure and shift master title up
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    # save figure
    fname_fig = os.path.join(path_output, 'figures', metric + '_per_sites.png')
    plt.savefig(fname_fig, dpi=200)
    logger.info('\nCreated: {}\n'.format(fname_fig))

    # plt.show()


def generate_level_evolution_pervendor(df_vendor, df_summary_vendor, metric, path_output):
    """
    Generate figure for each metric - level evolution (C2, C3, C4, C5) per ROI for individual vendors
    :param df_vendor: pandas df with aggregated data pervendor
    :param df_summary_vendor: pandas df with number of subjects and sites per vendor
    :param metric: currently processed qMRI metric (dti_fa)
    :param path_output: path where figures will be generated (manually entered -path-results path or current dir)
    :return:
    """

    # optimize font size
    FONTSIZE = 20
    TICKSIZE = 15
    LEGENDSIZE = 12
    # transparency of artists
    ALPHA = 0.6

    # Initialize figure for subplots
    # TODO - modify command below to allow working with different number of ROIs/labels
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    # Flatten 2D array into 1D to allow iteration by loop
    axs = axs.ravel()

    # add master title for whole figure
    fig.suptitle('{} for {} subjects across {} vendors'.
                 format(metric_to_title[metric],
                        df_summary_vendor['M'].sum() +
                        df_summary_vendor['F'].sum(),
                        len(vendor_to_color)),
                 fontsize=FONTSIZE, fontweight='bold')

    # loop across vendors (Siemens, Philips, ...)
    for vendor, row in df_vendor.iterrows():
        # loop across roi/labels (spinal cord, white matter, ...)
        for index, label in enumerate(roi_to_label.keys()):

            y = list()      # mean values
            e = list()      # sd values
            # loop across levels (C2, C3, ...)
            for level in levels_to_label.keys():
                # get mean value for given label (e.g, C2, ...) and given label/roi (e.g., spinal cord, ...)
                y.append(row[level,label][0])
                # get sd value for given label (e.g, C2, ...) and given label/roi (e.g., spinal cord, ...)
                e.append(row[level, label][1])

            # plot mean and sd values pervendor for each level (C2, C3, C4, C5)
            x = [float(key) for key in levels_to_label]  # individual levels - 2,3,4,5
            axs[index].errorbar(x,      # vertebral levels
                                y,      # mean values for currently processed qMRI metric (e.g., FA, ...)
                                e,      # sd values for currently processed qMRI metric (e.g., FA, ...)
                                marker=vendor_to_marker[vendor],  # change marker symbol based on vendor
                                markersize=10,     # size of marker symbol
                                capsize=5,         # the length of the error bar caps
                                alpha=ALPHA,         # transparency
                                color=vendor_to_color[vendor],
                                label=vendor)      # label for legend
            # rename xticks to C2, C3, C4, C5
            #axs[index].xticks(x, levels_to_label.values(), fontsize=TICKSIZE)
            # Rename x-tick labels to C2, C3, C4, C5
            plt.setp(axs[index], xticks=x, xticklabels=levels_to_label.values())
            # Increase size of x- and y-ticks
            plt.setp(axs[index].xaxis.get_majorticklabels(), fontsize=TICKSIZE)
            plt.setp(axs[index].yaxis.get_majorticklabels(), fontsize=TICKSIZE)
            # add grid
            axs[index].grid(axis='y', linestyle="--")
            # add title to individual subpolots (i.e., ROI/label)
            axs[index].set_title(roi_to_label[label], fontsize=FONTSIZE)
            # set y-label (FA, MD, ...) only once for each row
            # TODO - number of indexes will have to be fixed when number of ROIs/labels will be changed
            if index == 0 or index == 3:
                axs[index].set_ylabel(metric_to_label[metric], fontsize=TICKSIZE)


    # LEGEND - create custom legend
    # https://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    # https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    lines = list()      # initialize list for individual symbols in the legend
    labels = list()     # initialize list for individual text labels in the legend
    # loop across vendors (Siemens, Philips, ...)
    for num, vendor in enumerate(vendor_to_color):
        lines.append(Line2D([0], [0], color=vendor_to_color[vendor],
                            marker=vendor_to_marker[vendor],
                            markersize=LEGENDSIZE,
                            alpha=ALPHA,
                            linestyle=''))
        labels.append('{}: {} subjects, {} sites'.format(df_summary_vendor.index.values[num],
                                                         df_summary_vendor.iloc[num, 0],
                                                         df_summary_vendor.iloc[num, 1]))

    # Move subplots closer to each other
    plt.subplots_adjust(wspace=-0.5)
    plt.tight_layout()

    # Insert legend below subplots, NB - this line has to be below the plt.tight_layout()
    legend = fig.legend(lines, labels, loc='lower left', bbox_to_anchor=(0.2, 0),
                        bbox_transform=plt.gcf().transFigure, ncol=len(lines), fontsize=LEGENDSIZE)
    # Change box's frame color to black
    frame = legend.get_frame()
    frame.set_edgecolor('black')

    # tight layout of whole figure and shift master title up
    fig.tight_layout()
    fig.subplots_adjust(top=0.88, bottom=0.1)

    # save figure
    fname_fig = os.path.join(path_output, 'figures', metric + '_per_vendor.png')
    plt.savefig(fname_fig, dpi=200)
    logger.info('\nCreated: ' + fname_fig)


def format_pvalue(p_value, alpha=0.05, include_equal=True):
    """
    If p-value is lower than 0.05, change it to "<0.05", otherwise, round it to two decimals
    :param p_val: input p-value as a float
    :param alpha: significance level
    :param include_equal: include equal sign ('=') to pvalue (e.g., '=0.06') or not (e.g., '0.06')
    :return: p_val: processed p-value (replaced by "<0.05" or rounded to two decimals) as a str
    """
    if p_value < alpha:
        p_value = "<" + str(alpha)
    else:
        if include_equal:
            p_value = '=' + str(round(p_value, 3))
        else:
            p_value = str(round(p_value, 3))

    return p_value


def check_consistency(results_dict, path_participants, csv_file):
    """
    Check if number of subjects in participants.tsv file corresponds with number of subject in input .csv file
    :param results_dict: dict with input metric data based on input csv_file
    :param path_participants: path to participants.tsv file
    :param csv_file: currently processed .csv file (e.g., DWI_MD_perlevel.csv)
    :return: None:
    """

    participants_df = load_participants_file(path_participants)

    # loop across individual sites
    for site in results_dict.keys():
        # get number of subject from input csv file with metric (e.g., DWI_MD_perlevel.csv)
        num_of_sub_csv_file = len(results_dict[site]['processed_subjects'])
        # get number of subject from input tsv participants file (participants.tsv)
        num_of_sub_participants_file = len(participants_df[participants_df['institution_id'] == site])
        if num_of_sub_csv_file != num_of_sub_participants_file:
            logger.info('site {}: {} subjects were found in input {} file but {} subjects were found in participants.tsv file'.
                        format(site,
                               num_of_sub_csv_file,
                               csv_file,
                               num_of_sub_participants_file))


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
    Get subject ID from filename
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
    :param dict_exclude_subj: dictionary with subjects to exclude from analysis (due to bad data quality, etc.)
    :return: Bool
    """
    if metric in dict_exclude_subj.keys():
        # remove single subject (e.g., sub-geneva02)
        if subject in dict_exclude_subj[metric]:
            return True
        # remove all subjects for given site (e.g., fslAchieva)
        elif subject[4:-2] in dict_exclude_subj[metric]:    # subject[4:-2] extract only site from subject
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
        help='Path to .tsv file with participants characteristics (sex, age, ...).')
    parser.add_argument(
        '-config',
        required=False,
        metavar='<fname>',
        help='Path to yaml file listing subjects and sites that you want to exclude from the processing.')

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

    # create directory figures where created figures will be saved
    if not os.path.exists(os.path.join(path_output, 'figures')):
        os.makedirs(os.path.join(path_output, 'figures'))

    # create directory tables where created tables will be saved
    if not os.path.exists(os.path.join(path_output, 'tables')):
        os.makedirs(os.path.join(path_output, 'tables'))

    # fetch where participants.tsv file is located
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

        logger.info('\n{}\nProcessing: {}\n{}'.format(log_line, csv_file, log_line))
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

        # fetch all successfully processed subjects
        subjects_processed = [item for sublist in list(df['processed_subjects'].values) for item in sublist]
        logger.info('{} subjects were processed: {}\n'.format(len(subjects_processed), subjects_processed))
        check_consistency(results_dict, path_participants, csv_file)

        # ------------------------------------------------------------------
        # compute per vendor summary
        # ------------------------------------------------------------------
        df_vendor_mean, df_vendor_sd, df_vendor_mean_and_sd, df_vendor_cov, df_summary_vendor = summary_per_vendor(df, metric)
        #  Save mean_and_sd tables and cov as a .csv files
        fname_csv_per_vendor_mean_sd = os.path.join(os.getcwd(), 'tables', metric + '_mean_and_sd_per_vendors.csv')
        df_vendor_mean_and_sd.to_csv(fname_csv_per_vendor_mean_sd)
        logger.info('\nCreated: {}'.format(fname_csv_per_vendor_mean_sd))
        fname_csv_per_vendor_cov = os.path.join(os.getcwd(), 'tables', metric + '_cov_per_vendors.csv')
        df_vendor_cov.to_csv(fname_csv_per_vendor_cov)
        logger.info('\nCreated: {}'.format(fname_csv_per_vendor_cov))

        # ------------------------------------------------------------------
        # compute age and sex per vendor
        # ------------------------------------------------------------------
        df_age, df_sex = aggregate_age_and_sex_per_vendor(path_participants, subjects_processed)

        # Concatenate number of sites and subjects with sex and age pervendor
        df_summary_vendor = pd.concat([df_summary_vendor, df_age, df_sex], sort=False, axis=1)

        logger.info(df_summary_vendor.to_string(index=False))

        # ------------------------------------------------------------------
        # generate per VENDORs figure - level evolution per ROI for individual vendors
        # ------------------------------------------------------------------

        # concatenate mean and std values
        df_vendor = pd.concat([df_vendor_mean, df_vendor_sd], axis=1)
        generate_level_evolution_pervendor(df_vendor, df_summary_vendor, metric, path_output)

        # ------------------------------------------------------------------
        # generate per SITEs figure - level evolution per ROI for individual sites
        # ------------------------------------------------------------------

        #generate_level_evolution_persite(df, df_summary_vendor, metric, path_output)

        print('\n Analysis of {} is completed.'.format(csv_file))

if __name__ == "__main__":
    main()
