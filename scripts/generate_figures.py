#!/usr/bin/env python
#
# Generate figures for individual qMRI metrics
#
# Authors: Julien Cohen-Adad, Jan Valosek

import os
import argparse
import tqdm
import sys
import glob


# List of sites to exclude based on the metric
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

# FIGURE PARAMETERS
FONTSIZE = 15
TICKSIZE = 10
LABELSIZE = 15


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
        raise RuntimeError("Variable 'csv_files' is empty, i.e. no *.csv files were found in current directory.")

    print(csv_files)

if __name__ == "__main__":
    main()