#!/usr/bin/env python

"""
Introduction:
    Example generating eisa features

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    May 26, 2022

"""

import argparse
from subprocess import call
from sys import exit
import numpy as np
import pandas as pd

# Specify the path of the PDBbind dataset
data_folder = '/scratch/mra309/PdbbindDataSets/v2007'

# Specify the output directory for the features
out_dir = '/scratch/mra309/projects/eisa_score/features'

# Specify the csv data file containing 'PDBID' and 'pK' values
# from utils/csv_data_files
dataset_csv_file_path = '../utils/csv_data_files/PDBbindv2007_RefinedSet.csv'

df_pdbids = pd.read_csv(dataset_csv_file_path)
pdbids = df_pdbids['PDBID'].tolist()


def run_args_local(index):
    argPairs = {
        'kernel_type': 'exponential',
        'kernel_tau': 0.5,
        'kernel_power': 15.0,
        'cutoff': 7.0,
        'isovalue': 0.25,
        'mesh_size': 0.5,
        'pdbid': pdbids[index],
        'data_folder': data_folder,
        'out_dir': out_dir
    }

    args = ['./get_eisa_features_local.py']

    for k in argPairs.keys():
        args.append('--%s=%s' % (k, argPairs[k]))
    # print(args)

    rtnCode = call(args)

    if rtnCode != 0:
        print('Code failed!')
        exit(-1)


def run_args_global(index):
    argPairs = {
        'kernel_type': 'exponential',
        'kernel_tau': 0.5,
        'kernel_power': 15.0,
        'cutoff': 12.0,
        'mesh_size': 0.5,
        'pdbid': pdbids[index],
        'data_folder': data_folder,
        'out_dir': out_dir
    }

    args = ['./get_eisa_features_global.py']

    for k in argPairs.keys():
        args.append('--%s=%s' % (k, argPairs[k]))
    # print(args)

    rtnCode = call(args)

    if rtnCode != 0:
        print('Code failed!')
        exit(-1)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate EISA Features")

    parser.add_argument('--pdb_index', type=int, default=1,
                        help='index of pdbid')

    parser.add_argument('--surface_type', type=str, default='global',
                        help='index of pdbid')

    args = parser.parse_args()

    pdb_index = args.pdb_index
    surface_type = args.surface_type

    if surface_type[0] in ['g', 'G']:
        run_args_global(pdb_index)

    elif surface_type[0] in ['l', 'L']:
        run_args_local(pdb_index)

    else:
        print("Warning: invalid surface type!")
