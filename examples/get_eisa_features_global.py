#!/usr/bin/env python

"""
Introduction:
    Get eisa features using global surface area

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    May 26, 2022

"""

import argparse
import os
import sys
sys.path.append('/home/mra309/Projects/EISA-Score/src')
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2

import numba_utils_methods as nmb
from eisa_score_global import EISA_Score_Global_Surface
import time



def main(args):

	kernel_type = args.kernel_type
	kernel_tau = args.kernel_tau
	kernel_power = args.kernel_power
	cutoff = args.cutoff
	mesh_size = args.mesh_size

	data_folder = args.data_folder
	out_dir = args.out_dir 

	pdbid = args.pdbid




	eisa_class = EISA_Score_Global_Surface(data_folder,
	                                      pdbid, 
										  kernel_type=kernel_type,
		                                  kernel_tau=kernel_tau,
		                                  kernel_power=kernel_power,
		                                  cutoff=cutoff,
		                                  mesh_size=mesh_size)
	
	features = eisa_class.get_features()

	output_file_name = f'{out_dir}/global-tau-{kernel_tau}'\
	                      f'-power-{kernel_power}-cutoff-{cutoff}'\
	                     f'-pdbid-{pdbid}'

	np.save(output_file_name, features)


def parse_args(args):
	parser = argparse.ArgumentParser(description="Get EISA Features")

	parser.add_argument('--mesh_size', type=float, action='store',
	                    default=1.5, help='mesh size')

	parser.add_argument('--cutoff', type=float, action='store',
	                    default=5, help='cutoff value')

	parser.add_argument('--kernel_type', type=str, action='store',
	                    default='exponential', help='kernel type')

	parser.add_argument('--kernel_tau', type=float, action='store',
	                    default=1.0, help='scale parameter')

	parser.add_argument('--kernel_power', type=float, action='store',
	                    default=2.0, help='Kernel parameter')

	parser.add_argument('--out_dir', type=str, action='store',
	                    help='output directory')

	parser.add_argument('--data_folder', type=str, action='store',
	                    help='dataset folder path')

	parser.add_argument('--pdbid', type=str, action='store',
	                    help='pdb id of the molecule')

	args = parser.parse_args()

	return args


def cli_main():
	args = parse_args(sys.argv[1:])

	print(args)

	main(args)



if __name__ == "__main__":

	t0 = time.time()

	cli_main()

	print('Done!')
	print('Elapsed time: ', time.time()-t0)
