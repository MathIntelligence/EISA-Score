#!/usr/bin/env python

"""
Introduction:
    Get eisa features using local surface area

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    May 26, 2022

"""

import argparse
import os
import sys
sys.path.append('../src')
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2

from eisa_score_local import EISA_Score_Local_Surface
import time




def main(args):

	kernel_type = args.kernel_type
	kernel_tau = args.kernel_tau
	kernel_power = args.kernel_power
	cutoff = args.cutoff
	isovalue = args.isovalue
	mesh_size = args.mesh_size

	data_folder = args.data_folder
	out_dir = args.out_dir 

	dataset_csv_file = args.dataset_csv_file
	pdbid_index = args.pdbid_index

	df_pdbids = pd.read_csv(dataset_csv_file)
	pdbids = df_pdbids['PDBID'].tolist()
	
	pdbid = pdbids[pdbid_index]




	eisa_class = EISA_Score_Local_Surface(data_folder,
	                                      pdbid, 
										  kernel_type=kernel_type,
		                                  kernel_tau=kernel_tau,
		                                  kernel_power=kernel_power,
		                                  cutoff=cutoff,
		                                  isovalue=isovalue,
		                                  mesh_size=mesh_size)
	
	features = eisa_class.get_features()

	output_file_name = f'{out_dir}/local-type-{kernel_type}'\
	                     f'-tau-{kernel_tau}'\
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

	parser.add_argument('--isovalue', type=float, action='store',
	                    default=0.25, help='isovalue')

	parser.add_argument('--out_dir', type=str, action='store',
	                    help='output directory')

	parser.add_argument('--data_folder', type=str, action='store',
	                    help='dataset folder path')

	parser.add_argument('--dataset_csv_file', type=str, action='store',
	                    help='dataset csv file containing PDBID and pK values')

	parser.add_argument('--pdbid_index', type=int, action='store',
	                    help='index of PDBID of the molecule in the dataset')

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
