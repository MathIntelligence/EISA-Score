#!/usr/bin/env python

"""
Introduction:
    EISA score global surface area

Author:
    Masud Rana (masud.rana@uky.edu)
    
Date last modified:
    May 26, 2022

"""

import os
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2
from rdkit import Chem
from scipy.spatial.distance import cdist

import numba_utils_methods as nmb
from element_interactive_density import ElementInteractiveDensity

import time


class EISA_Score_Global_Surface():
    def __init__(self, path, pdbid, kernel_type='exponential',
                 kernel_tau=1.0, kernel_power=2.0, cutoff=12.0):

        self.path = path  # path to pdbbind dataset
        self.pdbid = pdbid
        self.cutoff = cutoff

        self.kernel_type = kernel_type
        self.kernel_tau = kernel_tau
        self.kernel_power = kernel_power

        self.protein_atom_type = ['C', 'N', 'O', 'S']
        self.ligand_atom_type = ['H', 'C', 'N',
                                 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']

        self.atom_type_radii = {
            'H': 1.2,
            'C': 1.7,
            'N': 1.55,
            'O': 1.52,
            'S': 1.8,
            'P': 1.8,
            'F': 1.8,
            'Cl': 1.47,
            'Br': 1.75,
            'I': 1.98
        }

        self.num_stat_measure = 6  # sum,mean,median,std,max,and min
        self.mesh_size = 0.5

        self.isovalue_list = np.arange(0.05, 0.8, 0.05)

    def sdf_to_df(self, sdf_file):

        m = Chem.MolFromMolFile(sdf_file, sanitize=False)
        m.UpdatePropertyCache(strict=False)

        lines = []
        for atom in m.GetAtoms():
            if atom.GetSymbol() in self.ligand_atom_type:
                entry = [int(atom.GetIdx())]
                entry.append(atom.GetSymbol())
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                lines.append(entry)

        df = pd.DataFrame(lines)
        df.columns = ["ATOM_INDEX", "ATOM_ELEMENT", "X", "Y", "Z"]

        if len(set(df["ATOM_ELEMENT"]) - set(self.ligand_atom_type)) > 0:
            print("WARNING: Ligand contains unsupported atom types."
                  "Only supported atom-type pairs are counted.")
        return(df)

    def pdb_to_df(self, pdb_file):

        ppdb = PandasPdb()
        ppdb = ppdb.read_pdb(pdb_file)
        ppdb_all_df = ppdb.df['ATOM']
        ppdb_df = ppdb_all_df[ppdb_all_df['element_symbol'].isin(
            self.protein_atom_type)]

        atom_index = ppdb_df['atom_number']
        atom_element = ppdb_df['element_symbol']
        x, y, z = ppdb_df['x_coord'], ppdb_df['y_coord'], ppdb_df['z_coord']
        df = pd.DataFrame.from_dict({'ATOM_INDEX': atom_index,
                                     'ATOM_ELEMENT': atom_element,
                                     'X': x, 'Y': y, 'Z': z})

        return df

    def get_features(self):

        protein_len = len(self.protein_atom_type)
        ligand_len = len(self.ligand_atom_type)

        num_stat_measure = self.num_stat_measure
        num_features = protein_len*ligand_len*num_stat_measure
        final_features_mat = np.zeros([1, num_features])

        folder = str(self.pdbid)

        p_file_path = f'{self.path}/{folder}/{folder}_protein.pdb'

        l_file_path = f'{self.path}/{folder}/{folder}_ligand.sdf'

        protein_df = self.pdb_to_df(p_file_path)

        ligand_df = self.sdf_to_df(l_file_path)

        # Main loop
        pair_wise_features = np.zeros(
            (ligand_len*protein_len, num_stat_measure))

        atomic_pair_count = 0

        for i, l_atom_type in enumerate(self.ligand_atom_type):

            ligand_atom_type = l_atom_type
            ligand_vdW = self.atom_type_radii[l_atom_type]
            ligand_df_slice = ligand_df[ligand_df['ATOM_ELEMENT']
                                        == ligand_atom_type]

            l_xyz = ligand_df_slice[['X', 'Y', 'Z']].values

            for j, p_atom_type in enumerate(self.protein_atom_type):

                protein_atom_type = p_atom_type
                protein_vdW = self.atom_type_radii[p_atom_type]
                protein_df_slice = protein_df[protein_df['ATOM_ELEMENT']
                                              == protein_atom_type]

                p_xyz = protein_df_slice[['X', 'Y', 'Z']].values

                dmat = cdist(l_xyz, p_xyz, metric='euclidean')

                index = np.where(dmat < self.cutoff)
                l_index = np.unique(index[0])
                p_index = np.unique(index[1])
                ei_xyz = np.concatenate((l_xyz, p_xyz[p_index, :]))

                if np.size(ei_xyz) != 0:
                    bf = 2.0
                    x_left = min(ei_xyz[:, 0])-bf
                    x_right = max(ei_xyz[:, 0])+bf
                    y_left = min(ei_xyz[:, 1])-bf
                    y_right = max(ei_xyz[:, 1])+bf
                    z_left = min(ei_xyz[:, 2])-bf
                    z_right = max(ei_xyz[:, 2])+bf

                    h = self.mesh_size

                    x = np.arange(x_left, x_right, h)
                    y = np.arange(y_left, y_right, h)
                    z = np.arange(z_left, z_right, h)

                    nx, ny, nz = len(x), len(y), len(z)

                    eid_class = ElementInteractiveDensity(kernel_type=self.kernel_type,
                                                          kernel_tau=self.kernel_tau,
                                                          kernel_power=self.kernel_power,
                                                          ligand_vdW=ligand_vdW,
                                                          protein_vdW=protein_vdW
                                                          )

                    rho = eid_class.main(nx, ny, nz, x, y, z,
                                         ei_xyz)

                    rho_max = max(rho.flatten())

                    if rho_max != 0:
                        rho_bar = rho/rho_max
                    else:
                        rho_bar = rho

                    isovalue_features = np.zeros(len(self.isovalue_list))

                    for isovalue_count, isovalue in enumerate(self.isovalue_list):

                        f = rho_bar - isovalue

                        N_x, N_y, N_z = nmb.normal_vector_components(
                            nx, ny, nz, h, f)
                        eisa = nmb.surface_area(
                            nx, ny, nz, N_x, N_y, N_z, x, y, z, h, f, isovalue)

                        isovalue_features[isovalue_count] = eisa

                    isovalue_stat_measure = [np.sum(isovalue_features),
                                             np.mean(isovalue_features),
                                             np.median(isovalue_features),
                                             np.std(isovalue_features),
                                             max(isovalue_features),
                                             min(isovalue_features)]

                    pair_wise_features[atomic_pair_count,
                                       :] = isovalue_stat_measure

                atomic_pair_count += 1

        final_features_mat = pair_wise_features.reshape(1, num_features)

        return final_features_mat
