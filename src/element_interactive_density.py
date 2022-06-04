#!/usr/bin/env python

"""
Introduction:
    EID: element interactive density

Author:
    Masud Rana (masud.rana@uky.edu)
    
Date last modified:
    May 26, 2022

"""

import numpy as np
from scipy.spatial.distance import cdist

import numba_utils_methods as nmb


class ElementInteractiveDensity:
    def __init__(self, kernel_type='exponential', kernel_tau=1.0,
                 kernel_power=2.0, ligand_vdW=1.7, protein_vdW=1.5):

        self.kernel_type = kernel_type
        self.kernel_tau = kernel_tau
        self.kernel_power = kernel_power

        self.ligand_vdW = ligand_vdW
        self.protein_vdW = protein_vdW

        self.kernel_function = self.build_kernel_function(kernel_type)

    def build_kernel_function(self, kernel_type):
        if kernel_type[0] in ['E', 'e']:
            return self.exponential_function
        elif kernel_type[0] in ['L', 'l']:
            return self.lorentz_function

    def exponential_function(self, atomic_distance):

        eta = self.kernel_tau * (self.ligand_vdW + self.protein_vdW)

        phi = np.exp(-(atomic_distance/eta) ** self.kernel_power)

        return np.round(phi, 5)

    def lorentz_function(self, atomic_distance):

        eta = self.kernel_tau * (self.ligand_vdW + self.protein_vdW)

        phi = 1 / (1 + atomic_distance/eta) ** self.kernel_power

        return np.round(phi, 5)

    def atomic_density(self, grid_point_x, grid_point_y, grid_point_z,
                       interactive_atoms_coordinates):

        cloudpoint = np.array([[grid_point_x, grid_point_y, grid_point_z]])

        #d = cdist(cloudpoint, interactive_atoms_coordinates, metric='euclidean')
        d = nmb.distance_matrix(cloudpoint, interactive_atoms_coordinates)

        phi = self.kernel_function(d)

        return phi.sum()

    def main(self, nx, ny, nz, x, y, z, ei_xyz):

        rtn = np.zeros([nx, ny, nz])

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    rtn[i, j, k] = self.atomic_density(
                        x[i], y[j], z[k], ei_xyz)

        return rtn

