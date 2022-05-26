"""
Introduction:
    Methods using numba python

Author:
    Masud Rana (masud.rana@uky.edu)

Date last modified:
    May 26, 2022

"""


import numpy as np

from numba import jit
from math import sqrt

# import os


@jit(nopython=True)
def distance_matrix(Lxyz, Pxyz):
    locsP = Pxyz
    locsL = Lxyz

    num_Patom = len(locsP)
    num_Latom = len(locsL)

    dmat = np.zeros((num_Latom, num_Patom))

    for i in range(num_Latom):
        for j in range(num_Patom):
            d = sqrt((locsL[i, 0]-locsP[j, 0])**2 +
                     (locsL[i, 1] - locsP[j, 1])**2 +
                     (locsL[i, 2]-locsP[j, 2])**2)
            dmat[i, j] = d

    return dmat


@jit(nopython=True)
def normal_vector_components(nx, ny, nz, h, f):
    N_x = np.zeros((nx-2, ny-2, nz-2))
    N_y = np.zeros((nx-2, ny-2, nz-2))
    N_z = np.zeros((nx-2, ny-2, nz-2))
    # epsilon = 1.0e-10
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                N_x[i-1, j-1, k-1] = (f[i+1, j, k]-f[i-1, j, k])/(2*h)
                N_y[i-1, j-1, k-1] = (f[i, j+1, k]-f[i, j-1, k])/(2*h)
                N_z[i-1, j-1, k-1] = (f[i, j, k+1]-f[i, j, k-1])/(2*h)

    return N_x, N_y, N_z


@jit(nopython=True)
def surface_area(nx, ny, nz, N_x, N_y, N_z, x, y, z, h, f, isovalue):

    area = 0.0

    # count = 0
    for i in range(nx-3):
        for j in range(ny-3):
            for k in range(nz-3):

                # x-component of the normal vector
                # at the intersection point (x_o,y_j,z_k)

                No_x = 0.0

                # y-component of the normal vector
                # at the intersection point (x_i,y_o,z_k)

                No_y = 0.0

                # z-component of the normal vector
                # at the intersection point (x_i,y_j,z_o)

                No_z = 0.0

                # Irregular grid points: points with neighbor
                # from the other side of the interface
                # Irregular grid points along x-mesh line:
                # f(x_i,y_j,z_k)*f(x_i+1,y_j,z_k)<=0
                # and f(x_i+1,y_j,z_k)*f(x_i-1,y_j,z_k)<0.

                if ( (f[i+1, j, k]*f[i, j, k]) <= 0
                        and f[i+1, j, k]*f[i-1, j, k] )< 0:

                    factor = (isovalue-f[i, j, k])/(f[i+1, j, k]-f[i, j, k])
                    No_x1 = factor*(N_x[i+1, j, k]-N_x[i, j, k]) + N_x[i, j, k]
                    No_x2 = factor*(N_y[i+1, j, k]-N_y[i, j, k]) + N_y[i, j, k]
                    No_x3 = factor*(N_z[i+1, j, k]-N_z[i, j, k]) + N_z[i, j, k]
                    norm = np.sqrt(No_x1**2+No_x2**2+No_x3**2)
                    No_x = No_x1/norm if norm != 0 else 0.0

                if ( (f[i, j+1, k]*f[i, j, k]) <= 0 
                       and f[i, j+1, k]*f[i, j-1, k] ) < 0:

                    factor = (isovalue-f[i, j, k])/(f[i, j+1, k]-f[i, j, k])
                    No_y1 = factor*(N_x[i, j+1, k]-N_x[i, j, k]) + N_x[i, j, k]
                    No_y2 = factor*(N_y[i, j+1, k]-N_y[i, j, k]) + N_y[i, j, k]
                    No_y3 = factor*(N_z[i, j+1, k]-N_z[i, j, k]) + N_z[i, j, k]
                    norm = np.sqrt(No_y1**2+No_y2**2+No_y3**2)
                    No_y = No_y2/norm if norm != 0 else 0

                if ( (f[i, j, k+1]*f[i, j, k]) <= 0 
                      and f[i, j, k+1]*f[i, j, k-1] ) < 0:

                    factor = (isovalue-f[i, j, k])/(f[i, j, k+1]-f[i, j, k])
                    No_z1 = factor*(N_x[i, j, k+1]-N_x[i, j, k]) + N_x[i, j, k]
                    No_z2 = factor*(N_y[i, j, k+1]-N_y[i, j, k]) + N_y[i, j, k]
                    No_z3 = factor*(N_z[i, j, k+1]-N_z[i, j, k]) + N_z[i, j, k]
                    norm = np.sqrt(No_z1**2+No_z2**2+No_z3**2)
                    No_z = No_z3/norm if norm != 0 else 0

                area = area + (abs(No_x) + abs(No_y) + abs(No_z))*h**2
                # count +=1

    return area
