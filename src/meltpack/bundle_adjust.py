"""
Functions to compute least-squares corrections to smoothly mosaic DEMs

Method:

    Defines a comparison matrix C (n x m) where m is the number of DEMs and n =
    binom(m, 2).

    Compute the vector matrix product CD, which is the DEM differences.

    Mask out parts of the system corresponding to non-overlapping DEMs and
    unconstrained DEMs, giving C' and CD'

    Compute diagonal matrix of weights W based on the time difference between
    DEMs. This makes only makes sense if the DEM differences are over the ice
    (and liable to change with time), and not in bedrock polygons.

    Solve the least squares system C'^T inv(W) C' dZ = C'T inv(W) CD' for the
    vertical correction vector dZ.
"""

import itertools
import numpy as np
import karta

from . import correlate
from . import utilities

def _comparison_matrix(n):
    c = list(itertools.combinations(range(n), 2))
    C = np.zeros([len(c), n])

    I = np.arange(len(c))
    J0 = np.array([i[0] for i in c])
    J1 = np.array([i[1] for i in c])
    C[I,J0] = 1.0
    C[I,J1] = -1.0
    return c, C

def compute_vertical_corrections(grids, min_pixel_overlap=100,
    weighting_func=None, polymasks=()):
    """ Perform pairwaise comparisons of a list of DEMs and return a dictionary
    of {index -> vertical correction} that minimizes the misfit between
    overlapping DEMs, subject to weights from *weighting_func(fnm1, fnm2)*.

    Parameters:
    -----------
    grids : list of filenames of RegularGrid-like instances
        Grids to compute corrections for. Grids must have the same extent and
        resolution.
    min_pixel_overlap : int, optional
        Minimum number of overlapping data pixels for computing grid
        relationships (default 100)
    weighting_func : callable
        Function that takes two of the input grids (either filenames or
        RegularGrid instances) and returns a weight. (default, uniforma
        weighting).
    polymasks : list or tuple of Polygon instances
        Polygons that define the data region to use for comparisons, for
        example to restrict grid correction to comparisons between stable
        bedrock polygons. If empty (default), all data pixels are considered
        valid for comparison.
    """
    if isinstance(grids[0], str):
        gridnames = True
    else:
        gridnames = False

    if weighting_func is None:
        weighting_func = lambda a,b: 1.0

    c, C = _comparison_matrix(len(grids))

    # Piecewise comparison (computes the matrix-vector product CD)
    CD = []
    #S = []
    for i, j in c:
        if gridnames:
            dem0 = karta.read_gtiff(grids[i])
            dem1 = karta.read_gtiff(grids[j])
            dem0.values[dem0.values<-1000000] = np.nan
            dem1.values[dem1.values<-1000000] = np.nan
        else:
            dem0 = grids[i]
            dem1 = grids[j]

        # If bedrock regions are provided, limit the overlapping pixel counts
        # the regions inside the bedrock masking polygons.
        # Iterate through all polygons, keeping track of pixels that are in at
        # least one polygon in the *mask_count* array variables.
        # Then, set all DEM pixel where *mask_count* is zero to NODATA
        if polymasks is not None and len(polymasks)!=0:
            mask_count0 = np.zeros(dem0.size, dtype=np.int16)
            mask_count1 = np.zeros(dem1.size, dtype=np.int16)
            for poly in polymasks:

                x = [a[0] for a in poly.get_vertices(dem0.crs)]
                y = [a[1] for a in poly.get_vertices(dem0.crs)]

                ny, nx = dem0.size
                msk0 = karta.raster.grid.mask_poly(x, y, nx, ny, dem0.transform)
                ny, nx = dem1.size
                msk1 = karta.raster.grid.mask_poly(x, y, nx, ny, dem1.transform)

                mask_count0[msk0] += 1
                mask_count1[msk1] += 1

            dem0.values[mask_count0==0] = dem0.nodata
            dem1.values[mask_count1==0] = dem1.nodata

        n = utilities.count_shared_pixels(dem0, dem1)
        if n >= min_pixel_overlap:
            CD.append(np.mean(utilities.difference_shared_pixels(dem0, dem1)))
        else:
            CD.append(np.nan)
        del dem0, dem1

    # Augment C by removing rows where there is no appreciable overlap
    Ca = C[~np.isnan(CD),:]             # Drop non-ovelapping relations
    CaD = np.array(CD)[~np.isnan(CD)]

    corrected_grids = [i for tf, i in zip(~np.all(Ca==0, axis=0), range(len(grids))) if tf]
    Ca = Ca[:,~np.all(Ca==0, axis=0)]   # Drop DEMs not involved in remaining relations

    # Compute weights
    Wdiag = np.zeros(len(Ca))
    for i,row in enumerate(Ca):
        grid1 = grids[corrected_grids[np.argmin(row)]]
        grid2 = grids[corrected_grids[np.argmax(row)]]
        Wdiag[i] = weighting_func(grid1, grid2)

    W = np.diag(Wdiag)
    Winv = np.linalg.inv(W)

    # Solve the weighted least-squares problem: C^T W^-1 C dz = -C^T W^-1 C D
    #dz = np.linalg.solve(np.dot(np.dot(Ca.T, Winv), Ca),
    #                     np.dot(np.dot(Ca.T, Winv), -CaD))

    A = np.dot(np.dot(Ca.T, Winv), Ca)
    RHS = np.dot(np.dot(Ca.T, Winv), -CaD)

    # Append a constraint to ensure that the mean of dz is zero
    n = A.shape[0]
    A = np.r_[np.c_[A, np.zeros(n)], np.ones([1,n+1])]
    RHS = np.r_[RHS, 0.0]

    # A^TA may be singular, so solve system with the pseudoinverse
    dz = np.dot(np.linalg.pinv(A), RHS)[:-1]
    return {i: _dz for i, _dz in zip(corrected_grids, dz)}

def compute_horizontal_corrections(grid_fnms, min_pixel_overlap=100,
    weighting_func=None, polymasks=()):
    """ Compute the least squares horizontal (dx, dy) corrections to align two
    images that share at least *min_pixel_overlap* pixels. """
    if weighting_func is None:
        weighting_func = lambda a,b: 1.0

    raise NotImplementedError()

    # c, C = _comparison_matrix(len(grid_fnms))

    # # Piecewise comparison
    # CD = []
    # for i, j in c:
    #     dem0 = karta.read_gtiff(grid_fnms[i])
    #     dem1 = karta.read_gtiff(grid_fnms[j])
    #     dem0.values[dem0.values<-1000000] = np.nan
    #     dem1.values[dem1.values<-1000000] = np.nan

    #     # If bedrock regions are provided, limit the overlapping pixel counts
    #     # the regions inside the bedrock masking polygons.
    #     # Iterate through all polygons, keeping track of pixels that are in at
    #     # least one polygon in the *mask_count* array variables.
    #     # Then, set all DEM pixel where *mask_count* is zero to NODATA
    #     if polymasks is not None and len(polymasks)!=0:
    #         mask_count0 = np.zeros(dem0.size, dtype=np.int16)
    #         mask_count1 = np.zeros(dem1.size, dtype=np.int16)
    #         for poly in polymasks:

    #             x = [a[0] for a in poly.get_vertices(dem0.crs)]
    #             y = [a[1] for a in poly.get_vertices(dem0.crs)]

    #             ny, nx = dem0.size
    #             msk0 = karta.raster.grid.mask_poly(x, y, nx, ny, dem0.transform)
    #             ny, nx = dem1.size
    #             msk1 = karta.raster.grid.mask_poly(x, y, nx, ny, dem1.transform)

    #             mask_count0[msk0] += 1
    #             mask_count1[msk0] += 1

    #         dem0.values[mask_count0==0] = dem0.nodata
    #         dem1.values[mask_count0==0] = dem1.nodata

    #     n = utilities.count_shared_pixels(dem0, dem1)
    #     if n >= min_pixel_overlap:
    #         CD.append(np.mean(utilities.difference_shared_pixels(dem0, dem1)))
    #     else:
    #         CD.append(np.nan)
    #     del dem0, dem1

    # # Augment C by removing rows where there is no appreciable overlap
    # Ca = C[~np.isnan(CD),:]             # Drop non-ovelapping relations
    # CaD = np.array(CD)[~np.isnan(CD)]

    # corrected_grids = [fnm for tf, fnm in zip(~np.all(Ca==0, axis=0), grid_fnms) if tf]
    # Ca = Ca[:,~np.all(Ca==0, axis=0)]   # Drop DEMs not involved in remaining relations

    # # Compute weights
    # Wdiag = np.zeros(len(Ca)+1)     # size +1 to leave room for constraint
    # for i,row in enumerate(Ca):
    #     fnm1 = corrected_grids[np.argmin(row)]
    #     fnm2 = corrected_grids[np.argmax(row)]
    #     Wdiag[i] = weighting_func(fnm1, fnm2)
    # Wdiag[-1] = 1.0

    # W = np.diag(Wdiag)
    # Winv = np.linalg.inv(W)

    # # Add a constraint to pin the solution in z (first DEM held fixed)
    # # otherwise there are infinite solutions dz + constant
    # Ca = np.vstack([Ca, np.r_[1, np.zeros(Ca.shape[1]-1)]])
    # CaD = np.r_[CaD, 0.0]

    # # Solve the least-squares problem
    # dz = np.linalg.solve(np.dot(np.dot(Ca.T, Winv), Ca),
    #                      np.dot(np.dot(Ca.T, Winv), -CaD))
    # return {fnm: _dz for fnm, _dz in zip(corrected_grids, dz)}

