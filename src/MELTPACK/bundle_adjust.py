"""
Functions to compute least-squares corrections to smoothly mosaic DEMs
"""

import itertools
import numpy as np
import karta

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

def count_overlapping_pixels(grid1, grid2):
    """ Returns the number of non-NaN pixels that overlap between two grids with
    the same size and geotransform.
    """
    return np.sum(~np.isnan(grid1.values) & ~np.isnan(grid2.values))

def compute_vertical_corrections(grid_fnms, min_pixel_overlap=1000, weighting_func=None):
    """ Perform pairwaise comparisons of a list of DEMs and return a dictionary
    of {filename -> vertical correction} that minimizes the misfit between
    overlapping DEMs, subject to weights from *weighting_func(fnm1, fnm2)*.

    By default, weights are uniform.
    """
    if weighting_func is None:
        weighting_func = lambda a,b: 1.0

    c, C = _comparison_matrix(len(grid_fnms))

    # Piecewise comparison
    CD = []
    #S = []
    for i, j in c:
        dem0 = karta.read_gtiff(grid_fnms[i])
        dem1 = karta.read_gtiff(grid_fnms[j])
        
        n = utilities.count_shared_pixels(dem0, dem1)
        if n >= min_pixel_overlap:
            CD.append(np.mean(utilities.difference_shared_pixels(dem0, dem1)))
            #S.append(np.nanstd(dem0.values-dem1.values))
        else:
            CD.append(np.nan)
            #S.append(np.nan)
        del dem0, dem1

    # Augment C by removing rows where there is no appreciable overlap
    Ca = C[~np.isnan(CD),:]             # Drop non-ovelapping relations
    CaD = np.array(CD)[~np.isnan(CD)]

    corrected_grids = [fnm for tf, fnm in zip(~np.all(Ca==0, axis=0), grid_fnms) if tf]
    Ca = Ca[:,~np.all(Ca==0, axis=0)]   # Drop DEMs not involved in remaining relations

    # Compute weights
    Wdiag = np.zeros(len(Ca))
    for i,row in enumerate(Ca):
        fnm1 = corrected_grids[np.argmin(row)]
        fnm2 = corrected_grids[np.argmax(row)]
        Wdiag[i] = weighting_func(fnm1, fnm2)

    W = np.diag(Wdiag)
    Winv = np.linalg.inv(W)

    # Solve the least-squares problem
    dz = np.linalg.solve(np.dot(np.dot(Ca.T, Winv), Ca),
                         np.dot(np.dot(Ca.T, Winv), -CaD))
    corrections = {fnm: _dz for fnm, _dz in zip(corrected_grids, dz)}
    return corrections
