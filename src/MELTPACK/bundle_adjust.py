"""
Functions to compute least-squares corrections to smoothly mosaic DEMs
"""

import itertools
import numpy as np
import karta

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

def compute_vertical_corrections(grid_fnms, min_pixel_overlap=1000):
    """ Perform pairwaise comparisons of a list of DEMs and return a dictionary
    of {filename -> vertical correction} that minimizes the misfit between
    overlapping DEMs.
    """
    c, C = _comparison_matrix(len(grid_fnms))

    # Piecewise comparison
    CD = []
    S = []
    for i, j in c:
        dem0 = karta.read_gtiff(os.path.join(basedir, dem_list[i]))
        dem1 = karta.read_gtiff(os.path.join(basedir, dem_list[j]))
        
        n = count_overlapping_pixels(dem0, dem1)
        if n >= min_pixel_overlap:
            CD.append(np.nanmean(dem0.values-dem1.values))
            S.append(np.nanstd(dem0.values-dem1.values))
        else:
            CD.append(np.nan)
            S.append(np.nan)
        del dem0, dem1

    # Augment C by removing rows where there is no appreciable overlap
    Ca = C[~np.isnan(CD),:]             # Drop non-ovelapping relations
    CaD = np.array(CD)[~np.isnan(CD)]

    corrected_grids = [fnm for tf, fnm in zip(~np.all(Ca==0, axis=0), grid_fnms) if tf]
    Ca = Ca[:,~np.all(Ca==0, axis=0)]   # Drop DEMs not involved in remaining relations

    # Solve the least-squares problem
    dz = np.linalg.solve(np.dot(Ca.T, Ca),
                         np.dot(Ca.T, -CaD))
    corrections = {fnm: _dz for fnm, _dz in zip(corrected_dems, dz)}
    return corrections
