""" Implementation of Gauss-Markov estimators """

import numpy as np
from scipy.spatial import cKDTree
from scipy import linalg, sparse
import scipy.sparse.linalg as splinalg

def _model_covariance_matrix(model, kd1, kd2, maxdist=1e3):
    """ Build a covariance matrix given two KD-trees and a structure function """
    D = kd1.sparse_distance_matrix(kd2, maxdist).tocsc()
    D.data = model(D.data)
    return D

def _error_variance(model, Rxy, Rxx_inv):
    """ Compute the error variance """
    eps = model(0) + np.ones(Rxy.shape[1])
    for i in range(len(eps)):
        eps[i] -= (Rxy[:,i] * Rxy[:,i].T * Rxx_inv).sum()
    return eps

def predict(model, Xi, X, Y, eps0=1e-1, maxdist=1e3, compute_error_variance=False):
    """ Return the Gauss-Markov minimum variance estimate for points *Xi* given
    data *Y* observed at *X*.

    Arguments:
    ----------
    model: function, returns the isotropic structure function given a distance
    Xi: np.ndarray, (n x 2)
    X: np.ndarray, (n x 2)
    Y: np.ndarray, (n)
    eps0: zero lag variance, or measurement error

    Returns:
    --------
    (np.ndarray, np.ndarray)
    Predictions and predicted variance

    Notes:
    ------
    Oceanographers call this optimal interpolation, or objective analysis.
    Geologists call this simple kriging.
    """
    # Matrix inversion and multiplication are somewhat slow
    # Error variance is extremely slow with the current algorithm
    Ym = Y.mean()
    Yd = Y-Ym
    
    kdx = cKDTree(X)
    kdxi = cKDTree(Xi)
    
    Rxx = _model_covariance_matrix(model, kdx, kdx, maxdist=maxdist) \
            + sparse.diags(eps0*np.ones_like(Y), 0)
    Rxy = _model_covariance_matrix(model, kdx, kdxi, maxdist=maxdist)

    Rxx_inv = splinalg.inv(Rxx)

    α = Rxx_inv*Rxy
    # all matrices are CSC at this point
    Yi = α.T*Yd + Ym
    if compute_error_variance:
        ϵi = _error_variance(model, Rxy, Rxx_inv)
    else:
        ϵi = np.nan*np.empty_like(Yi)
    return Yi, ϵi

