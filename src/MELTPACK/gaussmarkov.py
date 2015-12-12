# -*- coding: utf-8 -*-
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

def _uncertainty(Ryy, Rxy, Rxx_inv):
    """ Return the diagonal of the prediction uncertainty matrix.
    (DISEP Eqn 2.398) """
    P = Ryy - Rxy*Rxx_inv*Rxy.T
    return P.diagonal()

# def _error_variance(model, Rxy, Rxx_inv):
#     """ Compute the error variance (Brink version - slow) """
#     eps = model(0) + np.ones(Rxy.shape[1])
#     for i in range(len(eps)):
#         eps[i] -= (Rxy[:,i] * Rxy[:,i].T * Rxx_inv).sum()
#     return eps

def predict(model, Xi, X, Y, eps0=1e-1, maxdist=1e3, compute_uncertainty=False):
    """ Return the Gauss-Markov minimum variance estimate for points *Xi* given
    data *Y* observed at *X*.
    (DISEP Eqn 2.397)

    Arguments:
    ----------
    model: function, returns the isotropic structure function given a distance
    Xi: np.ndarray, (n x 2)
    X: np.ndarray, (n x 2)
    Y: np.ndarray, (n)
    eps0: zero lag variance, or measurement error
    compute_uncertainty: boolean, optional

    Returns:
    --------
    (np.ndarray, np.ndarray)
    Predictions and prediction uncertainty (variance)

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

    alpha = Rxx_inv*Rxy
    # all matrices are CSC
    Yi = alpha.T*Yd + Ym
    if compute_uncertainty:
        Ryy = _model_covariance_matrix(model, kdxi, kdxi, maxdist=maxdist)
        epsi = _uncertainty(Ryy, Rxy, Rxx_inv)
    else:
        epsi = np.nan*np.empty_like(Yi)
    return Yi, epsi

def _subset_data(X, Y, n):
    idx = np.random.random_integers(0, len(Y)-1, n)
    x_ = X[idx,:]
    y_ = Y[idx]
    return x_, y_

def data_covariance(X, Y, n=500, maxdist=1e3):
    """ Estimate a structure function for data *Y* at positions *X*.
    """
    x_, y_ = _subset_data(X, Y, n)
    kd = KDTree(x_)
    d = kd.sparse_distance_matrix(kd, maxdist)
    cov = np.dot(np.atleast_2d(y_).T, np.atleast_2d(y_))
    return np.asarray(d.todense()).ravel(), cov.ravel()

def data_variogram(X, Y, n=500, maxdist=1e3):
    """ Estimate a structure function for data *Y* at positions *X*.
    The variogram is defined as

        2γ(h) = 1/N(h) * Σ(z(x)-z(x+h))²
    """
    x_, y_ = _subset_data(X, Y, n)

    # Compute differences between data
    E = np.atleast_2d(np.ones(len(y_)))
    G = (np.dot(E.T, np.atleast_2d(y_)) - np.dot(np.atleast_2d(y_).T, E))**2

    # Compute pair-wise distances
    kd = KDTree(x_)
    d = kd.sparse_distance_matrix(kd, maxdist)
    return np.asarray(d.todense()).ravel(), np.abs(G.ravel())

