""" Functions for computing time derivatives on grids """

import karta
import numpy as np

def lagrangian_dhdt(g1, g2, u, v, timespan):
    """ Compute the Lagrangian time derivative of change between two grids.

    Arguments:
    ----------
    g1: RegularGrid
    g2: RegularGrid
    u: RegularGrid, horizontal velocity component
    v: RegularGrid, vertical velocity component
    timespan: datetime.timedelta, temporal baseline

    Returns:
    --------
    ndarray

        [x_0, x_1 .....
         y_0, y_1 .....
         dx_0, dx_1 ...
         dy_0, dy_1 ...
         z1_0, z1_1 ...
         z2_0, z2_1 ...]
    """

    # create a grid of points that are valid data in g1
    _ind = np.indices(g1.size)
    ind = _ind[:,~np.isnan(g1.values)]
    t1 = g1.transform
    x = t1[0] + t1[2]*ind[1]
    y = t1[1] + t1[3]*ind[0]

    # sample velocity at every point on the grid (one-sided eulerian integration)
    u = uvel.sample(x, y)
    v = vvel.sample(x, y)

    # compute displacements
    dx = u*timespan.days/365.0 
    dy = v*timespan.days/365.0 

    # sample g1 at original locations
    z1 = g1.sample(x, y)
    # sample g2 at displaced locations
    z2 = g2.sample(x+dx, y+dy)

    # limit results to where data exists in both grids
    m = np.isnan(z1) | np.isnan(z2)
    return np.vstack([x[~m], y[~m], dx[~m], dy[~m], z1[~m], z2[~m]])

