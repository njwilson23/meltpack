""" From a set of point data, compute bounding alpha shapes. """

import math
import karta
from shapely.geometry import MultiLineString
from shapely.ops import polygonize, cascaded_union
import scipy.spatial
import numpy as np

def alpha_shape(x, y, alpha):
    coords = np.c_[x, y]
    tri = scipy.spatial.Delaunay(coords)
    
    # edges = []
    edge_coords = []
    for ia, ib, ic in tri.vertices:
        ab = math.sqrt((coords[ia,0]-coords[ib,0])**2 + (coords[ia,1]-coords[ib,1])**2)
        bc = math.sqrt((coords[ib,0]-coords[ic,0])**2 + (coords[ib,1]-coords[ic,1])**2)
        ac = math.sqrt((coords[ia,0]-coords[ic,0])**2 + (coords[ia,1]-coords[ic,1])**2)
        
        semiperim = 0.5*(ab+bc+ac)
        area = math.sqrt(semiperim*(semiperim-ab)*(semiperim-bc)*(semiperim-ac))
        radius = 0.25*ab*bc*ac/area
        
        if radius < 1.0/alpha:
            # edges.append((tuple(coords[ia]), tuple(coords[ib])))
            # edges.append((tuple(coords[ib]), tuple(coords[ic])))
            # edges.append((tuple(coords[ia]), tuple(coords[ic])))
            edge_coords.append(coords[[ia, ib]])
            edge_coords.append(coords[[ib, ic]])
            edge_coords.append(coords[[ia, ic]])
            
    # edges = set(edges)
    m = MultiLineString(edge_coords)
    return cascaded_union(list(polygonize(m))), edge_coords

def restrict_grid(grid, x_obs, y_obs, alpha):
    """ Mask out grid beyond the (concave) region bounded by x_obs, y_obs
    Smaller *alpha* results in a smoother bounding shape.

    If the restriction step results in multiple disjoint shapes, currently only
    the largest is used. Extending to use all should be easy, if it's ever
    wanted.
    """
    p, _ = alpha_shape(x_obs, y_obs, alpha)
    try:
        areas = [_p.area for _p in p]
        kp = karta.vector.from_shape(p[areas.index(max(areas))])
    except TypeError:
        kp = karta.vector.from_shape(p)
    kp.crs = grid.crs
    return grid.mask_by_poly(kp)

