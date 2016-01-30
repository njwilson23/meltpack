from __future__ import division
from concurrent.futures import (ThreadPoolExecutor, as_completed, wait,
                                FIRST_COMPLETED, ALL_COMPLETED)
from multiprocessing import cpu_count
from math import log
import numpy as np
from scipy import signal
from karta import Point

from . import utilities

def _normalize_chip(chip):
    s = chip.std()
    if s == 0.0:
        return np.zeros_like(chip)
    else:
        return (chip-chip.mean())/chip.std()

def _autocorrelate(search_chip, ref_chip, mode="valid"):
    return signal.fftconvolve(search_chip, ref_chip[::-1,::-1], mode=mode)

def findpeak(c):
    """ Return in the integer row, column indices of the largest value in array
    *c*. """
    size = c.shape
    idx = np.argmax(c)
    i = idx//size[1]
    j = idx-i*size[1]
    return i, j

def findpeak_subpixel(c, size):
    """ Version of findpeak with subpixel accuracy, using two 1D gaussians

    Based on formulae in Debella-Gilo and Kaab (2011), "Sub-pixel precision
    image matching for measuring surface displacements on mass movements using
    normalized cross-correlation." """
    idx = np.argmax(c)
    i = idx//size[1]
    j = idx-i*size[1]
    
    cij = log(c[i,j])
    dx = (log(c[i,j-1]) - log(c[i,j+1])) / (2*log(c[i,j+1]) - 4*cij + 2*log(c[i,j-1]))
    dy = (log(c[i-1,j]) - log(c[i+1,j])) / (2*log(c[i+1,j]) - 4*cij + 2*log(c[i-1,j]))
    return i+dy, j+dx

def findoffset(size, peak):
    """ Given an array size and peak indices, return the offset from the array
    center """
    y = peak[0] - size[0]/2
    x = peak[1] - size[1]/2
    return x, y

def correlate_chips(search_chip, ref_chip, mode="valid"):
    """ Given two images, compute offsets using normalized cross-correlation.
    This is intended to be useful for image co-registration for feature
    tracking.

    Set `mode="same"` for image co-registration. Use `mode="valid"` for feature
    tracking.
    """
    c = _autocorrelate(_normalize_chip(search_chip), _normalize_chip(ref_chip), mode=mode)
    i, j = findpeak_subpixel(c)
    return findoffset(search_chip.shape, (i, j)), c[i,j]

def _do_correlation(searchimage, refimage, refcenter, ox, oy, dx, dy):
    """
    searchimage and refimage are numpy arrays
    refcenter is a tuple indicating the physical center of refimage
    ox and oy are offsets applied to the search image in physical units, which are added to the final displacements
    dx and dy are floating point grid spacings
    """
    (x, y), strength = correlate_chips(searchimage, refimage, mode="same")
    x_displ = x*dx+ox
    y_displ = y*dy+oy
    return refcenter, (x_displ, y_displ), strength

def correlate_scenes(scene1, scene2, uguess, vguess, dt, searchsize=(128, 128),
        refsize=(32, 32), resolution=(50.0, 50.0), nprocs=None):

    bboxc = utilities.overlap_bbox(scene1.data_bbox, scene2.data_bbox)
    scene1c = scene1.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    scene2c = scene2.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    dx, dy = scene2c.transform[2:4]
    ny, nx = scene2c.size

    points = []
    displs = []
    strengths = []

    if nprocs is None:
        nprocs = cpu_count()

    vel_bbox = uguess.bbox
    rhx = refsize[0]//2
    rhy = refsize[1]//2
    shx = searchsize[0]//2
    shy = searchsize[1]//2


    # compute reference chip centers
    xmin1, xmax1, ymin1, ymax1 = scene1c.extent
    xmin2, xmax2, ymin2, ymax2 = scene2c.extent
    xmin = max(xmin1, xmin2) + dx
    xmax = min(xmax1, xmax2) - dx
    ymin = max(ymin1, ymin2) + dy
    ymax = min(ymax1, ymax2) - dy
    x = np.arange(xmin, xmax, resolution[0])
    y = np.arange(ymin, ymax, resolution[1])
    Xref, Yref = np.meshgrid(x, y)
    Xref = Xref.ravel()
    Yref = Yref.ravel()

    # filter out nan locations
    v1 = scene1c.sample(Xref, Yref)
    v2 = scene2c.sample(Xref, Yref)
    mask = np.isnan(v1) | np.isnan(v2)
    Xref = Xref[~mask]
    Yref = Yref[~mask]

    # filter out locations beyond the velocity grid
    xmin, xmax, ymin, ymax = uguess.extent
    mask = (Xref<xmin) | (Xref>xmax) | (Yref<ymin) | (Yref>ymax)
    Xref = Xref[~mask]
    Yref = Yref[~mask]

    # get indices for ref centers
    Iref, Jref = scene1c.get_indices(Xref, Yref)

    # compute velocity guesses
    uref = uguess.sample(Xref, Yref)
    vref = vguess.sample(Xref, Yref)

    # compute expected offsets
    offx = np.round(uref*dt/dx).astype(np.int16)
    offy = np.round(vref*dt/dy).astype(np.int16)

    # extract chips and farm out to threadpool
    with ThreadPoolExecutor(nprocs) as executor:

        futures = []
        for xrefcenter, yrefcenter, ir, jr, ox, oy in zip(Xref, Yref, Iref, Jref, offx, offy):

            if len(futures) == 5000:
                for fut in as_completed(futures):
                    ref_center, displ, strength = fut.result()

                    points.append(ref_center)
                    displs.append(displ)
                    strengths.append(strength)
                futures = []

            refchip = scene1c.values[max(0, ir-rhy):min(ny-1, ir+rhy),
                                     max(0, jr-rhx):min(nx-1, jr+rhx)]
            searchchip = scene2c.values[max(0, ir+oy-shy):min(ny-1, ir+oy+shy),
                                        max(0, jr+ox-shx):min(nx-1, jr+ox+shx)]

            if ((searchchip.shape == searchsize) and (refchip.shape == refsize) and
                (not np.any(np.isnan(searchchip))) and (not np.any(np.isnan(refchip)))):

                fut = executor.submit(_do_correlation, searchchip, refchip,
                                      (xrefcenter, yrefcenter), ox*dx, oy*dy, dx, dy)
                futures.append(fut)

        for fut in as_completed(futures):
            ref_center, displ, strength = fut.result()

            points.append(ref_center)
            displs.append(displ)
            strengths.append(strength)

    return np.array(points), np.array(displs), np.array(strengths)
