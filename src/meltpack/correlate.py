from __future__ import division
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from math import log
import numpy as np
from scipy import signal

from . import utilities

def _normalize_chip(chip):
    s = chip.std()
    if s == 0.0:
        return np.zeros_like(chip)
    else:
        return (chip-chip.mean())/s

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

def findpeak_subpixel(c):
    """ Version of findpeak with subpixel accuracy, using two 1D gaussians

    Based on formulae in Debella-Gilo and Kaab (2011), "Sub-pixel precision
    image matching for measuring surface displacements on mass movements using
    normalized cross-correlation." """
    size = c.shape
    idx = np.argmax(c)
    i = idx//size[1]
    j = idx-i*size[1]
    if (i == 0) or (i == size[0]-1) or (j == 0) or (j == size[1]-1):
        return i, j
    else:
        cmin = c.min() - 1e-8
        cij = log(c[i,j]-cmin)
        try:
            dx = (log(c[i,j-1]-cmin) - log(c[i,j+1]-cmin)) / \
                    (2*log(c[i,j+1]-cmin) - 4*cij + 2*log(c[i,j-1]-cmin))
        except ZeroDivisionError:
            dx = 0.0
        try:
            dy = (log(c[i-1,j]-cmin) - log(c[i+1,j]-cmin)) / \
                    (2*log(c[i+1,j]-cmin) - 4*cij + 2*log(c[i-1,j]-cmin))
        except ZeroDivisionError:
            dy = 0.0
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

    Returns pixel offsets (tuple) and correlation strength relative to
    correlation standard deviation (float).
    """
    c = _autocorrelate(_normalize_chip(search_chip), _normalize_chip(ref_chip), mode=mode)
    i, j = findpeak_subpixel(c)
    #i, j = findpeak(c)
    cstd = c.std()
    cmean = c.mean()
    if cstd == 0.0:
        cstd = 1e9
    return findoffset(search_chip.shape, (i, j)), \
           (c[int(round(i)),int(round(j))]-cmean)/cstd

def _do_correlation(searchimage, refimage, refcenter, ox, oy, dx, dy):
    """
    searchimage and refimage are numpy arrays
    refcenter is a tuple indicating the physical center of refimage
    ox and oy are offsets applied to the search image in physical units, which are added to the final displacements
    dx and dy are floating point grid spacings
    """
    if np.any(np.isnan(searchimage)) or np.any(np.isnan(refimage)):
        return refcenter, (np.nan, np.nan), np.nan
    (x, y), strength = correlate_chips(searchimage, refimage, mode="same")
    x_displ = x*dx + ox
    y_displ = y*dy + oy
    return refcenter, (x_displ, y_displ), strength

def correlate_scenes(scene1, scene2, uguess, vguess, dt, searchsize=(128, 128),
        refsize=(32, 32), resolution=(50.0, 50.0), nprocs=None):
    """ Compute apparent offsets between two scenes at grid points.

    scene1 : karta.RegularGrid, earlier scene with features to match

    scene2 : karta.RegularGrid, later scene with features to match

    uguess : karta.RegularGrid, grid of expected horizontal displacement rate

    vguess : karta.RegularGrid, grid of expected vertical displacement rate

    dt : float, time offset between scenes in the same units as uguess/vguess

    searchsize : tuple(int, int), size of extracted search region. A larger
        search region permits greater deviation from the absolute velocity
        guess, but becomes less robust to deformation gradients.
        (default (128, 128))

    refsize : tuple(int, int), size of reference region to match within search
        region. Typically, refsize should be at least 1/4 the size of
        searchsize.
        (default (32, 32))

    resolution : tuple(float, float), resolution of sampling grid in projected
        units.

    nprocs : number of worker threads to launch
    """
    bboxc = utilities.overlap_bbox(scene1.data_bbox, scene2.data_bbox)
    scene1c = scene1.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    scene2c = scene2.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    dx, dy = scene2c.resolution
    ny, nx = scene2c.size

    points = []
    displs = []
    strengths = []

    if nprocs is None:
        nprocs = cpu_count()

    rhx = refsize[0]//2
    rhy = refsize[1]//2
    shx = searchsize[0]//2
    shy = searchsize[1]//2

    # compute reference chip centers
    xmin1, xmax1, ymin1, ymax1 = scene1c.extent
    xmin2, xmax2, ymin2, ymax2 = scene2c.extent
    xmin = max(xmin1, xmin2)
    xmax = min(xmax1, xmax2)
    ymin = max(ymin1, ymin2)
    ymax = min(ymax1, ymax2)
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
    vdx, vdy = uguess.resolution
    xmin, xmax, ymin, ymax = uguess.extent
    mask = (Xref<xmin+vdx) | (Xref>xmax-vdx) | (Yref<ymin+vdy) | (Yref>ymax-vdy)
    Xref = Xref[~mask]
    Yref = Yref[~mask]

    # get indices for ref centers
    Iref, Jref = scene1c.get_indices(Xref, Yref)

    # compute the *actual* gridded reference chip centers
    T = scene1c.transform
    Xref = T[0] + Jref*T[2] + Iref*T[4]
    Yref = T[1] + Iref*T[3] + Jref*T[5]

    # compute velocity guesses
    uref = uguess.sample(Xref, Yref)
    vref = vguess.sample(Xref, Yref)

    # compute expected offsets
    offx = np.round(uref*dt/dx).astype(np.int16)
    offy = np.round(vref*dt/dy).astype(np.int16)

    # extract chips
    searchchips = []
    refchips = []
    val1 = scene1c.values
    val2 = scene2c.values

    for ir, jr, ox, oy in zip(Iref, Jref, offx, offy):

        refchip = val1[max(0, ir-rhy):min(ny-1, ir+rhy),
                       max(0, jr-rhx):min(nx-1, jr+rhx)]
        searchchip = val2[max(0, ir+oy-shy):min(ny-1, ir+oy+shy),
                          max(0, jr+ox-shx):min(nx-1, jr+ox+shx)]
        refchips.append(refchip)
        searchchips.append(searchchip)

    # extract chips and farm out to threadpool
    with ThreadPoolExecutor(nprocs) as executor:

        futures = []
        for xr, yr, schip, rchip, ox, oy in zip(Xref, Yref, searchchips, refchips, offx, offy):

            if (schip.shape == searchsize) and (rchip.shape == refsize):
                fut = executor.submit(_do_correlation, schip, rchip,
                                      (xr, yr), ox*dx, oy*dy, dx, dy)
                futures.append(fut)

        for fut in as_completed(futures):
            ref_center, displ, strength = fut.result()

            points.append(ref_center)
            displs.append(displ)
            strengths.append(strength)

    return np.array(points), np.array(displs), np.array(strengths)

def correlate_scenes_at_points(scene1, scene2, uguess, vguess, dt, corrpoints,
        searchsize=(128, 128), refsize=(32, 32), nprocs=None):
    """ Compute apparent offsets between two scenes at grid points.

    scene1 : karta.RegularGrid, earlier scene with features to match

    scene2 : karta.RegularGrid, later scene with features to match

    uguess : karta.RegularGrid, grid of expected horizontal displacement rate

    vguess : karta.RegularGrid, grid of expected vertical displacement rate

    dt : float, time offset between scenes in the same units as uguess/vguess

    corrpoints : []karta.Point, list of points at which to compute a correlation
        vector

    searchsize : tuple(int, int), size of extracted search region. A larger
        search region permits greater deviation from the absolute velocity
        guess, but becomes less robust to deformation gradients.
        (default (128, 128))

    refsize : tuple(int, int), size of reference region to match within search
        region. Typically, refsize should be at least 1/4 the size of
        searchsize.
        (default (32, 32))

    nprocs : number of worker threads to launch
    """
    bboxc = utilities.overlap_bbox(scene1.data_bbox, scene2.data_bbox)
    scene1c = scene1.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    scene2c = scene2.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    dx, dy = scene2c.resolution
    ny, nx = scene2c.size

    points = []
    displs = []
    strengths = []

    if nprocs is None:
        nprocs = cpu_count()

    rhx = refsize[0]//2
    rhy = refsize[1]//2
    shx = searchsize[0]//2
    shy = searchsize[1]//2

    # compute reference chip centers
    Xref = np.array([pt.x for pt in corrpoints])
    Yref = np.array([pt.y for pt in corrpoints])

    # filter out nan locations
    v1 = scene1c.sample(Xref, Yref)
    v2 = scene2c.sample(Xref, Yref)
    mask = np.isnan(v1) | np.isnan(v2)
    Xref = Xref[~mask]
    Yref = Yref[~mask]

    # filter out locations beyond the velocity grid
    vdx, vdy = uguess.resolution
    xmin, xmax, ymin, ymax = uguess.extent
    mask = (Xref<xmin+vdx) | (Xref>xmax-vdx) | (Yref<ymin+vdy) | (Yref>ymax-vdy)
    Xref = Xref[~mask]
    Yref = Yref[~mask]

    # get indices for ref centers
    Iref, Jref = scene1c.get_indices(Xref, Yref)

    # compute the *actual* gridded reference chip centers
    T = scene1c.transform
    Xref = T[0] + Jref*T[2] + Iref*T[4]
    Yref = T[1] + Iref*T[3] + Jref*T[5]

    # compute velocity guesses
    uref = uguess.sample(Xref, Yref)
    vref = vguess.sample(Xref, Yref)

    # compute expected offsets
    offx = np.round(uref*dt/dx).astype(np.int16)
    offy = np.round(vref*dt/dy).astype(np.int16)

    # extract chips
    searchchips = []
    refchips = []
    val1 = scene1c.values
    val2 = scene2c.values

    for ir, jr, ox, oy in zip(Iref, Jref, offx, offy):

        refchip = val1[max(0, ir-rhy):min(ny-1, ir+rhy),
                       max(0, jr-rhx):min(nx-1, jr+rhx)]
        searchchip = val2[max(0, ir+oy-shy):min(ny-1, ir+oy+shy),
                          max(0, jr+ox-shx):min(nx-1, jr+ox+shx)]
        refchips.append(refchip)
        searchchips.append(searchchip)

    # farm out to threadpool
    with ThreadPoolExecutor(nprocs) as executor:

        futures = []
        for xr, yr, schip, rchip, ox, oy in zip(Xref, Yref, searchchips, refchips, offx, offy):

            if (schip.shape == searchsize) and (rchip.shape == refsize):
                fut = executor.submit(_do_correlation, schip, rchip,
                                      (xr, yr), ox*dx, oy*dy, dx, dy)
                futures.append(fut)

        for fut in as_completed(futures):
            ref_center, displ, strength = fut.result()

            points.append(ref_center)
            displs.append(displ)
            strengths.append(strength)

    return np.array(displs), np.array(strengths), np.dstack(refchips), np.dstack(searchchips)
