from __future__ import division
from concurrent.futures import (ThreadPoolExecutor, as_completed, wait,
                                FIRST_COMPLETED, ALL_COMPLETED)
from multiprocessing import cpu_count
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
    i, j = findpeak(c)
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
        refsize=(32, 32), resolution=(128, 128), nprocs=None):

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

    #overlap = ((searchsize[0]-resolution[0]), (searchsize[1]-resolution[1]))
    #pad = ((searchsize[0]-refsize[0])//2, (searchsize[1]-refsize[1])//2)
    vel_bbox = uguess.bbox
    hx = searchsize[0]//2
    hy = searchsize[1]//2

    with ThreadPoolExecutor(nprocs) as executor:

        futures = []
        for refchunk in scene1c.aschunks(refsize, (resolution[0]-refsize[0], resolution[1]-refsize[1]), copy=False):

            if len(futures) == 5000:
                for fut in as_completed(futures):
                    ref_center, displ, strength = fut.result()

                    points.append(ref_center)
                    displs.append(displ)
                    strengths.append(strength)
                futures = []

            min_valid = refsize[0]*refsize[1]
            if (refchunk.size == refsize) and not (np.any(np.isnan(refchunk.values))):
                refimage = refchunk.values
                bbox = refchunk.bbox
                xrefcenter = 0.5*(bbox[0]+bbox[2])
                yrefcenter = 0.5*(bbox[1]+bbox[3])

                if ((xrefcenter > vel_bbox[0]+hx) and (xrefcenter < vel_bbox[2]-hx) and
                    (yrefcenter > vel_bbox[1]+hy) and (yrefcenter < vel_bbox[3]-hy)):

                    # Choose a search chunk based on velocity guess
                    refcenter = (xrefcenter, yrefcenter)
                    refcenter_pt = Point(refcenter, crs=refchunk.crs)
                    u = uguess.sample(refcenter_pt)[0]
                    v = vguess.sample(refcenter_pt)[0]
                    xdispl_expected = u*dt
                    ydispl_expected = v*dt
                    try:
                        i, j = scene2c.get_indices(xrefcenter+xdispl_expected,
                                                   yrefcenter+ydispl_expected)
                    except ValueError:
                        continue
                    searchimage = scene2c.values[max(0, i-hy):min(ny-1, i+hy),
                                                 max(0, j-hx):min(nx-1, j+hx)]

                    # Submit job
                    if not np.any(np.isnan(searchimage)):
                        fut = executor.submit(_do_correlation, searchimage, refimage,
                                              refcenter, xdispl_expected, ydispl_expected,
                                              dx, dy)
                        futures.append(fut)

        # for refchunkbig, searchchunk in zip(scene1c.aschunks(searchsize, overlap, copy=False),
        #                                     scene2c.aschunks(searchsize, overlap, copy=False)):

        #     if searchchunk.size == (searchsize):

        #         if len(futures) == 2000:
        #             for fut in as_completed(futures):
        #                 ref_center, displ, strength = fut.result()

        #                 points.append(ref_center)
        #                 displs.append(displ)
        #                 strengths.append(strength)
        #             futures = []

        #         searchimage = searchchunk.values
        #         refimage = refchunkbig.values[pad[0]:pad[0]+refsize[0], pad[1]:pad[1]+refsize[1]]

        #         bbox = refchunkbig.bbox
        #         refcenter = (bbox[0]+dx*(pad[0]+0.5*refsize[0]),
        #                      bbox[1]+dy*(pad[1]+0.5*refsize[1]))

        #         fut = executor.submit(_do_correlation, searchimage, refimage,
        #                               refcenter, dx, dy)
        #         futures.append(fut)

        for fut in as_completed(futures):
            ref_center, displ, strength = fut.result()

            points.append(ref_center)
            displs.append(displ)
            strengths.append(strength)

    return np.array(points), np.array(displs), np.array(strengths)
