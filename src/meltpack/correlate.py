from __future__ import division
from concurrent.futures import (ThreadPoolExecutor, as_completed, wait,
                                FIRST_COMPLETED, ALL_COMPLETED)
from multiprocessing import cpu_count
import numpy as np
from scipy import signal

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

def _do_correlation(searchimage, refimage, refcenter, dx, dy):
    (x, y), strength = correlate_chips(searchimage, refimage, mode="same")
    x_displ = x*dx
    y_displ = y*dy
    return refcenter, (x_displ, y_displ), strength

def correlate_scenes(scene1, scene2, searchsize=(256, 256), refsize=(32, 32),
        resolution=(128, 128), nprocs=None):

    bboxc = utilities.overlap_bbox(scene1.data_bbox, scene2.data_bbox)
    scene1c = scene1.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    scene2c = scene2.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    dx, dy = scene2c.transform[2:4]

    points = []
    displs = []
    strengths = []

    if nprocs is None:
        nprocs = cpu_count()

    overlap = ((searchsize[0]-resolution[0]), (searchsize[1]-resolution[1]))
    dx, dy = scene1c.transform[2:4]
    pad = ((searchsize[0]-refsize[0])//2, (searchsize[1]-refsize[1])//2)

    with ThreadPoolExecutor(nprocs) as executor:

        futures = []
        for refchunkbig, searchchunk in zip(scene1c.aschunks(searchsize, overlap),
                                            scene2c.aschunks(searchsize, overlap)):

            if searchchunk.size == (searchsize):

                if len(futures) == 2000:
                    for fut in as_completed(futures):
                        ref_center, displ, strength = fut.result()

                        points.append(ref_center)
                        displs.append(displ)
                        strengths.append(strength)
                    futures = []

                searchimage = searchchunk.values
                refimage = refchunkbig.values[pad[0]:pad[0]+refsize[0], pad[1]:pad[1]+refsize[1]]

                bbox = refchunkbig.bbox
                refcenter = (bbox[0]+dx*(pad[0]+0.5*refsize[0]),
                             bbox[1]+dy*(pad[1]+0.5*refsize[1]))

                fut = executor.submit(_do_correlation, searchimage, refimage,
                                      refcenter, dx, dy)
                futures.append(fut)

        for fut in as_completed(futures):
            ref_center, displ, strength = fut.result()

            points.append(ref_center)
            displs.append(displ)
            strengths.append(strength)

    return points, displs, strengths
