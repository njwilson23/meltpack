from __future__ import division
from concurrent.futures import ThreadPoolExecutor, as_completed
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

def _do_correlation(search_img, ref_img, center_coords):
    (x, y), strength = correlate_chips(search_img.values, ref_img.values)

    # compute the offsets in geographical units
    ref_bbox = ref_img.bbox
    ref_center = (0.5 * (ref_bbox[0] + ref_bbox[2]),
                  0.5 * (ref_bbox[1] + ref_bbox[3]))
    dx, dy = search_img.transform[2:4]
    xg = center_coords[0] + x*dx
    yg = center_coords[1] + y*dy
    x_displ = ref_center[0] - xg
    y_displ = ref_center[1] - yg
    return ref_center, (x_displ, y_displ), strength

def correlate_scenes(scene1, scene2, search_size=(256, 256), ref_size=(32, 32),
        resolution=(16, 16)):

    bboxc = utilities.overlap_bbox(scene1.data_bbox, scene2.data_bbox)
    scene1c = scene1.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    scene2c = scene2.clip(bboxc[0], bboxc[2], bboxc[1], bboxc[3])
    # scene1c = scene1
    # scene2c = scene2
    dx, dy = scene2c.transform[2:4]

    points = []
    displs = []
    strengths = []

    nprocs = cpu_count()
    overlap = (ref_size[0]-resolution[0], ref_size[1]-resolution[1])

    for chunk1, chunk2 in zip(scene1c.aschunks(search_size, ref_size),
                              scene2c.aschunks(search_size, ref_size)):

        bbox1 = chunk1.bbox
        bbox2 = chunk2.bbox
        search_center = (0.5 * (bbox1[0] + bbox1[2]),
                         0.5 * (bbox1[1] + bbox1[3]))

        futures = []
        with ThreadPoolExecutor(nprocs) as executor:
            for ref_chip in chunk2.aschunks(ref_size, overlap):
                fut = executor.submit(_do_correlation, chunk1, ref_chip,
                                      search_center)
                futures.append(fut)

            for fut in as_completed(futures):
                ref_center, displ, strength = fut.result()

                points.append(ref_center)
                displs.append(displ)
                strengths.append(strength)

    return points, displs, strengths
