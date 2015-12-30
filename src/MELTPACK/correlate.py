from __future__ import division
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
import numpy as np
from scipy import signal

def _normalize_chip(chip):
    return (chip-chip.mean())/chip.std()

def _autocorrelate(search_chip, ref_chip):
    return signal.fftconvolve(search_chip, ref_chip[::-1,::-1], mode="same")

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

def correlate_chips(search_chip, ref_chip):
    """ Given two images, compute offsets using normalized cross-correlation.
    This is intended to be useful for image co-registration for feature
    tracking.
    """
    c = _autocorrelate(_normalize_chip(search_chip), _normalize_chip(ref_chip))
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

def correlate_scenes(scene1, scene2,
        search_chip_size=(256, 256),
        ref_chip_size=(32, 32),
        resolution=(16, 16)):

    dx, dy = scene2.transform[2:4]

    points = []
    displs = []
    strengths = []

    nprocs = cpu_count()
    overlap = (ref_chip_size[0]-resolution[0], ref_chip_size[1]-resolution[1])

    for chunk1, chunk2 in zip(scene1.aschunks(search_chip_size),
                              scene2.aschunks(search_chip_size)):

        bbox1 = chunk1.bbox
        bbox2 = chunk2.bbox
        search_center = (0.5 * (bbox1[0] + bbox1[2]),
                         0.5 * (bbox1[1] + bbox1[3]))

        futures = []
        with ThreadPoolExecutor(nprocs) as executor:
            for ref_chip in chunk2.aschunks(ref_chip_size, overlap):
                fut = executor.submit(_do_correlation, chunk1, ref_chip,
                                      search_center)
                futures.append(fut)

            for fut in as_completed(futures):
                ref_center, displ, strength = fut.result()

                points.append(ref_center)
                displs.append(displ)
                strengths.append(strength)

    return points, displs, strengths
