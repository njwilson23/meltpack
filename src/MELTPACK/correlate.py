from __future__ import division
from scipy import signal

def _normalize_chip(chip):
    return (chip-chip.mean())/chip.std()

def _autocorrelate(search_chip, ref_chip):
    return signal.fftconvolve(search_size, ref_chip[::-1,::-1], mode="same")

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
    return findoffset(search_chip.shape, ref_chip.shape, findpeak(c))

