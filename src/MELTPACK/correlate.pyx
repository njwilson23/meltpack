import cython
cimport numpy as np
import numpy as np


cdef extern:
    void gcorr(int* images, int* imager,
               float* npls, float* nplr, float* csmin,
               int* mfit,
               float* ddmx, float* ioffrq, float* nomoff,
               int* iacrej,
               float* streng, float* bfoffs, float* tlerrs, float* ddact)

def correlate(int[:,::1] search_img, int[:,::1] ref_img,
              float min_corr_strength=2.0, int fit_method=1, float max_dist=-1,
              float[::1] max_srch_offset=None, float[::1] nominal_offset=None):
    """ Compute correlation offsets between two images.

    Arguments
    ---------
    search_img: 2d integer array
    ref_img: 2d integer array
    min_corr_strength: (default 2.0)
    fit_method: (default 1)
    max_dist: (default -1)
    max_srch_offset: (default [64, 64])
    nominal_offset: (default [47, 47])

    Returns
    -------
    Length 5 array containing

        `okparam`, `corr_strength`, `best_offsets(2)`, `error_offsets(2)`
    """
    if max_srch_offset is None:
        max_srch_offset = np.array([64.0, 64.0])
    if nominal_offset is None:
        nominal_offset = np.array([47.0, 47.0])

    cdef float[2] search_size, ref_size

    cdef int okparam
    cdef float corr_strength
    cdef float[2] best_offsets
    cdef float[2] error_offsets
    cdef float diag_disp
    cdef float[6] result

    search_size = search_img.shape
    ref_size = ref_img.shape

    okparam = 1
    best_offsets[0] = 0.0
    best_offsets[1] = 0.0
    error_offsets[0] = 0.0
    error_offsets[1] = 0.0

    # From NSIDC, implemented in fortran
    gcorr(&search_img[0][0], &ref_img[0][0], &search_size[0], &ref_size[0],
        &min_corr_strength, &fit_method, &max_dist,
        &max_srch_offset[0], &nominal_offset[0],
        &okparam, &corr_strength, &best_offsets[0], &error_offsets[0],
        &diag_disp)

    result[0] = okparam
    result[1] = corr_strength
    result[2:4] = best_offsets
    result[4:] = error_offsets
    return result

# def mass_correlate():
#     return
