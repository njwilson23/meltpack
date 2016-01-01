# import cython
cimport numpy as np
import numpy as np

cdef extern:
    void gcorr_(float* images, float* imager,
               float* npls, float* nplr, float* csmin,
               int* mfit,
               float* ddmx, float* ioffrq, float* nomoff,
               int* iacrej,
               float* strength, float* bfoffs, float* tlerrs, float* ddact)


# Tie point active flag accept/reject codes
# -----------------------------------------
# TPL_STATE_OK_MANUAL       0   Accept manually extracted tie point
# TPL_STATE_OK_AUTO         1   Accept auto registered tie point
# TPL_STATE_BAD_EDGE        2   Reject; correlation peak too near edge
#                               of search area.
# TPL_STATE_BAD_SUBPEAK     3   Reject; subsidiary peak comparable in
#                               strength to main peak
# TPL_STATE_BAD_PEAK        4   Reject; strength of peak below minimum
#                               specified by user
# TPL_STATE_BAD_DISP        5   Reject; diagonal displacement from
#                               nominal location exceeds maximum
#                               specified by user
# TPL_STATE_BAD_PARAM       6   Correlation not attempted; error in
#                               correlation parameters
# TPL_STATE_BAD_USER        7   Reject; manually by user
# TPL_STATE_BAD_MODEL       8   Reject; automatically by the modeling
#                               process
# TPL_STATE_BAD_USER_MODEL  9   Reject; manually by user during the
#                               modeling process
# TPL_STATE_OK_NOW         10   Prev. rejected tie point re-accepted
#                               manually by user.  Must add original
#                               active code of tie point (0 - 9) to
#                               get actual flag value (10 - 19).

def correlate(float[:,::1] srch_img, float[:,::1] ref_img,
              float min_corr_strength=2.0, int fit_method=1, float max_dist=-1,
              float[:] max_srch_offset=None, float[:] nominal_offset=None):
    """ Compute correlation offsets between two images.

    Arguments
    ---------
    srch_img: 2d float array
    ref_img: 2d float array
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
    cdef float[2] srch_size, ref_size
    cdef int okparam
    cdef float corr_strength
    cdef float[:] best_offsets
    cdef float[:] error_offsets
    cdef float diag_disp
    cdef float[:] result

    error_offsets = np.empty(2, dtype=np.float32)
    best_offsets = np.empty(2, dtype=np.float32)
    result = np.empty(7, dtype=np.float32)

    srch_size[0] = srch_img.shape[1]
    srch_size[1] = srch_img.shape[0]
    ref_size[0] = ref_img.shape[1]
    ref_size[1] = ref_img.shape[0]

    if max_srch_offset is None:
        max_srch_offset = np.array([64.0, 64.0], dtype=np.float32)
        # max_srch_offset[0] = 64.0
        # max_srch_offset[1] = 64.0
    if nominal_offset is None:
        nominal_offset = 0.5*np.array([srch_size[0]-ref_size[0],
                                       srch_size[1]-ref_size[1]], dtype=np.float32)
        # nominal_offset[0] = 0.5*(srch_size[0]-ref_size[0])
        # nominal_offset[1] = 0.5*(srch_size[1]-ref_size[1])
    okparam = 1
    best_offsets[0] = 0.0
    best_offsets[1] = 0.0
    error_offsets[0] = 0.0
    error_offsets[1] = 0.0

    # From NSIDC, implemented in fortran
    # Name mangled according to the GNU/Intel trailing underscore convention
    gcorr_(&srch_img[0][0], &ref_img[0][0], srch_size, ref_size,
        &min_corr_strength, &fit_method, &max_dist,
        &max_srch_offset[0], &nominal_offset[0],
        &okparam, &corr_strength, &best_offsets[0], &error_offsets[0], &diag_disp)

    result[0] = okparam
    result[1] = corr_strength
    result[2:4] = best_offsets
    result[4:6] = error_offsets
    result[6] = diag_disp
    return result

