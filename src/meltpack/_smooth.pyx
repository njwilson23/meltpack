import numpy as np
cimport numpy as np
from libc.math cimport isnan
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def smooth5(np.ndarray[DTYPE_t, ndim=2] img, int niter=1):
    """ Apply a five-point smoothing kernel to *img* *niter* times.
    Arguments:
    img: np.ndarray
    niter: int
    """
    cdef int i
    for i in range(niter):
        img = _smooth5(img)
    return img

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _smooth5(np.ndarray[DTYPE_t, ndim=2] img):

    cdef int nx, ny
    cdef int i, j
    cdef double runsum
    cdef int count
    cdef np.ndarray[DTYPE_t, ndim=2] out

    ny = img.shape[0]
    nx = img.shape[1]
    out = np.zeros([ny, nx], dtype=DTYPE)

    for i in range(1, ny-1):
        for j in range(1, nx-1):
            count = 0
            runsum = 0.0

            if not isnan(img[i,j]):
                runsum += img[i,j]
                count += 1

                if not isnan(img[i-1,j]):
                    runsum += img[i-1,j]
                    count += 1
                if not isnan(img[i,j-1]):
                    runsum += img[i,j-1]
                    count += 1
                if not isnan(img[i,j+1]):
                    runsum += img[i,j+1]
                    count += 1
                if not isnan(img[i+1,j]):
                    runsum += img[i+1,j]
                    count += 1

            if count == 0:
                out[i,j] = np.nan
            else:
                out[i,j] = runsum/count

    out[0,:] = np.nan
    out[ny-1,:] = np.nan
    out[:,0] = np.nan
    out[:,nx-1] = np.nan
    return out

