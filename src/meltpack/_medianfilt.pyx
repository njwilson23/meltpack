import numpy as np
cimport numpy as np
from libc.math cimport isnan
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def medianfilt(np.ndarray[DTYPE_t, ndim=2] img):
    """ Perform a nodata-aware median filtering operation on a numpy grid. """

    cdef int nx, ny
    cdef int i, j
    cdef double runsum
    cdef int count
    cdef np.ndarray[DTYPE_t, ndim=2] out

    ny = img.shape[0]
    nx = img.shape[1]
    out = np.zeros([ny, nx], dtype=DTYPE)

    arr = np.zeros(9, dtype=DTYPE)
    cdef double[:] buf = arr

    for i in range(1, ny-1):
        for j in range(1, nx-1):
            count = 0
            runsum = 0.0

            if not isnan(img[i,j]):
                count += 1
                buf[0] = img[i,j]

                if not isnan(img[i-1,j-1]):
                    buf[count] = img[i-1,j-1]
                    count += 1
                if not isnan(img[i-1,j]):
                    buf[count] = img[i-1,j]
                    count += 1
                if not isnan(img[i-1,j+1]):
                    buf[count] = img[i-1,j+1]
                    count += 1
                if not isnan(img[i,j-1]):
                    buf[count] = img[i,j-1]
                    count += 1
                if not isnan(img[i,j+1]):
                    buf[count] = img[i,j+1]
                    count += 1
                if not isnan(img[i+1,j-1]):
                    buf[count] = img[i+1,j-1]
                    count += 1
                if not isnan(img[i+1,j]):
                    buf[count] = img[i+1,j]
                    count += 1
                if not isnan(img[i+1,j+1]):
                    buf[count] = img[i+1,j+1]
                    count += 1

            if count == 0:
                out[i,j] = np.nan
            else:
                out[i,j] = np.median(arr[:count])

    out[0,:] = img[0,:]
    out[ny-1,:] = img[ny-1,:]
    out[:,0] = img[:,0]
    out[:,nx-1] = img[:,nx-1]
    return out
