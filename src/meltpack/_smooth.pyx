import numpy as np
cimport numpy as np
from libc.math cimport isnan
cimport cython

DOUBLE = np.float64
ctypedef np.float64_t DOUBLE_t

INT = np.int16
ctypedef np.int16_t INT_t

def smooth5_int(np.ndarray[INT_t, ndim=2] img, int niter=1, int nodata=-1):
    """ Apply a five-point smoothing kernel to integer *img* *niter* times.
    Arguments:
    img: np.ndarray[np.int16]
    niter: int
    nodata: int
        value used as a NODATA flag
    """
    cdef int i
    for i in range(niter):
        img = _smooth5_int(img, nodata)
    return img

def smooth5(np.ndarray[DOUBLE_t, ndim=2] img, int niter=1):
    """ Apply a five-point smoothing kernel to *img* *niter* times.
    Arguments:
    img: np.ndarray[np.float64]
    niter: int
    """
    cdef int i
    for i in range(niter):
        img = _smooth5_double(img)
    return img

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _smooth5_int(INT_t [:,:] img, int nodata):

    cdef int nx, ny
    cdef int i, j
    cdef double runsum
    cdef int count
    cdef np.ndarray[INT_t, ndim=2] out

    ny = img.shape[0]
    nx = img.shape[1]
    out = np.zeros([ny, nx], dtype=INT)
    cdef INT_t [:,:] out_view = out

    for i in range(1, ny-1):
        for j in range(1, nx-1):
            count = 0
            runsum = 0.0

            if img[i,j] != nodata:
                runsum += img[i,j]
                count += 1

                if img[i-1,j] != nodata:
                    runsum += img[i-1,j]
                    count += 1
                if img[i,j-1] != nodata:
                    runsum += img[i,j-1]
                    count += 1
                if img[i,j+1] != nodata:
                    runsum += img[i,j+1]
                    count += 1
                if img[i+1,j] != nodata:
                    runsum += img[i+1,j]
                    count += 1

            if count == 0:
                out_view[i,j] = nodata
            else:
                out_view[i,j] = round_int(runsum/count)

    out_view[0,:] = nodata
    out_view[ny-1,:] = nodata
    out_view[:,0] = nodata
    out_view[:,nx-1] = nodata
    return out

cdef int round_int(double a):
    if a%1 < 0.5:
        return int(a//1)
    else:
        return int(a//1)+1

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _smooth5_double(DOUBLE_t [:,:] img):

    cdef int nx, ny
    cdef int i, j
    cdef double runsum
    cdef int count
    cdef np.ndarray[DOUBLE_t, ndim=2] out

    ny = img.shape[0]
    nx = img.shape[1]
    out = np.zeros([ny, nx], dtype=DOUBLE)
    cdef DOUBLE_t [:,:] out_view = out

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
                out_view[i,j] = np.nan
            else:
                out_view[i,j] = runsum/count

    out_view[0,:] = np.nan
    out_view[ny-1,:] = np.nan
    out_view[:,0] = np.nan
    out_view[:,nx-1] = np.nan
    return out

