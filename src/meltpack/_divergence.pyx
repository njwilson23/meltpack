import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def divergence(DTYPE_t dx, DTYPE_t dy,
               np.ndarray[DTYPE_t, ndim=2] h,
               np.ndarray[DTYPE_t, ndim=2] u,
               np.ndarray[DTYPE_t, ndim=2] v):
    """ Return the divergence of a vector field using upwind finite volumes.
    
    Arguments
    ---------
    dx: float
    dy: float
    h: np.ndarray
    u: np.ndarray
    v: np.ndarray
    """
    cdef unsigned int ny, nx
    cdef unsigned int i, j
    cdef DTYPE_t u_, v_
    cdef np.ndarray[DTYPE_t, ndim=2] div
    cdef np.ndarray[DTYPE_t, ndim=2] xfluxes, yfluxes

    # Allocate intermediate and output arrays
    ny = h.shape[0]
    nx = h.shape[1]
    xfluxes = np.zeros([ny, nx-1], dtype=DTYPE)
    yfluxes = np.zeros([ny-1, nx], dtype=DTYPE)
    div = np.nan*np.zeros([ny-2, nx-2], dtype=DTYPE)

    # Compute fluxes in x direction
    for i in range(ny):
        for j in range(nx-1):
            u_ = 0.5*(u[i,j] + u[i,j+1])
            if u_ > 0.0:
                xfluxes[i,j] = h[i,j]*u[i,j]
            elif u_ < 0.0:
                xfluxes[i,j] = h[i,j+1]*u[i,j+1]
            else:
                xfluxes[i,j] = 0.5*(h[i,j]*u[i,j] + h[i,j+1]*u[i,j+1])

    # Compute fluxes in y direction
    for i in range(ny-1):
        for j in range(nx):
            v_ = 0.5*(v[i,j] + v[i+1,j])
            if v_ > 0.0:
                yfluxes[i,j] = h[i,j]*v[i,j]
            elif v_ < 0.0:
                yfluxes[i,j] = h[i+1,j]*v[i+1,j]
            else:
                yfluxes[i,j] = 0.5*(h[i,j]*v[i,j] + h[i+1,j]*v[i+1,j])

    # Apply fluxes
    for i in range(1, ny-2):
        for j in range(1, nx-2):

            div[i,j] = (-xfluxes[i+1,j] + xfluxes[i+1,j+1])/dx \
                     + (-yfluxes[i,j+1] + yfluxes[i+1,j+1])/dy

    return div
