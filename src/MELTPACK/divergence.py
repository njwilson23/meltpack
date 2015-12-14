""" Use optimization techniques to compute flux or velocity divergence from data.
"""

from ._divergence import divergence

# import karta
# import numpy as np
# import meltfuncs
# import scipy.optimize as scopt

# def _devec_grid(grid):
#     return grid.values.ravel(order='F').astype(np.float64)
# 
# def solve_fluxdiv(h, u, v, alpha, options):
#     """
#     Arguments:
#     ----------
#     h, u, v: karta.RegularGrid, thickness and velocity components
#     alpha: float, smoothing parameter
#     options: dict, solver-specific options
#     """
# 
#     ny, nx = h.size
#     dx, dy = h.transform[2:4]
#     h_ = _devec_grid(h)
#     u_ = _devec_grid(u)
#     v_ = _devec_grid(v)
#     mask_ = (~np.isnan(h_) & ~np.isnan(u_) & ~np.isnan(v_)).astype(np.int32)
#     epsh_ = 5*np.ones_like(h)
#     epsu_ = 1*np.ones_like(h)
#     epsv_ = 1*np.ones_like(h)
# 
#     def f(v):
#         hg_ = v[:nx*ny]
#         ug_ = v[nx*ny:2*nx*ny]
#         vg_ = v[2*nx*ny:]
#         return meltfuncs.functional(nx, ny, dx, dy, mask_, h_, u_, v_,
#                                     epsh_, epsu_, epsv_, hg_, ug_, vg_, alpha)
# 
#     def g(v):
#         hg_ = v[:nx*ny]
#         ug_ = v[nx*ny:2*nx*ny]
#         vg_ = v[2*nx*ny:]
#         dJ = meltfuncs.gradient(nx, ny, dx, dy, mask_, h_, u_, v_,
#                                 epsh_, epsu_, epsv_, hg_, ug_, vg_, alpha)
#         return np.hstack(dJ)
# 
#     guess = np.hstack([h_, u_, v_])
#     result = scopt.minimize(f, jac=g, method="Newton-CG", **options)
#     return result
# 
# def solve_veldiv(h, u, v, alpha, options):
#     """
#     Arguments:
#     ----------
#     h, u, v: karta.RegularGrid, thickness and velocity components
#     alpha: float, smoothing parameter
#     options: dict, solver-specific options
#     """
# 
#     ny, nx = h.size
#     dx, dy = h.transform[2:4]
#     h_ = _devec_grid(h)
#     u_ = _devec_grid(u)
#     v_ = _devec_grid(v)
#     mask_ = (~np.isnan(h_) & ~np.isnan(u_) & ~np.isnan(v_)).astype(np.int32)
#     epsh_ = 5*np.ones_like(h)
#     epsu_ = 1*np.ones_like(h)
#     epsv_ = 1*np.ones_like(h)
# 
#     def f(v):
#         hg_ = v[:nx*ny]
#         ug_ = v[nx*ny:2*nx*ny]
#         vg_ = v[2*nx*ny:]
#         return meltfuncs.functional_veldiv(nx, ny, dx, dy, mask_, h_, u_, v_,
#                                            epsh_, epsu_, epsv_, hg_, ug_, vg_, alpha)
# 
#     def g(v):
#         hg_ = v[:nx*ny]
#         ug_ = v[nx*ny:2*nx*ny]
#         vg_ = v[2*nx*ny:]
#         dJ = meltfuncs.gradient_veldiv(nx, ny, dx, dy, mask_, h_, u_, v_,
#                                        epsh_, epsu_, epsv_, hg_, ug_, vg_, alpha)
#         return np.hstack(dJ)
# 
#     guess = np.hstack([h_, u_, v_])
#     result = scopt.minimize(f, jac=g, method="Newton-CG", **options)
#     return result

