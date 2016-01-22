import MELTPACK.correlate as mpc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches

x, y = np.meshgrid(np.linspace(0, 1, 128), np.linspace(0, 1, 128))

im0 = (30*(0.25-(x-0.5)**2) + 20*(0.25-(y-0.5)**2)).astype(np.float32)
im0[40:85,60:66] = 11.0

xoff = 44
yoff = 78
nx = 32
ny = 32

ims = np.ascontiguousarray(im0[:,:])
imr = np.ascontiguousarray(im0[yoff:yoff+ny, xoff:xoff+nx])

# subimage is called with the signature
#   subimage(x0, y0,
#            out_buffer, out_width, out_height,
#            in_buffer, in_width, in_height)

plt.subplot(1, 2, 1)
plt.imshow(ims, vmin=3, vmax=12, origin="bottom")
rect = matplotlib.patches.Rectangle((xoff-1, yoff-1), nx, ny, edgecolor="k", facecolor="none")
plt.gca().add_patch(rect)
plt.title("search image")
plt.subplot(1, 2, 2)
plt.imshow(imr, vmin=3, vmax=12, origin="bottom")
plt.title("reference image")

res = mpc.correlate(ims, imr, min_corr_strength=2.0, max_dist=80.0,
                    max_srch_offset=np.array([64.0, 64.0], np.float32),
                    fit_method=1)
print("Exit code:     %d" % res[0])
print("Corr strength: %.5f" % res[1])
print("Best fit:      %.3f, %.3f" % (res[2], res[3]))
print("Error est.:    %.3f, %.3f" % (res[4], res[5]))
print("Diag. displ.:  %.3f" % res[6])

plt.show()
