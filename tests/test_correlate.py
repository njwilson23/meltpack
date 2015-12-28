import MELTPACK.correlate as mpc
import numpy as np
import matplotlib.pyplot as plt

x, y = np.meshgrid(np.linspace(0, 1), np.linspace(0, 1))

im0 = (30*(0.25-(x-0.5)**2) + 20*(0.25-(y-0.5)**2)).astype(np.float32)

im1 = np.ascontiguousarray(im0[5:35, 7:37])
im2 = np.ascontiguousarray(im0[10:40, 10:40])

print(im1[0,0], im1[0,1])

# subimage is called with the signature
#   subimage(x0, y0,
#            out_buffer, out_width, out_height,
#            in_buffer, in_width, in_height)

plt.subplot(1, 2, 1)
plt.imshow(im1, vmin=6, vmax=12)
plt.subplot(1, 2, 2)
plt.imshow(im2, vmin=6, vmax=12)

res = mpc.correlate(im1, im2, min_corr_strength=1.0, max_dist=10.0)
print(res)

plt.show()
