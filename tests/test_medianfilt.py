import numpy as np
import matplotlib.pyplot as plt
from meltpack.filt import medianfilt

orig = np.zeros([100, 100], dtype=np.float64)
orig[30:70,30:70] = 1.0
orig[60:75,45:55] = np.nan

orig += np.random.rand(100, 100)

sm = orig.copy()
for i in range(3):
    sm = medianfilt(sm)

plt.subplot(1, 2, 1)
plt.imshow(orig)
plt.subplot(1, 2, 2)
plt.imshow(sm)
plt.show()
