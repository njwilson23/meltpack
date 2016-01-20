import numpy as np
import matplotlib.pyplot as plt
from meltpack.smooth import smooth5

orig = np.zeros([100, 100], dtype=np.float64)
orig[30:70,30:70] = 1.0
orig[60:75,45:55] = np.nan
sm = smooth5(orig, 1000)

plt.subplot(1, 2, 1)
plt.imshow(orig)
plt.subplot(1, 2, 2)
plt.imshow(sm)
plt.show()
