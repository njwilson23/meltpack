from MELTPACK.gaussmarkov import predict
import numpy as np
import matplotlib.pyplot as plt

n = 50
X, Y = np.meshgrid(np.linspace(0, 10, n), np.linspace(0, 10, n))
X = X.ravel()
Y = Y.ravel()
Z = 0.05*(X-2)**2 + (Y-3)
idx = np.random.randint(0, n**2, n**2//100)
Xo = X[idx]
Yo = Y[idx]
Zo = Z[idx]

# Structure function
def f(x):
    V = 10
    L = 10.0
    return V * np.exp(-x**2/L)

Zi, εi = predict(f, np.c_[X, Y], np.c_[Xo, Yo], Zo)

# Visualize results
fig = plt.figure()

ax_true = fig.add_subplot(2, 2, 1)
c = plt.scatter(X, Y, c=Z, marker="o")
plt.scatter(Xo, Yo, c="k", marker="o")
plt.colorbar(c)

ax_est = fig.add_subplot(2, 2, 2)
c = plt.scatter(X, Y, c=Zi, marker="o", vmin=Z.min(), vmax=Z.max())
plt.scatter(Xo, Yo, c="k", marker="o")
plt.colorbar(c)

ax_est = fig.add_subplot(2, 2, 3)
c = plt.scatter(X, Y, c=εi, marker="o")
plt.scatter(Xo, Yo, c="k", marker="o")
plt.colorbar(c)

ax_est = fig.add_subplot(2, 2, 4)
c = plt.scatter(X, Y, c=np.abs(Z-Zi), marker="o")
plt.scatter(Xo, Yo, c="k", marker="o")
plt.colorbar(c)
