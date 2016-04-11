import numpy as np
import matplotlib.pyplot as plt
from math import *

def fluxes(h):
    # (-1)flux(i-1) + (2 + h^2/L^2)flux(i) (-1)flux(x)
    # (-1)x + (2 + h^2/L^2)y + (-1)z
    # Generate Matrix
    a = 4   #cm
    D = 1   #cm
    E = 0.2 #cm^-1
    s = 8   #n/(cm^3*s)
    L = sqrt(D / E)
    n = int(2*a / h)
    A = np.matrix(np.zeros((n + 1, n + 1)))
    x = np.repeat(0.0, n + 1)
    y = np.repeat(0.0, n + 1)
    z = np.repeat(0.0, n + 1)
    S = []
    X = np.repeat(0.0, n + 1)
    xran = np.linspace(-a, a, num=n + 1)
    for i in range(0, n + 1):
        A[i, i] = 2 + (h / L) ** 2
        S.append((h**2) * 8/D)
        y[i] = 2 + (h / L) ** 2
        X[i] = 0.0
        if i <= n - 1:
            A[i + 1, i] = -1.0
            A[i, i + 1] = -1.0
            x[i] = -1.0
            z[i] = -1.0
    # Solves Matrix
    xc, yc, zc, Sc = x, y, z, S     # copy the array

    for i in range(1, n):
        mc = xc[i]/yc[i-1]
        yc[i] = yc[i] - mc*zc[i-1]
        Sc[i] = Sc[i] - mc*Sc[i-1]
    xc[-1] = Sc[-1] / yc[-1]
    for i in range(n - 1, -1, -1):
        xc[i] = (Sc[i] - zc[i] * xc[i + 1]) / yc[i]
    flux = []
    for i in range(len(xran)):
        temp = ((exp(xran[i]/L)+exp(-xran[i]/L))/(exp(a/L)+exp(-a/L)))*-s/E + s/E
        flux.append(temp)
    xc[0] = 0
    err = []
    for i in range(len(xc)):
        if flux[i] != 0.0:
            err.append((xc[i]-flux[i])/flux[i])
        else:
            err.append(0)
    err = 1-max(err)
    return err, n

nerr = []
n = []
xr = [1, 0.5, 0.1, 0.05, 0.01]
for i in xr:
    nerr.append(fluxes(i)[0])
    n.append(fluxes(i)[1])
fig = plt.figure(figsize=(16, 9))
print("The Relative Error Decreases Drastically with More Meshes")
plt.xlabel('Number of Meshes', size=18)
plt.ylabel('Maximum Relative Error', size=18)
plt.title('Maximum Error vs Number of Meshes', size=22)
plt.plot(n, nerr, 'bo')
plt.show()
