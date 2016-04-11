import numpy as np
from scipy import *
import matplotlib.pyplot as plt

# introduce all my variables and creates my matrix
errk = 1
errf = 1
tol = 1e-4
iters = 0
a = 4.0
D = 1.0
E = 0.7
vE = 0.6
h = 0.1
k = 0.1
n = int((a - (-1 * a)) / h)
nps = n + 1
flux = np.transpose(np.matrix(np.ones(n - 1)))
flux /= np.linalg.norm(flux)
a = [-1] * int(n - 2)
b = [2 + (E * h ** 2 / D)] * int(n - 1)
c = [-1] * int(n - 2)
A = np.matrix(np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1))
Q = (h ** 2 * vE / D) * flux


# my gauss-seidel solver
def GS(A, b, x, tol):
    Y = np.tril(A)
    Z = A - Y
    Max_Iter = 10000
    for i in range(Max_Iter):
        x_prev = x
        x = np.dot(np.linalg.inv(Y), b - np.dot(Z, x))
        error = np.linalg.norm(x - x_prev) / np.linalg.norm(x)
        if error <= tol:
            return x


# my power iteration solver
while errk > tol or errf > tol:
    iters += 1
    b = (1 / k) * Q
    fluxnew = GS(A, b, flux, 10E-6)
    Qnew = (h ** 2 * vE / D) * fluxnew
    knew = k * sum(Qnew) / sum(Q)
    errk = abs(knew - k)
    errf = np.linalg.norm(fluxnew - flux)
    k = knew
    flux = fluxnew
    Q = Qnew

print('Eigenvalue K = ' + str(k))
print('Number of Power Iterations = ', str(iters))
a = 4
xr = []
for i in range(1, n):
    xr.append((-1 * a + h * i))
fig = plt.figure(figsize=(16, 9))
plt.xlabel('x [cm]', size=18)
plt.ylabel(r'Flux [#/[$cm^{3}s$]]', size=18)
plt.title('Solution of Diffusion Equation Through Power Iteration', size=22)
plt.plot(xr, flux, "bx")
plt.show()
