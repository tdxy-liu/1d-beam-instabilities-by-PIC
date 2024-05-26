from numpy import append, zeros, dot
from numpy.linalg import inv
from sympy import solve


def Field(charge_density, Ng, dx, L, Ax, kap):
    Phi = zeros(Ng - 1)
    Ax = Ax.todense()
    result = dot(inv(Ax), (-charge_density[0:Ng - 1] * (dx ** 2)))
    for i in range(0, Ng - 1):
        Phi[i] = result[0, i]
    Phi = append(Phi, 0)
    Eg = (append(Phi[Ng - 1], Phi[0:Ng - 1]) - append(Phi[1:Ng], Phi[0])) / (2 * dx)
    return Phi, Eg