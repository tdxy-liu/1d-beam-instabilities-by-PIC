from numpy import real, zeros
from numpy.fft import fft, ifft


def Field(charge_density, Ng, dx, L, Ax, kap):
    rho_k = fft(charge_density)
    phi_k = zeros(Ng, dtype=complex)
    for i in range(0, len(rho_k)):
        phi_k[i] = rho_k[i] / (kap[i] ** 2)
    Phi = real(ifft(phi_k))
    Phi[Ng - 1] = 0
    Eg_k = -1j * kap * phi_k
    Eg = real(ifft(Eg_k))
    return Eg, Phi
