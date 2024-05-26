from matplotlib.pyplot import scatter, pause, ion, show, cla
from numpy.fft import fftshift
from scipy.sparse import spdiags

from Field_FDM import Field
from MotionV import MotionV
#import field
#import interpolation
from numpy import arange, ones, transpose, pi, linspace, zeros, sum
from Charge_density import Charge_density
from MotionX import MotionX
from Charge import Charge
from InitialLoading import InitialLoading
from AuxVector import AuxVector
from InterpolationF import InterpolationF
from PeriodicBC import PeriodicBC
from MotionV_initial import MotionV_initial

L = 64
dt = 0.1
time = 100
Nt = int(time / dt)
Ng = 256

#beam 1
N1 = 10000
V01 = 5
Vth1 = 1
QM1 = -1
XP1 = 0
Mode1 = 0
WP = 1

#beam 2
N2 = 10000
V02 = -5
Vth2 = 1
QM2 = -1
XP2 = 0
Mode2 = 0

#Ion
IB = 20000

#Size of the cell
dx = L / Ng

#Size of the time step
itime = arange(0, time, dt)

#Charge density
Q1, Q2, rho_back = Charge(QM1, QM2, IB, N1, N2, L, WP)

#Initial conditions
xp1, vp1 = InitialLoading(N1, V01, Vth1, XP1, Mode1, L)
xp2, vp2 = InitialLoading(N2, V02, Vth2, XP2, Mode2, L)

#Auxiliarity vectors
p1 = AuxVector(N1)
p2 = AuxVector(N2)

#Poisson's equation preparation
un = ones((Ng - 1, 1))
Ax = spdiags((un * [1, -2, 1]).transpose(), [-1, 0, 1], Ng - 1, Ng - 1)
kap = (2 * pi / L) * linspace(-Ng / 2, Ng / 2 - 1, 256)
kap = transpose(kap)
kap = fftshift(kap)
kap[0] = 1

mat1 = 0
mat2 = 0
Eg = 0
Phi = 0

#Main loop
mom = zeros(Nt)
E_kin = zeros(Nt)
E_pot = zeros(Nt)

for it in range(0, Nt - 1):
    mom[it] = (Q1 / QM1) * (sum(vp1)) + (Q2 / QM2) * (sum(vp2))
    E_kin[it] = 0.5 * abs(Q1 / QM1) * (sum(vp1 ** 2)) + 0.5 * abs(Q2 / QM2) * (sum(vp2 ** 2))
    E_pot[it] = 0.5 * sum(Eg ** 2) * dx
    vp1 = MotionV(vp1, QM1, mat1, Eg, N1, dt)
    vp2 = MotionV(vp2, QM2, mat2, Eg, N2, dt)

    #Updating positions
    xp1 = MotionX(xp1, vp1, dt)
    xp2 = MotionX(xp2, vp2, dt)

    #Periodic boundary conditions
    xp1 = PeriodicBC(xp1, 0, L)
    xp2 = PeriodicBC(xp2, 0, L)

    #Interpolation functions
    mat1 = InterpolationF(dx, Ng, xp1, N1, p1)
    mat2 = InterpolationF(dx, Ng, xp2, N2, p2)

    #Charge density
    rho1 = Charge_density(Q1, mat1, dx)
    rho2 = Charge_density(Q2, mat2, dx)
    rhot = rho1 + rho2 + rho_back

    #Field equations
    Phi, Eg = Field(rhot, Ng, dx, L, Ax, kap)

    #Updating velocity
    vp1 = MotionV(vp1, QM1, mat1, Eg, N1, dt)
    vp2 = MotionV(vp2, QM2, mat2, Eg, N2, dt)

    ion()
    scatter(xp1, vp1, s=2)
    scatter(xp2, vp2, s=2)
    pause(0.1)
    cla()

show()
