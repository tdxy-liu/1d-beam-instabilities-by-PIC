import numpy
from numpy import linspace
from numpy.random import randn


def InitialLoading(N_particles, V0, Vth, Amplitude, Mode, L):
    position = numpy.transpose(linspace(0, L - L / N_particles, N_particles))
    random = list(randn(N_particles, 1))
    velocity = numpy.zeros(N_particles)
    for i in range(0, len(velocity)):
        velocity[i] = Vth * random[i]
    velocity = velocity + V0
    if Amplitude != 0:
        position = position + Amplitude * numpy.cos(2 * numpy.pi * Mode * (position / L))

    return position, velocity
