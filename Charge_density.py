import numpy


def Charge_density(charge, interp, dx):
    rho = numpy.zeros(256)
    sum = numpy.sum(interp, axis=0)
    for i in range(0, 256):
        rho[i] = (charge / dx) * sum[i]
    return rho
