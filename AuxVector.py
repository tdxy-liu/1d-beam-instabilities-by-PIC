from numpy import linspace, append


def AuxVector(N_particles):
    p = linspace(0, N_particles - 1, N_particles)
    p = append(p, p)
    return p
