from numpy import floor, append, zeros
from scipy.sparse import csr_matrix, coo_matrix

from PeriodicBC import PeriodicBC


def InterpolationF(Dx, Ng, position, N_particles, aux_vector):
    if N_particles == 0:
        interp = 0

    else:
        im = zeros(len(position))
        ip = zeros(len(position))
        Fim = zeros(len(position))
        Fip = zeros(len(position))
        for i in range(0, len(position)):
            im[i] = floor(position[i] / Dx) - 1
            ip[i] = im[i] + 1
        Project = append(im, ip)
        Project = PeriodicBC(Project, 0, Ng)
        for i in range(0, len(position)):
            Fim[i] = 1 - abs((position[i] / Dx) - im[i])
            Fip[i] = 1 - Fim[i]
        Fraction = append(Fim, Fip)
        interp = zeros((N_particles, Ng))
        for i in range(0, len(aux_vector)):
            interp[int(aux_vector[i])][int(Project[i])] = Fraction[i]

    return interp
