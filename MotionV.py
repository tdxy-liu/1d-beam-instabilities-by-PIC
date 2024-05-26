import numpy
from numpy import dot


def MotionV(velocity, qm, interp, efieldg, n_particles, dt):
    efieldg = numpy.transpose(numpy.mat(efieldg))
    if type(interp) == int:
        sum = interp * efieldg
        if n_particles == 0:
            velocity = 0
        else:
            for i in range(0, len(velocity)):
                velocity[i] = velocity[i] + 0.5 * qm * float(sum[0, 0]) * dt
    else:
        sum = dot(interp, efieldg)
        if n_particles == 0:
            velocity = 0
        else:
            for i in range(0, len(velocity)):
                velocity[i] = velocity[i] + 0.5 * qm * float(sum[i, 0]) * dt
    return velocity
