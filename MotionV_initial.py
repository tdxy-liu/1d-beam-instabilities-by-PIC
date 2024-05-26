import numpy
def MotionV_initial(velocity, qm, interp, efieldg, n_particles, dt):
    if n_particles == 0:
        velocity = 0
    else:
        velocity = velocity + 0.5 * qm * interp * efieldg * dt
    return velocity
