import numpy as np
def MotionX(position, velocity, dt):
    for i in range(0, len(position)):
        position[i] = position[i] + velocity[i] * dt
    return position
