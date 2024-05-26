def PeriodicBC(variable, LB, UB):
    # Periodic boundary conditions
    for index, value in enumerate(variable):
        if value < LB:
            variable[index] = value + UB
        elif value > UB:
            variable[index] = value - UB
        else:
            continue
    return variable
