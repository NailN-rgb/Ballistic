import math
import numpy as np
import matplotlib.pyplot as plt

from Back_Problem_Solver import BackSolver

# initialize the back problem variables
# need to find angle

# impact point

x_true, y_true = 1, 1

# shoot point coordinates
x_0, y_0 = 0, 0

# Energy and speed module

kinetic_energy = 100
v = 25

# air resistance coefficient

k = 0.5

# ground level function

ground = lambda x: max(0, 1 - (x-5)**2)

tol = 0.1

ball = BackSolver(x_true, y_true, x_0, y_0, v, k, kinetic_energy, ground, tol)
ball.solve()
ball.PreProccess()




