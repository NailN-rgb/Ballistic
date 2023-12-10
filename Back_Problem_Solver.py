import math

import numpy as np

from Ballistic import Ballistic

import matplotlib.pyplot as plt

# variable for diverging angle segment values
N_angle = 5


class BackSolver:
    iter_count = 30

    def __init__(self, x_t, y_t, x_0, y_0, v, k, energy, ground, tol):
        self.x_true = x_t
        self.y_true = y_t
        self.x_start = x_0
        self.y_start = y_0
        self.v = v
        self.k = k
        self.kinetic_energy = energy
        self.ground = ground
        self.tol = tol
        self.angle_max = math.radians(90)
        self.angle_min = math.radians(0)

    # return sum of the deviation squares
    def squares_errors(self, x, y):
        return math.pow(x - self.x_true, 2) + math.pow(y - self.y_true, 2)

    def angle_comparison(self):
        angles = np.linspace(self.angle_min, self.angle_max, N_angle)
        errors = np.zeros([N_angle, 2])

        for i in range(len(angles)):
            # solve for each angle
            ball = Ballistic(self.kinetic_energy, 0.9, angles[i], self.k, self.v, self.x_start, self.y_start, 1000,
                             self.ground)
            ball.numerical_solution()

            try:
                ball.reflection_initialize()
            except:
                ball.x0, ball.y0 = 0, 0

            errors[i, 0] = angles[i]
            errors[i, 1] = self.squares_errors(ball.x0, ball.y0)

        # sort by errors
        errors = errors[errors[:, 1].argsort()]

        # take two angles with minimal error

        angle_value1 = errors[0, 0]
        angle_value2 = errors[1, 0]

        new_angle_min = min(angle_value1, angle_value2)
        new_angle_max = max(angle_value1, angle_value2)

        if new_angle_min == self.angle_min and new_angle_max == self.angle_max:
            print("Here is a more than 1 solution")
            if angle_value1 < angle_value2:
                new_angle_min = angle_value1
                new_angle_max = angle_value1 + (self.angle_max - self.angle_min) / N_angle
            else:
                new_angle_max = angle_value1
                new_angle_min = angle_value1 - (self.angle_max - self.angle_min) / N_angle

        self.angle_min = new_angle_min
        self.angle_max = new_angle_max

        return errors[0, 0], errors[0, 1]

    # for angle calculation
    def solve(self):
        for i in range(self.iter_count):
            ang, err = self.angle_comparison()

            if err <= self.tol:
                print("Angle Found")
                break
            # else:
            # print(res)

        print("Angle = " + str(math.degrees(ang)))

    def PreProccess(self):
        x1 = []
        x2 = []
        y = []

        for i in range(0, 90):
            ang = math.radians(i)
            ball = Ballistic(self.kinetic_energy, 0.9, ang, self.k, self.v, self.x_start, self.y_start, 1000,
                             self.ground)
            ball.numerical_solution()

            try:
                ball.reflection_initialize()
            except:
                ball.x0, ball.y0 = ball.x[-1], ball.y[-1]

            x1.append(ball.x0)
            x2.append(ball.y0)
            y.append(i)

        plt.figure(0)
        plt.plot(y, x1)
        plt.title("X1-Y Coordinates")
        plt.grid()
        plt.xlabel("theta")
        plt.ylabel("X")

        plt.figure(1)
        plt.plot(y, x2)
        plt.title("X2-Y Coordinates")
        plt.grid()
        plt.xlabel("theta")
        plt.ylabel("Y")
        plt.show()