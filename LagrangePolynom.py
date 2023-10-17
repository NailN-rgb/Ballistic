import numpy as np
import math
import matplotlib.pyplot as plt


def polynom_coefficients(_pl, _xi):
    n = int(_xi.shape[0])
    coefficients = np.empty((n, 2))

    for i in range(n):
        if i == _pl:
            coefficients[i][0] = float('inf')
            coefficients[i][1] = float('inf')
        else:
            coefficients[i][0] = 1 / (_xi[_pl] - _xi[i])
            coefficients[i][1] = -_xi[i] / (_xi[_pl] - _xi[i])
    filt_coefficients = np.empty((n - 1, 2))

    j = 0

    for i in range(n):
        if coefficients[i][0] != float('inf'):
            filt_coefficients[j][0] = coefficients[i][1]
            filt_coefficients[j][1] = coefficients[i][0]
            j += 1

    return filt_coefficients


def polynomal_l(_xi):
    n = int(_xi.shape[0])
    pli = np.zeros((n, n))

    for pl in range(n):
        coefficients = polynom_coefficients(pl, _xi)
        for i in range(1, n - 1):
            if i == 1:
                pli[pl][0] = coefficients[i - 1][0] * coefficients[i][0]
                pli[pl][1] = coefficients[i - 1][1] * coefficients[i][0] + \
                             coefficients[i][1] * coefficients[i - 1][0]
                pli[pl][2] = coefficients[i - 1][1] * coefficients[i][1]

            else:
                clone_pli = np.zeros(n)
                for val in range(n):
                    clone_pli[val] = pli[pl][val]
                zeros_pli = np.zeros(n)
                for j in range(n - 1):
                    product_1 = clone_pli[j] * coefficients[i][0]
                    product_2 = clone_pli[j] * coefficients[i][1]
                    zeros_pli[j] += product_1
                    zeros_pli[j + 1] += product_2
                for val in range(n):
                    pli[pl][val] = zeros_pli[val]
    return pli


def get_polynom(_xi, _yi):
    n = int(_xi.shape[0])
    polyn_l = polynomal_l(_xi)

    for i in range(n):
        for j in range(n):
            polyn_l[i][j] *= _yi[i]
    L = np.zeros(n)
    for i in range(n):
        for j in range(n):
            L[i] += polyn_l[j][i]
    return L


class LagrangePolynom:
    z = 0
    L = 0

    def __init__(self, x, y, new_x):
        self.x = x
        self.y = y
        self.new_x = new_x
        self.L = get_polynom(self.x, self.y)

    def recalculate_in_new_nodes(self):
        res = np.zeros(len(self.new_x))

        for i in range(len(self.L)):
            for j in range(len(self.new_x)):
                res[j] += self.L[i] * math.pow(self.new_x[j], i)

        return res

    def visualize(self):
        plt.plot(self.x, self.y, 'o', self.new_x, self.recalculate_in_new_nodes())
        plt.grid(True)
        plt.show()
