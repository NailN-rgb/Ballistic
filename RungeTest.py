import math
import numpy as np
import matplotlib.pyplot as plt


def Runge_Kutta(aValue, x_mesh, equation, h_x, N_x):
    solution = [aValue]

    h2 = (x_mesh[1] - x_mesh[0]) / 2
    for i in range(0, N_x):
        k1 = equation(x_mesh[i])
        k2 = equation(x_mesh[i] + h2)
        k3 = equation(x_mesh[i] + h2)
        k4 = equation(x_mesh[i] + 2 * h2)
        solution.append(solution[-1] + (k1 + 2 * (k2 + k3) + k4) * h_x / 6)

    return solution


N = 25
a = 0
b = 1
h = (b - a) / N

equation_1 = lambda x: math.cos(x)
eq_1_an_solution = lambda x: math.sin(x)
equation_2 = lambda x: math.exp(x) * math.sin(math.sqrt(1 - x))


x1 = np.linspace(a, b, N)
x2 = np.linspace(a, b, N * 2)

sol1 = Runge_Kutta(eq_1_an_solution(a), x1, equation_1, h, N)
sol2 = Runge_Kutta(eq_1_an_solution(a), x2, equation_1, h / 2, N * 2)

res = []
for i in range(len(sol2)):
    if i % 2 == 0:
        idx = int(i/2)
        res.append(abs(sol2[i] - sol1[idx]))

# an_sol = eq_1_an_solution(x1)
# for i in range(len(an_sol)):
#     print(abs(an_sol[i] - sol1[i]))
#print(abs(sol2[-3] - sol1[-2]) / (h ** 4))
print(max(res)/(h**4))
