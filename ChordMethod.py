import math
import matplotlib.pyplot as plt


def chord_method(func, x0, x_i, N=100, eps=1e-4):
    x_next = x_i - (func(x_i) * (x_i - x0)) / (func(x_i) - func(x0))

    for i in range(N):
        x_i = x_next
        x_next = x_i - (func(x_i) * (x_i - x0)) / (func(x_i) - func(x0))

        if abs(x_next - x_i) < eps:
            print(str(i) + " iteration")
            break

        if i == N - 1:
            print("Cycle overfilled")

    return x_next


f = lambda x: math.exp(-x) - x
x0 = 1
xi = 20

print(chord_method(f, x0, xi))