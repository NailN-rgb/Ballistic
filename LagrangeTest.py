import numpy as np
import math
from LagrangePolynom import LagrangePolynom
import matplotlib.pyplot as plot


# x = np.array([2, 5, -6, 7, 4, 3, 8, 9, 1, -2], dtype=float)
# y = np.array([-1, 77, -297, 249, 33, 9, 389, 573, -3, -21], dtype=float)

def testfunc1(x_points):
    y_points = []
    for i in range(len(x_points)):
        y_points.append((x_points[i] + math.sin(x_points[i])) ** 2)
    return y_points


x = np.linspace(-5, 5, 3)
y = testfunc1(x)
xnew = np.linspace(np.min(x), np.max(x), 100)
pol = LagrangePolynom(x, y, xnew)
print(pol.L)
pol.visualize()
