import math
import pylab
import  numpy as np
from LagrangePolynom import LagrangePolynom
from ChordMethod import chord_method


class Ballistic:
    g = 9.81
    x, y, vx, vy = [], [], [], []
    angle, energy, v0, x0, y0 = [], [], [], [], []

    def __init__(self, energy, reflect_param, alpha, k, m, x0, y0, N, ground):
        self.angle = alpha
        self.energy = energy
        self.v0 = math.sqrt(2 * energy / m)
        self.reflection_param = reflect_param
        self.alpha = alpha
        self.k = k
        self.m = m
        self.x0 = x0
        self.y0 = y0
        self.N = N
        self.ground = ground

    def explicit_solution(self, T):
        expl_sol_x = lambda t: self.x0 + (self.m / self.k) * self.v0 * math.cos(self.alpha) - \
                               (1 - math.exp(-(self.k / self.m) * t))

        expl_sol_y = lambda t: self.y0 + (self.m / self.k) * (
                self.v0 * math.sin(self.alpha) + (self.g * self.m) / self.k) \
                               * (1 - math.exp(-(self.k / self.m) * t)) - (self.g * self.m) / self.k * t
        return expl_sol_x(T), expl_sol_y(T)

    def numerical_solution(self):
        t = 3
        h = t / self.N

        self.x.append(self.x0)
        self.y.append(self.y0)
        self.vx.append(self.v0 * math.cos(self.alpha))
        self.vy.append(self.v0 * math.sin(self.alpha))

        Fres = self.k * self.v0
        beta = self.alpha + math.pi

        for i in range(self.N - 1):
            k11 = self.vx[i] * h
            k12 = self.vy[i] * h
            k13 = (-Fres * math.cos(beta)) / self.m * h
            k14 = (-self.g - (Fres * math.sin(beta)) / self.m) * h

            k21 = (self.vx[i] + k13 / 2) * h
            k22 = (self.vy[i] + k14 / 2) * h
            k23 = ((-Fres * math.cos(beta)) / self.m) * h
            k24 = (-self.g - (Fres * math.sin(beta)) / self.m) * h

            k31 = (self.vx[i] + k23 / 2) * h
            k32 = (self.vy[i] + k24 / 2) * h
            k33 = ((-Fres * math.cos(beta)) / self.m) * h
            k34 = (-self.g - (Fres * math.sin(beta)) / self.m) * h

            k41 = (self.vx[i] + k33 / 2) * h
            k42 = (self.vy[i] + k34 / 2) * h
            k43 = ((-Fres * math.cos(beta)) / self.m) * h
            k44 = (-self.g - (Fres * math.sin(beta)) / self.m) * h

            self.x.append(self.x[i] + (k11 + 2 * k21 + 3 * k31 + k41) / 6)
            self.y.append(self.y[i] + (k12 + 2 * k22 + 3 * k32 + k42) / 6)
            self.vx.append(self.vx[i] + (k13 + 2 * k23 + 3 * k33 + k43) / 6)
            self.vy.append(self.vy[i] + (k14 + 2 * k24 + 3 * k34 + k44) / 6)

            Fres = self.k * math.sqrt(math.pow(self.vx[-1], 2) + math.pow(self.vy[-1], 2))
            beta = math.atan(self.vy[-1] / self.vx[-1])

            if self.y[-1] <= self.ground(self.x[-1]):
                break

        new_v = math.sqrt(self.vx[-1]**2 + self.vy[-1]**2)
        self.energy = (self.m * new_v ** 2)/2
        return self.x, self.y, self.vx, self.vy

    def reflection_initialize(self):
        h = 0.1
        interpolate_nodes, interpolate_values = self.x[-5:-1], self.x[-5:-1]
        new_nodes_to_interpolation = np.linspace(max(interpolate_nodes), min(interpolate_nodes), 20)
        lagrange_polynom = LagrangePolynom(np.array(interpolate_nodes), np.array(interpolate_values), new_nodes_to_interpolation)
        pol = lagrange_polynom.recalculate_in_new_nodes()
        res_polynom = lambda x: pol[0] + x * pol[1] + x ** 2 * pol[2] + x ** 3 * pol[3] + x ** 4 * pol[4]
        res_func = lambda x: pol[0] + x * pol[1] + x ** 2 * pol[2] + x ** 3 * pol[3] + x ** 4 * pol[4] - self.ground(x)
        root = chord_method(res_func, max(interpolate_nodes), min(interpolate_nodes))

        derivative_in_root_point_ground_function = (self.ground(root + h) - self.ground(root)) / h
        derivative_in_root_point_polynom = (res_polynom(root + h) - res_polynom(root) / h)

        angle_between_funcs = math.atan((derivative_in_root_point_polynom - derivative_in_root_point_ground_function)/\
                          (1 + derivative_in_root_point_polynom * derivative_in_root_point_ground_function))

        angle_between_ground_and_x_coordline = math.atan(derivative_in_root_point_ground_function)

        angle_falling = angle_between_funcs + angle_between_ground_and_x_coordline - 90

        new_angle = 90 - angle_falling
        new_x0, new_y0 = root, self.ground(root)

        self.init_after_reflection(new_angle, new_x0, new_y0)

    def init_after_reflection(self, new_angle, x0, y0):
        self.angle = new_angle
        self.energy = self.energy * self.reflection_param
        self.v0 = math.sqrt(2 * self.energy / self.m)
        self.x0 = x0
        self.y0 = y0

    def visualize(self):
        pylab.subplot(2, 1, 1)
        pylab.plot(self.x, self.y)
        pylab.title("График в координатах Х - У")

        # pylab.subplot(2, 1, 2)
        # pylab.plot(self.vx, self.vy)
        # pylab.plot("График в координатах скоростей")

        pylab.show()
