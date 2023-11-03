import math
import pylab
import numpy as np
from LagrangePolynom import LagrangePolynom
from ChordMethod import chord_method


def find_polynom_derevative(pol):
    # pol: array of coefficients of polynom.
    # Only for polynom of 4-th order
    new_coefficients = []
    for i in range(len(pol)):
        new_coefficients.append(pol[i] * i)

    return lambda x: new_coefficients[1] + x * new_coefficients[2] + x ** 2 * new_coefficients[3]


class Ballistic:
    g = 9.81
    # arrays for local calculation variables
    x, y, vx, vy = [], [], [], []
    # initial value variables
    angle, energy, v0, x0, y0 = [], [], [], [], []
    # arrays for global trajectory coordinates
    global_x, global_y = [], []

    def __init__(self, energy, reflect_param, alpha, k, v, x0, y0, N, ground):
        self.angle = alpha
        self.energy = energy
        self.v0 = v
        self.reflection_param = reflect_param
        self.alpha = alpha
        self.k = k
        self.m = (2 * energy)/pow(v, 2)
        self.x0 = x0
        self.y0 = y0
        self.N = N
        self.ground = ground

        # Need to be here

        self.x = [] if self.x else []
        self.y = [] if self.y else []
        self.vx = [] if self.vx else []
        self.vy = [] if self.vy else []

    def explicit_solution(self, T):
        expl_sol_x = lambda t: self.x0 + (self.m / self.k) * self.v0 * math.cos(self.alpha) - \
                               (1 - math.exp(-(self.k / self.m) * t))

        expl_sol_y = lambda t: self.y0 + (self.m / self.k) * (
                self.v0 * math.sin(self.alpha) + (self.g * self.m) / self.k) \
                               * (1 - math.exp(-(self.k / self.m) * t)) - (self.g * self.m) / self.k * t
        return expl_sol_x(T), expl_sol_y(T)

    def numerical_solution(self):
        # how to pre-calculate t - value?
        t = 3
        h = t / self.N

        self.x.append(self.x0)
        self.y.append(self.y0)
        self.vx.append(self.v0 * math.cos(self.alpha))
        self.vy.append(self.v0 * math.sin(self.alpha))

        FRes = self.k * self.v0  # Resulting power of resistance
        beta = self.alpha + math.pi  # Resulting power vector directed against the vector of flight

        # Runge - Kutta method main iteration process
        for i in range(self.N - 1):
            k11 = self.vx[i] * h
            k12 = self.vy[i] * h
            k13 = (-FRes * math.cos(beta)) / self.m * h
            k14 = (-self.g - (FRes * math.sin(beta)) / self.m) * h

            k21 = (self.vx[i] + k13 / 2) * h
            k22 = (self.vy[i] + k14 / 2) * h
            k23 = ((-FRes * math.cos(beta)) / self.m) * h
            k24 = (-self.g - (FRes * math.sin(beta)) / self.m) * h

            k31 = (self.vx[i] + k23 / 2) * h
            k32 = (self.vy[i] + k24 / 2) * h
            k33 = ((-FRes * math.cos(beta)) / self.m) * h
            k34 = (-self.g - (FRes * math.sin(beta)) / self.m) * h

            k41 = (self.vx[i] + k33 / 2) * h
            k42 = (self.vy[i] + k34 / 2) * h
            k43 = ((-FRes * math.cos(beta)) / self.m) * h
            k44 = (-self.g - (FRes * math.sin(beta)) / self.m) * h

            self.x.append(self.x[i] + (k11 + 2 * k21 + 2 * k31 + k41) / 6)
            self.y.append(self.y[i] + (k12 + 2 * k22 + 2 * k32 + k42) / 6)
            self.vx.append(self.vx[i] + (k13 + 2 * k23 + 2 * k33 + k43) / 6)
            self.vy.append(self.vy[i] + (k14 + 2 * k24 + 2 * k34 + k44) / 6)

            # Recalculation of resistance vector length and his direction
            FRes = self.k * math.sqrt(math.pow(self.vx[-1], 2) + math.pow(self.vy[-1], 2))
            beta = math.atan(self.vy[-1] / self.vx[-1])

            # if point under the ground line we exit from solving process
            if self.y[-1] <= self.ground(self.x[-1]):
                break

        # Recalculate the velocity vector direction and kinetic energy value
        new_v = math.sqrt(self.vx[-1] ** 2 + self.vy[-1] ** 2)
        self.energy = (self.m * new_v ** 2) / 2

        # Write trajectory points to global arrays
        self.global_x += self.x
        self.global_y += self.y

        # Зачем тут ретурн?
        return self.x, self.y, self.vx, self.vy

    def reflection_initialize(self):
        # selection points for polynom building
        interpolate_nodes, interpolate_values = self.x[-5:-1], self.y[-5:-1]
        new_nodes_to_interpolation = np.linspace(max(interpolate_nodes), min(interpolate_nodes), 20)

        lagrange_polynom = LagrangePolynom(np.array(interpolate_nodes), np.array(interpolate_values),
                                           new_nodes_to_interpolation)

        pol = lagrange_polynom.L  # returns coefficients of built Lagrange polynom
        res_polynom = lambda x: pol[0] + x * pol[1] + x ** 2 * pol[2] + x ** 3 * pol[3]
        res_func = lambda x: res_polynom(x) - self.ground(x)
        root = chord_method(res_func, max(interpolate_nodes), min(interpolate_nodes))

        polynom_derivative_func = find_polynom_derevative(pol)
        derivative_in_root_point_polynom = polynom_derivative_func(root)

        h = 0.0001  # step for derivative calculation
        derivative_in_root_point_ground_function = (self.ground(root + h) - self.ground(root)) / h

        angle_between_funcs = math.atan((derivative_in_root_point_polynom - derivative_in_root_point_ground_function) / \
                                        (
                                                1 + derivative_in_root_point_polynom * derivative_in_root_point_ground_function))

        angle_between_ground_and_x_coordline = math.atan(derivative_in_root_point_ground_function)

        # test
        angle_falling = (math.atan(derivative_in_root_point_polynom))

        if angle_between_funcs < 0:
            angle_between_funcs += math.pi

        new_angle = 180 - math.degrees(angle_between_funcs)

        new_x0, new_y0 = root, self.ground(root)

        self.init_after_reflection(math.radians(new_angle), new_x0, new_y0)

    # only for 4 pointed polynom

    def init_after_reflection(self, new_angle, x0, y0):
        self.alpha = new_angle
        self.energy = self.energy * self.reflection_param
        self.v0 = math.sqrt(2 * self.energy / self.m)
        self.x0 = x0
        self.y0 = y0
        self.x, self.y, self.vx, self.vy = [], [], [], []

    def visualize(self):
        pylab.subplot(2, 1, 1)
        pylab.plot(self.global_x, self.global_y)

        x_mesh = np.linspace(min(self.global_x), max(self.global_x), 100)
        y_mesh = []
        for i in range(len(x_mesh)):
            y_mesh.append(self.ground(x_mesh[i]))
        pylab.plot(x_mesh, y_mesh, 'r')

        pylab.title("График в координатах Х - У")

        # pylab.subplot(2, 1, 2)
        # pylab.plot(self.vx, self.vy)
        # pylab.plot("График в координатах скоростей")

        pylab.show()
