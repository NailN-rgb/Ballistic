import math

from Ballistic import Ballistic

# Здесь инициализируются параметры задачи выполняются пероразования
# Точное решение задачи +
# Метод рунге Кутты(частный случай) +
# Вывод результата +
# Интерполяция +
# Метод секущих +
# Отражение снаряда

kinetic_energy = 20
reflection_energy_loss_param = 0.9
alpha = math.radians(38.6)  # angle
k = 1  # air resistance coeficient
m = 1  # mass
x0, y0 = 0, 0  # start coordinates
ground = lambda x: -x  # ground level function

ball = Ballistic(kinetic_energy, reflection_energy_loss_param, alpha, k, m, x0, y0, 100, ground)

for i in range(10):
    ball.numerical_solution()
    ball.reflection_initialize()

    if ball.energy < 10:
        break

ball.visualize()
