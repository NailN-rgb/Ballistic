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
ground = lambda x: (x**2)/10   # ground level function

ball = Ballistic(kinetic_energy, reflection_energy_loss_param, alpha, k, m, x0, y0, 1000, ground)

for i in range(10):
    ball.numerical_solution()
    print(math.degrees(ball.alpha))
    print("Energy:" + str(ball.energy))
    try:
        ball.reflection_initialize()
    except:
        print("Процесс прерван на " + str(i) + " итерации")
    # ball.visualize()

ball.visualize()
