import math

from Ballistic import Ballistic

# Здесь инициализируются параметры задачи выполняются пероразования
# Точное решение задачи +
# Метод рунге Кутты(частный случай) +
# Вывод результата +
# Интерполяция +
# Метод секущих +
# Отражение снаряда +
# Сходимость Рунге-Кутты -
# Обратная задача

kinetic_energy = 100
reflection_energy_loss_param = 0.9
alpha = math.radians(27)  # angle: 38.6 main value
k = 1  # air resistance coefficient
v = 10  # speed
x0, y0 = 0, 0  # start coordinates
ground = lambda x: max(0, 1 - (x - 3)**2)  # ground level function (x**2)/10 tested

ball = Ballistic(kinetic_energy, reflection_energy_loss_param, alpha, k, v, x0, y0, 1000, ground)

for i in range(1):
    ball.numerical_solution()
    print(math.degrees(ball.alpha))
    print("Energy:" + str(ball.energy))
    try:
        ball.reflection_initialize()
    except:
        print("Процесс прерван на " + str(i) + " итерации")
    # ball.visualize()

ball.visualize()
