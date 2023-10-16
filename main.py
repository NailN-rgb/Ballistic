import math
import numpy as np
import matplotlib.pyplot as plt


v = 20  # start speed
alpha = math.radians(38.6)  # angle

k = 0.3 # air resistance coeficient
m = 1 # mass
x0 = 0  # start coordinates
y0 = 0  #
yB = 0  # last point height

# Explicit calculations (Test data)
y = lambda t: y0 + v * t * math.sin(alpha) - g * (t ** 2) / 2
x = lambda t: x0 + v * t * math.cos(alpha)

# Point A(max height)
tA = (v * math.sin(alpha)) / g
yA = y(tA)
xA = x(tA)

# Point B(max length)
tB = (v * math.sin(alpha) + math.sqrt(math.pow(v, 2) * math.pow(math.sin(alpha), 2) - 2 * g * (y0 + yB))) / g
yB2 = y(tB)
xB = x(tB)

print("Тестовое решение(формульное)")
print("Максимальная высота подъема: " + str(yA) + ", при координате x: " + str(xA) + ". Время: " + str(tA))
print("Ожимаемое значение высоты: " + str(yB) + ", действительная: " + str(yB2))
print("Дальность полета: " + str(xB) + ", время полета: " + str(tB))

# Numerical Calculations

tstart = 0
tend = 2.5
N = 100
h = (tend - tstart) / N
t = np.linspace(tstart, tend, N)


# Initial Conditions
x = []
y = []
vx = []
vy = []
x.append(x0)
y.append(y0)
vx.append(v * math.cos(alpha))
vy.append(v * math.sin(alpha))


for i in range(N - 1):
    k11 = vx[i] * h
    k12 = vy[i] * h
    k13 = 0
    k14 = -g * h

    k21 = (vx[i] + k13 / 2) * h
    k22 = (vy[i] + k14 / 2) * h
    k23 = 0
    k24 = -g * h

    k31 = (vx[i] + k23 / 2) * h
    k32 = (vy[i] + k24 / 2) * h
    k33 = 0
    k34 = -g * h

    k41 = (vx[i] + k33 / 2) * h
    k42 = (vy[i] + k34 / 2) * h
    k43 = 0
    k44 = -g * h

    x.append(x[i] + (k11 + 2 * k21 + 3 * k31 + k41) / 6)
    y.append(y[i] + (k12 + 2 * k22 + 3 * k32 + k42) / 6)
    vx.append(vx[i] + (k13 + 2 * k23 + 3 * k33 + k43) / 6)
    vy.append(vy[i] + (k14 + 2 * k24 + 3 * k34 + k44) / 6)

    if y[-1] <= yB:
        break

print("\nЭксперементальные данные:")
print("Дальность полета: " + str(x[len(x) - 1]) + ". Максимальная высота: " + str(max(y)) + ". Время полета: " +
       str(t[len(x)]))

xr, yr, vxr, vyr = [], [], [], []

xr.append(x0)
yr.append(y0)
vxr.append(v * math.cos(alpha))
vyr.append(v * math.sin(alpha))

Fres = k * v
beta = alpha + math.pi

for i in range(N - 1):
    k11 = vxr[i] * h
    k12 = vyr[i] * h
    k13 = (-Fres * math.cos(beta)) / m * h
    k14 = (-g - (Fres * math.sin(beta)) / m) * h

    k21 = (vxr[i] + k13 / 2) * h
    k22 = (vyr[i] + k14 / 2) * h
    k23 = ((-Fres * math.cos(beta)) / m) * h
    k24 = (-g - (Fres * math.sin(beta)) / m) * h

    k31 = (vxr[i] + k23 / 2) * h
    k32 = (vyr[i] + k24 / 2) * h
    k33 = ((-Fres * math.cos(beta)) / m) * h
    k34 = (-g - (Fres * math.sin(beta)) / m) * h

    k41 = (vxr[i] + k33 / 2) * h
    k42 = (vyr[i] + k34 / 2) * h
    k43 = ((-Fres * math.cos(beta)) / m) * h
    k44 = (-g - (Fres * math.sin(beta)) / m) * h

    xr.append(xr[i] + (k11 + 2 * k21 + 3 * k31 + k41) / 6)
    yr.append(yr[i] + (k12 + 2 * k22 + 3 * k32 + k42) / 6)
    vxr.append(vxr[i] + (k13 + 2 * k23 + 3 * k33 + k43) / 6)
    vyr.append(vyr[i] + (k14 + 2 * k24 + 3 * k34 + k44) / 6)

    Fres = k * math.sqrt(math.pow(vxr[-1], 2) + math.pow(vyr[-1], 2))
    beta = math.atan(vyr[-1]/vxr[-1])

    if yr[-1] <= yB:
        break

# Theoretical Calculations
xtheor = lambda t: v * math.cos(alpha) * (m / k) * (1 - math.exp(-k * t))
ytheor = lambda t: (v*math.sin(alpha) * (m/k) + m**2 * g / k**2) * (1 - math.exp(-k * t)) - m*g/k * t


print("\nС учетом сопротивления воздуха:")
print("Теоретически дальность полета: " + str(xtheor(t[len(xr)])) + ". Высота в последний момент времени: " +
      str(ytheor(t[len(xr)])))
print("Практически. Дальность полета: " + str(xr[-1]) + ". Максимальная высота: " + str(max(yr)) + ". Время полета: " +
      str(t[len(xr)]))

plt.plot(x, y)
plt.plot(xr, yr, color = 'tab:red')
plt.show()



