u"""Решение параболической задачи."""
from __future__ import print_function  # Форматный вывод
import numpy
import math
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt


def tridiagonalMatrixAlgorithm(a, b, c, d):
    u"""Решение трехдиагональной матрицы методом прогонки.

    На вход поступают массивы a,b,c - под главной диагональю, главная диагональ
    над главной диагональю, d - массив свободных членов
    Возвращает массив решений X
    """
    n = b.shape[0]
    alpha = numpy.empty(n)
    alpha[1] = -c[1]/b[1]
    beta = numpy.empty(n)
    beta[1] = d[1]/b[1]
    for i in range(2, n-1):
        alpha[i] = (-c[i] / (b[i] + a[i] * alpha[i-1]))
        beta[i] = ((d[i] - a[i] * beta[i-1]) / (b[i] + a[i] * alpha[i-1]))
    x = numpy.empty(n)
    x[n-1] = ((d[n-1]-a[n-1]*beta[n-2])/(b[n-1]+a[n-1]*alpha[n-2]))
    for i in range(n-2, -1, -1):
        x[i] = alpha[i] * x[i+1] + beta[i]
    return x


def parabol_1(a, t, h, l, T):
    u"""
    t,h - шаги, T,l -границы. func1 начально краевые условия по t.

    Возвращает КРС с начальными условиями, КРС с аналитическими значениями
    количество временных слоев, количество разбиений
    """
    def func1(t, x, l):
        u"""Значение искомой функции на границах."""
        if x == 0:
            return 0
        if x == l:
            return 0
        if t == 0:
            return math.sin(2 * math.pi * x)
        print(t, x, l)

    def analytical(t1, x1, a):
        u"""Аналитическое решение уравнения."""
        return (math.exp(-4 * math.pi ** 2 * a * t1) *
                math.sin(2 * math.pi * x1))

    #  Получение количества ячеек по T
    n = math.ceil(T/t) + 1
    # Получение количества ячеек по X
    m = math.ceil(l/h) + 1
    table = numpy.ones([n, m], dtype=float)
    analyticalTable = numpy.ones([n, m], dtype=float)
    x = numpy.empty(n)
    for i in range(n):
        for j in range(m):
            analyticalTable[i, j] = float(analytical(i * t, j * h, a * a))
            # начальные условия по X
            table[0, j] = func1(0 * t, j * h, l)
        # Начальные условия по t
        table[i, 0] = func1(i * t, 0 * h, l)
        table[i, m-1] = func1(i * t, (m-1) * h, l)
        x[i] = h * i
    y = numpy.empty(m)
    # Не внесено в верхний цикл, так как он выполняется большое количество раз
    for i in range(m):
        y[i] = t * i
    return table, analyticalTable, x, y


def explicit_finite_difference(a, t, h, l, T):
    u"""Решение параболического уравнения с помощью явной схемы.

    Принимает на вход значения a, сеточные шаги и границы условия
    """
    print("explicit")
    sigma = a ** 2 * t / h ** 2
    u, analyticalTable, n, m = parabol_1(a, t, h, l, T)
    if sigma > 0.5:
        print("Схема неустойчива")
        return u, analyticalTable, n, m
    for k in range(1, len(n)):
        for j in range(1, len(m)-1):
            u[k, j] = (sigma * u[k-1, j+1] + (1 - 2 * sigma) * u[k-1, j] +
                       sigma * u[k, j-1])
    for i in range(len(n)):
        for j in range(len(m)):
            print("%-11.4f" % (u[i, j]), end=' ')
        print()
    print()
    n, m = numpy.meshgrid(n, m)
    return u, analyticalTable, n, m


def implicit_finite_difference(a1, t, h, l, T):
    u"""Решение параболического уравнения с помощью неявной схемы.

    Принимает на вход значения a, сеточные шаги и границы условия
    """
    u, analyticalTable, x, y = parabol_1(a1, t, h, l, T)
    print("implicit")
    n = len(x)
    m = len(y)
    a = numpy.zeros(n-1)
    b = numpy.zeros(n-1)
    c = numpy.zeros(n-1)
    d = numpy.zeros(n-1)
    sigma = float(a1 * a1 * t / (h * h))
    # вынесено из цикла по слоям, т.к. во время итераций не изменяется
    for j in range(1, n-1):  # Кроме границ,
        # Создаем трехдиагональную матрицу для всех временных слоев
        if j > 1:
            a[j] = sigma
        if (j < (n-2)):
            c[j] = sigma
        b[j] = - (1 + 2 * sigma)
    for k in range(0, m-1):  # Временные слои
        for j in range(1, n-1):  # Кроме границ
            # Создаем трехдиагональную матрицу для всех временных слоев
            if j in range(1, n-2):
                d[j] = -u[k][j]
        d[1] = -(u[k][1] + sigma * u[k+1][0])
        d[n-2] = -(u[k][n-2] + sigma * u[k+1][n-1])
        x1 = tridiagonalMatrixAlgorithm(a, b, c, d)
        # Подставляем найденные задачи в КРС
        for j in range(1, n-1):
            u[k+1][j] = x1[j]
    for i in range(len(x)):
        for j in range(len(y)):
            print("%-11.4f" % (u[i, j]), end=' ')
        print()
    print()
    x, y = numpy.meshgrid(x, y)
    return u, analyticalTable, x, y


def crankNicolson(a1, t, h, l, T):
    u"""Решение параболического уравнения с помощью явно-неявной схемы.

    Принимает на вход значения a, сеточные шаги и границы условия
    """
    u, analyticalTable, x, y = parabol_1(a1, t, h, l, T)
    print("crankNicolson")
    n = len(x)
    m = len(y)
    a = numpy.zeros(n-1)
    b = numpy.zeros(n-1)
    c = numpy.zeros(n-1)
    d = numpy.zeros(n-1)
    sigma_2 = float(a1 * a1 * t / (h * h * 2))
    # вынесено из цикла по слоям, т.к. во время итераций не изменяется
    for j in range(1, n-1):  # Кроме границ
        # Создаем трехдиагональную матрицу для всех временных слоев
        if j > 1:
            a[j] = sigma_2
        if (j < (n-2)):
            c[j] = sigma_2
        b[j] = - (sigma_2 * 2 + 1)
    for k in range(0, m-1):  # Временные слои
        for j in range(1, n-1):  # Кроме границ
            # Создаем трехдиагональную матрицу для всех временных слоев
            if j in range(1, n-2):
                d[j] = (-u[k, j] - sigma_2 * (u[k, j+1] - 2 * u[k, j]
                        + u[k, j-1]))
        d[1] = ((-u[k][1] - sigma_2 * (u[k, 2] - 2 * u[k, 1] + u[k, 0])
                - sigma_2 * u[k+1][0]))
        d[n-2] = ((-u[k][n-2] - sigma_2 * (u[k, n-1] - 2 * u[k, n-2] +
                  u[k, n-3]) - sigma_2 * u[k+1][n-1]))
        x1 = tridiagonalMatrixAlgorithm(a, b, c, d)
        # Подставляем найденные задачи в КРС
        for j in range(1, n-1):
            u[k+1][j] = x1[j]
    for i in range(len(x)):
        for j in range(len(y)):
            print("%-11.4f" % (u[i, j]), end=' ')
        print()
    print()
    x, y = numpy.meshgrid(x, y)
    return u, analyticalTable, x, y


a = 0.025
t = 0.1
h = 0.1
T = 1.
L = 1.
print("a =", a, "t =", t, "h =", h)
u1, analyticalT, xgrid, ygrid = explicit_finite_difference(a, t, h, T, L)
print("Погрешность решения:")
deltaExplict = abs(analyticalT - u1)
for i in range(len(xgrid)):
    for j in range(len(ygrid)):
        print("%-11.4f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-11.6f" % deltaExplict.max())
print()

fig = pylab.figure()
axes = Axes3D(fig)
axes.plot_surface(xgrid, ygrid, u1, rstride=1, cstride=1, cmap=cm.jet)
plt.title(u"Явная схема")
pylab.show()
fig.savefig("parabЯвная.png")

fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Аналитическое решение")
axes.plot_surface(xgrid, ygrid, analyticalT, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("parabАналит.png")


u, analyticalT, xgrid, ygrid = implicit_finite_difference(a, t, h, T, L)
print("Погрешность решения:")
deltaExplict = abs(analyticalT - u)
for i in range(len(xgrid)):
    for j in range(len(ygrid)):
        print("%-11.4f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-11.6f" % deltaExplict.max())
print()
fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Неявная схема")
axes.plot_surface(xgrid, ygrid, u, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("parabНеявная.png")

u, analyticalT, xgrid, ygrid = crankNicolson(a, t, h, T, L)
print("Погрешность решения:")
deltaExplict = abs(analyticalT - u)
for i in range(len(xgrid)):
    for j in range(len(ygrid)):
        print("%-11.4f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-11.6f" % deltaExplict.max())
print()
fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Метод Кранка-Николсона")
axes.plot_surface(xgrid, ygrid, u, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("parabКранка-Николсона.png")


print("analytical")
for i in range(len(xgrid)):
    for j in range(len(ygrid)):
        print("%-11.4f" % (analyticalT[i, j]), end=' ')
    print()
print()
