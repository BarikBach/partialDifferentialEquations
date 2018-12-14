u"""Решение гиперболической задачи."""
from __future__ import print_function  # Форматный вывод
import numpy
import math
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt


def giperbol(a, t, h, l, T):
    u"""Создание начальной КРС, и аналитического решения.

    Вход: t,h - шаги, T,l -границы.
    Возвращает КРС с начальными условиями, КРС с аналитическими значениями
    вектор временных слоев, вектор разбиений
    """
    def func1(a, x, t, l, tau, h):
        u"""Значение искомой функции на границах."""
        if x == 0:
            return - math.sin(a * t)
        if x == l:
            return math.sin(a * t)
        if t == 0:
            return math.sin(x)
        if t == tau:
            return - a * tau * math.cos(x) + math.sin(x)
        return 10

    def analytical(t1, x1, a):
        u"""Аналитическое решение."""
        return math.sin(x1 - a * t1)
    # Получение количества ячеек по T
    n = math.ceil(T/t) + 1
    # Получение количества ячеек по X
    m = math.ceil(l/h) + 1
    table = numpy.ones([n, m], dtype=float)
    analyticalTable = numpy.ones([n, m], dtype=float)
    T = numpy.ones(n)
    X = numpy.ones(m)
    for i in range(n):  # T
        for j in range(m):  # X
            table[i, j] = func1(a, j * h, i * t, l, t, h)
            analyticalTable[i, j] = analytical(i * t, j * h, a)
            T[i] = i * t
        X[j] = j * h
    return table, analyticalTable, T, X


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
        alpha[i] = (- c[i] / (b[i] + a[i] * alpha[i-1]))
        beta[i] = ((d[i] - a[i] * beta[i-1]) / (b[i] + a[i] * alpha[i-1]))
    x = numpy.empty(n)
    x[n-1] = ((d[n-1]-a[n-1]*beta[n-2])/(b[n-1]+a[n-1]*alpha[n-2]))
    for i in range(n-2, -1, -1):
        x[i] = alpha[i] * x[i+1] + beta[i]
    return x


def explicit_finite_difference(a1, t, h, T, l):
    u"""Решение гиперболического уравнения с помощью явной схемы.

    Принимает на вход значения a, сеточные шаги и границы условия
    """
    u, analytical, xgird, ygird = giperbol(a1, t, h, l, T)
    n = len(xgird)
    m = len(ygird)
    sigma = float(a1**2 * t**2 / (h**2))
    if sigma >= 1:
        print("Схема неустойчива")
    print("explicit")
    for k in range(1, n-1):  # T
        for j in range(1, m-1):  # X
            u[k + 1, j] = (sigma * (u[k, j-1] - 2 * u[k, j] + u[k, j+1]) -
                           (u[k-1, j] - 2 * u[k, j]))
    for i in range(n):
        for j in range(m):
            print("%-9.4f" % (u[i, j]), end=" ")
        print()
    print()
    print("Analitycal")
    for i in range(n):
        for j in range(m):
            print("%-9.4f" % (analytical[i, j]), end=' ')
        print()
    xgird, ygird = numpy.meshgrid(xgird, ygird)
    return u, analytical, xgird, ygird


def implict_finite_difference(a1, t, h, T, l):
    u""" Решение гиперболического уравнения с помощью неявной схемы.

    Принимает на вход значения a, сеточные шаги и границы условия
    """
    u, analytical, x, y = giperbol(a1, t, h, l, T)
    print("implicit")
    n = len(x)
    m = len(y)
    a = numpy.zeros(n-1)
    b = numpy.zeros(n-1)
    c = numpy.zeros(n-1)
    d = numpy.zeros(n-1)
    sigma = float(a1**2 * t**2 / (h**2))
    # вынесено из цикла по слоям, т.к. во время итераций не изменяется
    for j in range(1, n-1):  # Кроме границ,
        # Создаем трехдиагональную матрицу для всех временных слоев
        if j > 1:
            a[j] = sigma
        if (j < (n-2)):
            c[j] = sigma
        b[j] = - (1 + 2 * sigma)
    for k in range(1, m-1):  # Временные слои
        for j in range(1, n-1):  # Кроме границ
            # Создаем трехдиагональную матрицу для всех временных слоев
            if j in range(1, n-2):
                d[j] = -2 * u[k][j] + u[k-1][j]
        d[1]   = (-2 * u[k][1]   + u[k-1][1])   - sigma * u[k+1][0]
        d[n-2] = (-2 * u[k][n-2] + u[k-1][n-2]) - sigma * u[k+1][n-1]
        x1 = tridiagonalMatrixAlgorithm(a, b, c, d)
        # Подставляем найденные задачи в КРС
        for j in range(1, n-1):
            u[k+1][j] = x1[j]
    for i in range(len(x)):
        for j in range(len(y)):
            print("%-9.4f" % (u[i, j]), end=' ')
        print()
    print()
    print("Analitycal")
    for i in range(n):
        for j in range(m):
            print("%-9.4f" % (analytical[i, j]), end=' ')
        print()
    x, y = numpy.meshgrid(x, y)
    return u, analytical, x, y


a = 0.025
t = 0.1
h = 0.1 * math.pi
T = 1.
L = math.pi
print("a =", a, "t =", t, "h =", h)

u, anT, x, y = explicit_finite_difference(a, t, h, T, L)
print("Погрешность решения:")
deltaExplict = abs(anT - u)
for i in range(len(x)):
    for j in range(len(y)):
        print("%-9.4f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-9.6f" % deltaExplict.max())
print()
fig = pylab.figure()
axes = Axes3D(fig)
axes.plot_surface(x, y, u, rstride=1, cstride=1, cmap=cm.jet)
plt.title(u"Явная схема")
pylab.show()
fig.savefig("gipЯвная")

u, anT, x, y = implict_finite_difference(a, t, h, T, L)
print("Погрешность решения:")
deltaExplict = abs(anT - u)
for i in range(len(x)):
    for j in range(len(y)):
        print("%-9.4f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-9.6f" % deltaExplict.max())
fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Неявная схема")
axes.plot_surface(x, y, u, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("gipНеявная")

fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Аналитическое решение")
axes.plot_surface(x, y, anT, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("gipАналит")
