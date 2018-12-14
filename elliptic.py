u"""Решение эллептической задачи."""
from __future__ import print_function  # Форматный вывод
import numpy
import math
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt


def ellips(x, y, X, Y):
    u"""
    x, y - шаги, X, Y -границы.

    Возвращает КРС с начальными условиями, КРС с аналитическими значениями
    вектор временных слоев, вектор разбиений
    """
    def func1(x, y, X, Y):
        u"""Значение искомой функции на границах."""
        if x == 0:
            return y
        if x == X:
            return 1 + y
        if y == 0:
            return x
        if y == Y:
            return 1 + x
        return 10

    def analytical(x, y):
        u"""Аналитическое рещение."""
        return x + y
    n = math.ceil(X/x) + 1
    m = math.ceil(Y/y) + 1
    table = numpy.ones([n, m], dtype=float)
    analyticalTable = numpy.ones([n, m], dtype=float)
    xgird = numpy.ones(n)
    ygird = numpy.ones(m)
    for i in range(n):  # T
        for j in range(m):  # X
            table[i, j] = func1(i * x, j * y, X, Y)
            analyticalTable[i, j] = analytical(i * x, j * y)
            ygird[j] = j * y
        xgird[i] = i * x
    return table, analyticalTable, xgird, ygird


def libman(x, y, X, Y, eps):
    u"""Решение эллептического уравнения.

    Использован центрально-семмитричный шаблон и метод Либмана
    Принимает на вход значения a, сеточные шаги и границы условия
    """
    def interpolation(x1, x2, y1, y2, h):
        x = numpy.ones(math.ceil((x2 - x1) / h)+1)
        for i in range(len(x)):
            x[i] = y1 + (((x1 + h * i) - x1) / (x2 - x1) * (y2 - y1))
        return x
    u, analytical, xgird, ygird = ellips(x, y, X, Y)
    n = len(xgird)
    m = len(ygird)
    for i in range(1, n-1):
        u[i] = interpolation(xgird[0], xgird[n-1], u[i][0], u[i][m-1], y)
    for k in range(100):  # Максимальное количество итераций 1000
        newU = u.copy()
        maxE = 0.0
        for i in range(1, n-1):
            for j in range(1, m-1):
                newU[i, j] = 0.25 * (u[i+1, j] + u[i-1, j] +
                                     u[i, j+1] + u[i, j-1])
                maxE = max(maxE, abs(u[i, j] - newU[i, j]))
        u = newU
        if maxE < eps:
            pass
            break
    else:
        print("Заданная точность не достигнута")
        return
    for i in range(n):
        for j in range(m):
            print("%-6.2f" % (u[i, j]), end=' ')
        print()
    print()
    print("Analytic")
    for i in range(n):
        for j in range(m):
            print("%-6.2f" % (analytical[i, j]), end=' ')
        print()
    print()
    xgird, ygird = numpy.meshgrid(xgird, ygird)
    return u, analytical, xgird, ygird


u, anT, x, y = libman(0.1, 0.1, 1., 1., 0.00001)
deltaExplict = abs(anT - u)
for i in range(len(x)):
    for j in range(len(y)):
        print("%-6.2f" % (deltaExplict[i, j]), end=' ')
    print()
print("Максимальная погрешность: ", "%-6.2f" % deltaExplict.max())
fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Центрально-симметричный шаблон.")
axes.plot_surface(x, y, anT, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("ellipsЦСШ.png")

fig = pylab.figure()
axes = Axes3D(fig)
plt.title(u"Аналитическое решение.")
axes.plot_surface(x, y, anT, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
fig.savefig("ellipsАналит.png")
