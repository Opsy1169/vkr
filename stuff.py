import math as mp

from mpl_toolkits.mplot3d import Axes3D
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
import cmath
import scipy
import scipy.fftpack as sp
from scipy.integrate import quad
import mpl_toolkits.mplot3d as mpl
import pylab


def sinpow2():
    x = np.linspace(-2 * np.pi, 2 * np.pi, 100)
    y = np.sin(x) ** 2
    plt.plot(x, y)
    plt.show()


t = 5


def func(x):
    return cmath.exp(1j * (x ** 4) + t * 1j * (x ** 2))


Aiparam = 1;


def Ai(x):
    return cmath.exp(I * Aiparam * (x ** 3))

def AiI(x):
    return exp(I * Aiparam * (x ** 3) + I*x)


def plotSmthn():
    x = np.linspace(-0.5, 0.5, 700)
    # z = list(map(func, x))
    z = list(map(Ai, x))
    y = np.fft.fft(z);
    print('----------------------------------------------')
    yphasenp = list(map(cmath.phase, y))

    yabsnp = list(map(abs, y))
    z = np.fft.ifft(y)
    q = sp.fft(z)
    qphase = list(map(cmath.phase, q))
    qabs = list(map(abs, q))

    zphase = list(map(cmath.phase, z))
    # plt.plot(x, zphase)
    # plt.grid()
    # plt.show()
    #
    # plt.plot(x, yphasenp, x, yabsnp)
    # plt.grid()
    # plt.title('t = 1')
    # plt.legend([ 'phase', 'module'])
    # plt.show()

    plt.plot(x, qphase, x, qabs)
    plt.grid()
    plt.title('t = 1')
    plt.legend(['phase', 'module'])
    plt.show()


def Pearcey():
    a = 1
    b = 2
    Nx = 601
    Nu = 601

    x = np.linspace(-3, 3, Nx)
    F = [0] * Nu
    u = np.linspace(-3, 3, Nu)
    k = 2 * np.pi / (0.000633 * 1000)
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += func((x[i] + x[i + 1]) / 2) * cmath.exp((-I*k * (x[i] + x[i + 1]) * u[j]) / 2) * (x[i + 1] - x[i])
        F[j] = F[j] * (6) / (Nx)
    print("done")
    Fabs = list(map(abs, F))
    # another = list(map(np.real, F))
    # anotheranother =list(map(np.imag, F))

    # F2dAbs = []
    #
    # for i in range(Nu):
    #     row = [j * Fabs[i] for j in Fabs]
    #     F2dAbs.append(row)
    #
    # # print(F2dAbs[50][50])
    # arrayF2dAbs = np.array(F2dAbs)
    # print(arrayF2dAbs)
    # x, y = np.meshgrid(u, u)
    #
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, arrayF2dAbs, cmap=plt.cm.jet)
    # pylab.show()
    # print("another done")
    plt.plot(u, Fabs)
    plt.grid()
    plt.show()
    # plt.plot(u, another)
    # plt.grid()
    # plt.show()
    # plt.plot(u, anotheranother)
    # plt.grid()
    # plt.show()


def Airy():
    a = 3
    b = 1
    Nx = 401
    Nu =301

    x = np.linspace(-2, 2, Nx)
    F = [0] * Nu
    u = np.linspace(-0.1, 4, Nu)
    k = 2*np.pi/(0.000633*1000)
    print("x: ")
    # print(x)
    print("u: ")
    # print(u)
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Ai((x[i] + x[i + 1]) / 2) * cmath.exp((-I * k * ((x[i] + x[i + 1])) * u[j])/2) * (x[i + 1] - x[i])#сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
            # print(F[j])
        F[j] = F[j] * (2 * a)/(Nx)
        # print(F[j])

    # print(F)
    print("done")
    Fabs = list(map(abs, F))

    # F2dAbs = []
    #
    # for i in range(Nu):
    #     row = [j * Fabs[i] for j in Fabs]
    #     F2dAbs.append(row)
    #
    # # print(F2dAbs[50][50])
    # arrayF2dAbs = np.array(F2dAbs)
    # print(arrayF2dAbs)
    # x, y = np.meshgrid(u, u)
    #
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, arrayF2dAbs, cmap=plt.cm.jet)
    # pylab.show()
    # print("another done")
    plt.plot(u, Fabs)
    plt.grid()
    plt.show()

def InitialAiryPhaze():
    x = np.linspace(-2, 2, 400)
    airyinit = list(map(Ai, x))
    airyinitphase = list(map(np.real, airyinit))
    plt.plot(x, airyinitphase)
    plt.grid()
    plt.show()

def lists():
    a = [1, 2, 3]
    newa = []
    for j in range(3):
        row = [i * a[j] for i in a]
        print(row)
        newa.append(row)

    print(newa)

#-------------------------------------------------------------------------------------- trying to do the same thing with scipy.integrate
# def AiryFourier(x, s):
#     return exp( 1j* Aiparam * (x ** 3) + 1j* x * s)
# def func(s):
#     def realAi(x):
#         return scipy.real(AiryFourier(x, s))
#     return quad(AiryFourier, -1, 1, args=(s))[0]
#
#
# vec_func = np.vectorize(func)
# vec_func(np.arange(-0.1, 0.7, 200))
# plotSmthn()

Pearcey()
# Airy()
# InitialAiryPhaze()
# print(Ai(2))
