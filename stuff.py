

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
import mpmath as mp


def sinpow2():
    x = np.linspace(-2 * np.pi, 2 * np.pi, 100)
    y = np.sin(x) ** 2
    plt.plot(x, y)
    plt.show()


t = 5

output = 0

def func(x):
    return cmath.exp(1j * (x ** 4) + t * 1j * (x ** 2))


Aiparam = 5


def Ai(x):
    return cmath.exp(I * Aiparam * (x ** 3))

def AiI(x):
    return exp(I * Aiparam * (x ** 3) + I*x)


def plotSmthn():
    x = np.linspace(-0.5, 2, 200)
    z = list(map(Ai, x))
    q = sp.fft(z)
    qshift = np.fft.fftshift(q)
    qphase = list(map(cmath.phase, q))
    qabs = list(map(abs, q))
    qabsshift = list(map(abs, qshift))
    zphase = list(map(cmath.phase, z))
    plt.plot( x, qabs, x, qabsshift )
    plt.grid()
    plt.title('t = 1')
    plt.legend(['phase', 'module'])
    plt.show()


def Pearcey():
    a = 1
    b = 2
    Nx = 201
    Nu = 201

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

    F2dAbs = []

    for i in range(Nu):
        row = [j * Fabs[i] for j in Fabs]
        F2dAbs.append(row)

    # print(F2dAbs[50][50])
    arrayF2dAbs = np.array(F2dAbs)
    # print(arrayF2dAbs)
    x, y = np.meshgrid(u, u)
    #
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(x, y, arrayF2dAbs, cmap=plt.cm.jet)
    pylab.show()
    # print("another done")
    # plt.plot(u, Fabs)
    # plt.grid()
    # plt.show()
    # plt.plot(u, another)
    # plt.grid()
    # plt.show()
    # plt.plot(u, anotheranother)
    # plt.grid()
    # plt.show()


def Airy():
    a = 3
    b = 1
    Nx = 201
    Nu =101

    x = np.linspace(-1, 3, Nx)
    F = [0] * Nu
    u = np.linspace(-0.1, 4, Nu)
    k = 2*np.pi/(0.000633*1000)

    # print(x)

    # print(u)
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Ai((x[i] + x[i + 1]) / 2) * cmath.exp((-I * k * ((x[i] + x[i + 1])) * u[j])/2) * (x[i + 1] - x[i])#сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
            # print(F[j])
        # F[j] = F[j] * (2 * a)/(Nx)
        # print(F[j])

    G = [0]*Nu
    for i in range(Nu):
        G[i] = Float(abs(mp.quad(lambda x: Ai(x) * cmath.exp(-I * k * x * u[i]), [-1, 3])), 15)
        print(G[i])
    print(F)
    # print(F)
    print("airy")
    print(k, F[0], u[0])
    Fabs = list(map(abs, F))
    print(Fabs)

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
    plt.plot(u, Fabs, u, G)
    plt.grid()
    plt.show()



def InitialAiryPhaze():
    x = np.linspace(-2, 2, 400)
    airyinit = list(map(Ai, x))
    airyinitphase = list(map(np.real, airyinit))
    plt.plot(x, airyinitphase)
    plt.grid()
    plt.show()


def quadint():
    Nx = 101
    Nu = 101
    x = np.linspace(-2, 2, Nx)
    Fr = [0] * Nu
    Fi = [0] * Nu
    F = [0] * Nu
    u = np.linspace(-0.1, 4, Nu)
    k = 2 * np.pi / (0.000633 * 1000)
    output = 0


    def funcforquad(x, u, k):
        return (Ai(x) * cmath.exp((I * k * (x) * u)))
    def funcforquadreal(x, u, k):
        return scipy.real(Ai(x) * cmath.exp((I * k * (x) * u)))
    def funcforquadim(x, u, k):
        return scipy.imag(Ai(x) * cmath.exp((I * k * (x) * u)))

    whole = funcforquad(x[90], u[90], 9)
    real = funcforquadreal(x[90], u[90], 9)
    imag = funcforquadim(x[90], u[90], 9)
    print(scipy.integrate.quad(funcforquadreal, 0, 1, args=(1, 9))[0])
    print(scipy.integrate.quad(funcforquadim, 0, 1, args=(1, 9))[0])




def anotherintegrate():
    Nx = 101
    Nu = 301
    x = np.linspace(-2, 2, Nx)
    Fr = [0] * Nu
    Fi = [0] * Nu
    F = [0] * Nu
    mp.dps = 30
    u = np.linspace(-0.1, 4, Nu)
    k = 2 * np.pi / (0.000633 * 1000)
    for i in range(Nu):
        F[i] = Float(abs(mp.quad(lambda x: Ai(x)*cmath.exp(-I*9.92*x*u[i]), [-1, 3])), 15)
        print(F[i])
    print(F)

    # Fabs = list(map(abs, F))
    plt.plot(u, F)
    plt.grid()
    plt.show()

plotSmthn()
# Pearcey()
# Airy()
# InitialAiryPhaze()
# print(Ai(2))
# sinpow2()
# quadint()
# a()
# anotherintegrate()
