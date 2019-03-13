

from mpl_toolkits.mplot3d import Axes3D
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
import cmath
from scipy import special
import scipy.fftpack as sp
from scipy.integrate import quad
import mpl_toolkits.mplot3d as mpl
import pylab
import mpmath as mp
import time


def sinpow2():
    x = np.linspace(-2 * np.pi, 2 * np.pi, 100)
    y = np.sin(x) ** 2
    plt.plot(x, y)
    plt.show()


pe_param = 5

output = 0

def Pe(x):
    return cmath.exp(1j * (x ** 4) + pe_param * 1j * (x ** 2))


Aiparam = 10

sigma = 0.1
def Gauss(x):
    return cmath.exp(-x**2/sigma**2)

def Ai(x):
    return cmath.exp(I * Aiparam * (x ** 3))

def AiI(x):
    return exp(I * Aiparam * (x ** 3) + I*x)


def plotSmthn():
    lenx = 7000
    x = np.linspace(-0.5, 2, lenx).tolist()
    zer = np.zeros(3200).tolist()
    # zer = []
    x = np.array(zer + x + zer)
    z = list(map(Pe, x))
    z = np.fft.fftshift(z)
    q = np.fft.fft(z)
    qshift = np.fft.fftshift(q)
    qphase = list(map(cmath.phase, q))
    qabs = list(map(abs, q))

    qabsshift = list(map(abs, qshift))
    qabsshift = qabsshift[len(zer):len(zer)+lenx]
    zphase = list(map(cmath.phase, z))
    plt.plot(x[len(zer):len(zer)+lenx], qabsshift )
    plt.grid()
    plt.title('t = 1')
    plt.legend(['phase', 'module'])
    plt.show()


def Pearcey():
    a = 1
    b = 2
    Nx = 801
    Nu = 401
    u_left = -4
    u_right = 4
    x_left = -4
    x_right = 4
    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633 * 1300)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = Pe((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Fmid[i] * cmath.exp((-I*k * (x[i] + x[i + 1]) * u[j]) / 2)
        # F[j] = F[j] * (x_right-x_left) / (Nx)
    print("done")
    F = list(map(lambda x: x*(x_right-x_left) / (Nx), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print("time:" + str(end-start))
    # another = list(map(np.real, F))
    # anotheranother =list(map(np.imag, F))

    # F2dAbs = []
    #
    # for i in range(Nu):
    #     row = [j * Fabs[i] for j in Fabs]
    #     F2dAbs.append(row)

    # print(F2dAbs[50][50])
    # arrayF2dAbs = np.array(F2dAbs)
    # print(arrayF2dAbs)
    # x, y = np.meshgrid(u, u)
    #
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, arrayF2dAbs, cmap=plt.cm.jet)
    # pylab.savefig("plots/pe3d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param( " + str(pe_param) + ").png")
    # pylab.show()
    # print("another done")
    plt.plot(u, Fabs)
    plt.grid()
    plt.savefig("plots/pe_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ").png")
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
    Nu = 401
    u_left = -0.1
    u_right = 7
    x_left = -2
    x_right = 2

    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
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
    # mp.dps = 50
    # for i in range(Nu):
    #     G[i] = Float(abs(mp.quad(lambda x: Ai(x) * cmath.exp(-I * k * x * u[i]), [-1, 3])), 35)
    #     print(G[i])
    print(F)
    # print(F)
    print("airy")
    print(k, F[0], u[0])
    Fabs = list(map(abs, F))
    print(Fabs)

    F2dAbs = []

    for i in range(Nu):
        row = [j * Fabs[i] for j in Fabs]
        F2dAbs.append(row)

    # print(F2dAbs[50][50])
    arrayF2dAbs = np.array(F2dAbs)
    print(arrayF2dAbs)
    x, y = np.meshgrid(u, u)

    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(x, y, arrayF2dAbs, cmap=plt.cm.jet)
    pylab.savefig("plots/ai2d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(Aiparam) + ").png")
    pylab.show()
    # print("another done")
    # plt.plot(u, Fabs, u, G)
    # plt.grid()
    # plt.legend(['Numeric', 'mp.quad'])
    # plt.savefig("plots/ai_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ").png")
    # plt.show()

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


    def Peforquad(x, u, k):
        return (Ai(x) * cmath.exp((I * k * (x) * u)))
    def Peforquadreal(x, u, k):
        return scipy.real(Ai(x) * cmath.exp((I * k * (x) * u)))
    def Peforquadim(x, u, k):
        return scipy.imag(Ai(x) * cmath.exp((I * k * (x) * u)))

    whole = Peforquad(x[90], u[90], 9)
    real = Peforquadreal(x[90], u[90], 9)
    imag = Peforquadim(x[90], u[90], 9)
    print(scipy.integrate.quad(Peforquadreal, 0, 1, args=(1, 9))[0])
    print(scipy.integrate.quad(Peforquadim, 0, 1, args=(1, 9))[0])

def anotherintegrate():
    Nx = 101
    Nu = 501
    x = np.linspace(-2, 2, Nx)
    Fr = [0] * Nu
    Fi = [0] * Nu
    F = [0] * Nu
    mp.dps = 70
    u = np.linspace(-0.1, 4, Nu)
    k = 2 * np.pi / (0.000633 * 1000)
    for i in range(Nu):
        F[i] = Float(abs(mp.quad(lambda x: Ai(x)*cmath.exp(-I*9.92*x*u[i]), [-1, 3])), 50)
        print(F[i])
    print(F)

    # Fabs = list(map(abs, F))
    plt.plot(u, F)
    plt.grid()
    plt.show()


def fresnel():
    a = 3
    b = 1
    Nx = 801
    Nu = 401
    u_left = -4
    u_right = 4
    x_left = -4
    x_right = 4

    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633)
    z = 650

    # print(x)

    # print(u)
    a = (1.777514461805397-1.777514461805397j)
    start = time.clock()
    Fmid = [0]*(Nx-1)
    for i in range(Nx - 1):
        Fmid[i] = Gauss((x[i] + x[i + 1]) / 2)* (x[i + 1] - x[i])
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * cmath.exp((I * k /(2*z))* ((x[i] + x[i + 1])/2 - u[n])**2)   # сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье

        # F[n] *=
            # print(F[j])
    F = list(map(lambda x: x*cmath.exp(1j * k * z) * np.sqrt(1j * k / (2 * np.pi * z)), F))
    print(np.sqrt(-1j * k / (2 * np.pi * z)))
    print('a')
    # F = np.fft.fftshift(F)
    # G = [0] * Nu
    # mp.dps = 50
    # for i in range(Nu):
    #     G[i] = Float(abs(mp.quad(lambda x: Ai(x) * cmath.exp(-I * k * x * u[i]), [-1, 3])), 35)
    #     print(G[i])

    # print(F)
    # print("airy")
    # print(k, F[0], u[0])
    Fabs = list(map(abs, F))
    end = time.clock()
    oldtime = "time:193.04681520659014"
    secondtime = "time:97.04681520659014"
    thirdtime = "105"
    print('time:' + str(end-start))
    # print(Fabs)
    plt.plot(u, Fabs)
    plt.grid()
    plt.savefig("plots/fresnelGauss_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(Aiparam) + ")_z(" + str(z) + ").png")
    plt.show()

def scipyfresnel():
    x = np.linspace(-30, 5, 1001)
    # x1 = list(map(abs, x))
    z = special.airy(x)[0]
    # my = list(map(Ai, x))
    zabs = list(map(abs, z))
    zphase = list(map(cmath.phase, z))
    # zphasemy = list(map(cmath.phase, my))
    # f = special.fresnel(z)
    # f = list(map(abs, f))
    plt.plot(x, zabs, x, zphase)
    plt.grid()
    plt.show()

# plotSmthn()
# Pearcey()
# Airy()
# InitialAiryPhaze()
# print(Ai(2))
# sinpow2()
# quadint()
# a()
# anotherintegrate()
fresnel()
# scipyfresnel()
