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
from functools import reduce

import beams


def fourierArr(Nx, x_right, x_left, f, u, func):
    x = np.linspace(x_left, x_right, Nx)
    # u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633 * f)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    output = []
    temp = -1j*k/2
    # temp = 1j
    # multiplier = (x_right - x_left) / (Nx)
    multiplier = 1
    for j in range(len(u)):
        val = 0
        u_j = u[j]
        for i in range(Nx - 1):
            val += Fmid[i] * np.exp((temp * (x[i] + x[i + 1]) * u_j))
            # F[j] = F[j] * (x_right-x_left) / (Nx)
        output.append(val * multiplier)
    Fabs = list(map(abs, output))
    plt.plot(u, Fabs)
    return np.array(output)


def fourierDot(Nx, x_right, x_left, f, u, func):
    x = np.linspace(x_left, x_right, Nx)
    # u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633 * f)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    output = 0
    temp = -1j * k / 2
    for i in range(Nx - 1):
        output += Fmid[i] * np.exp( temp* (x[i] + x[i + 1]) * u)
        # F[j] = F[j] * (x_right-x_left) / (Nx)
    # output = output * (x_right - x_left) / (Nx)
    return output


def Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, is2d, func):
    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633 * f)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])

    temp = -1j * k / 2
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Fmid[i] * np.exp((temp * (x[i] + x[i + 1]) * u[j]))
        # F[j] = F[j] * (x_right-x_left) / (Nx)
    # F = list(map(lambda x: x * (x_right - x_left) / (Nx), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print("time:" + str(end - start))
    if (is2d):
        F2dAbs = []
        for i in range(Nu):
            row = [j * Fabs[i] for j in Fabs]
            F2dAbs.append(row)

        arrayF2dAbs = np.array(F2dAbs)
        z, y = np.meshgrid(u, u)

        fig = pylab.figure()
        axes = Axes3D(fig)
        axes.plot_surface(z, y, arrayF2dAbs, cmap=plt.cm.jet)
        # pylab.savefig("plots/pe3d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param( " + str(pe_param) + ").png")
        pylab.show()
        print("another done")
    else:
        plt.plot(u, Fabs)
        plt.xlabel('x, mm')
        plt.ylabel('y, mm')
        plt.grid()
        plt.savefig(
            "plots/pe_odd_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(
                x_left) + ", " + str(
                x_right) + ", " + str(Nx) + ")_param(" + str() + ").png")
        plt.show()
    return Fabs


def Fourier2d():
    init, x, y = beams.getInitPe2d()
    start = time.clock()
    # init, x, y = getInitAi2d()
    uleft = -0.5
    uright = 0.5
    vleft = -0.5
    vright = 0.5

    nv = 71
    nu = 71
    koef = 2 * np.pi / (0.000633 * 300)
    u = np.linspace(uleft, uright, nu)
    v = np.linspace(vleft, vright, nv)
    output = []
    len_x = len(init[0])
    len_y = len(init)
    koef_mul_i = -1j * koef
    for i in range(nu):
        outputline = []
        u_i = u[i]
        for j in range(nv):
            integrateelem = 0
            v_j = v[j]
            for k in range(len_x):  # x
                x_k = x[k]
                x_u = x_k * u_i
                for p in range(len_y):  # y
                    integrateelem += init[k][p] * np.exp(koef_mul_i * (x_u + y[p] * v_j))
            outputline.append(integrateelem)
        output.append(outputline)
    # output = np.array(output)
    for i in range(nu):
        for j in range(nv):
            a = output[i][j]
            output[i][j] = abs(a)

    end = time.clock()
    print("TIME: ", end - start)
    output = np.array(output)
    u, v = np.meshgrid(u, v)
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(u, v, output, cmap=plt.cm.jet)
    pylab.show()


