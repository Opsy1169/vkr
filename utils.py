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
import fourier
import fresnel


def plotSmthn():
    lenx = 7000
    x = np.linspace(-0.5, 2, lenx).tolist()
    zer = np.zeros(3200).tolist()
    # zer = []
    x = np.array(zer + x + zer)
    z = list(map(beams.Pe, x))
    z = np.fft.fftshift(z)
    q = np.fft.fft(z)
    qshift = np.fft.fftshift(q)
    qphase = list(map(cmath.phase, q))
    qabs = list(map(abs, q))

    qabsshift = list(map(abs, qshift))
    qabsshift = qabsshift[len(zer):len(zer) + lenx]
    zphase = list(map(cmath.phase, z))
    plt.plot(x[len(zer):len(zer) + lenx], qabsshift)
    plt.grid()
    plt.title('t = 1')
    plt.legend(['phase', 'module'])
    plt.show()


def Pearcey():
    Nx = 5001
    Nu = 2001
    u_left = -3
    u_right = 3
    x_left = -7
    x_right = 7
    F = [0] * Nu
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    Fabs = fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, False, beams.Pe)
    maxVal = np.amax(Fabs)
    mean = reduce(lambda x1, x2: x1 + x2, Fabs) / len(Fabs)
    print('maxval:', str(maxVal))
    print('mean:', str(mean))


def halfPearcey():
    Nx = 8001
    Nu = 2001
    u_left = -4
    u_right = 5
    x_left = -5
    x_right = 5
    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    # k = 2 * np.pi / (0.000633)
    z = 0.95
    xo = 1
    zc = 1
    l = 0.6
    k = 2 * np.pi / (l)
    #
    x0 = l
    zc = (2 * k * l ** 2)
    z = (2 * k * l ** 2) * 0.9

    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += np.exp((1j * ((x[i] + x[i + 1]) / 2) ** 2 * u[j] / (x0 * (1 - z / zc) ** (1 / 2))) + 1j * (
                        (x[i] + x[i + 1]) / 2) ** 4) * ((x[i + 1]) - x[i])

    print("done")
    F = list(map(lambda x: x / ((1 - z / zc) ** (1 / 4)), F))
    Fabs = list(map(abs, F))
    plt.plot(u, Fabs)
    plt.grid()
    # plt.savefig("plots/pe_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(
    #     x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ").png")
    plt.show()


def Airy():
    Nx = 3501
    Nu = 2001
    u_left = -.1
    u_right = 15
    x_left = -4
    x_right = 4
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, False, beams.Ai)


def PearceyOdd():
    Nx = 2501
    Nu = 2001
    u_left = -.1
    u_right = 3
    x_left = -2
    x_right = 2
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, False, beams.PeOdd)


def InitialAiryPhaze():
    x = np.linspace(-2, 2, 400)
    airyinit = list(map(beams.Ai, x))
    airyinitphase = list(map(np.real, airyinit))
    plt.plot(x, airyinitphase)
    plt.grid()
    plt.show()



def acceleration():
    Nx = 3001
    Nu = 1001
    u_left = -6
    u_right = 6
    x_left = -2
    x_right = 2
    u = np.linspace(u_left, u_right, Nu)

    z_quan = 4
    repeat_time = 40
    z_0 = 100
    z_step = 10
    v = np.linspace(z_0, z_0 + z_step * (z_quan - 1), z_quan)
    output = []
    start = time.time()
    for i in range(z_quan):
        outputLine = fresnel.fresnel(i * z_step + z_0, Nx, Nu, x_right, x_left, u_right, u_left, beams.Ai)
        output.append(outputLine)
        # for j in range(repeat_time):
        #     output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end - start))
    u, v = np.meshgrid(u, v)
    # fig = plt.figure()
    # axes = Axes3D(fig)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.view_init(100, 100)
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')
    # axes = fig.add_subplot(122, projection='3d')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary)
    pylab.show()



def accelerationVer3():
    Nx = 3001
    Nu = 5001
    u_left = -1
    u_right = 9
    x_left = -4
    x_right = 4
    u = np.linspace(u_left, u_right, Nu)


    z_quan = 1
    z_0 = 400
    z_step = 40
    v = np.linspace(z_0, z_0 + z_step * (z_quan - 1), z_quan)


    out_u = np.linspace(u_left, u_right,  Nu)
    output = []
    start = time.time()
    func = beams.Ai
    new_u = []
    for i in range(Nu - 1):
        new_u.append((u[i + 1] + u[i]) / 2)
    Ffour = fourier.fourierArr(Nx, x_right, x_left, 300, new_u, func)
    Fabs = list(map(lambda x: abs(x), Ffour))
    plt.plot(new_u, Fabs)
    plt.xlabel('x, mm')
    plt.ylabel('y, mm')
    plt.grid()
    plt.show()


    for i in range(z_quan):
        outputLine = fresnel.fresnelArr(i * z_step + z_0, new_u, out_u,  Ffour)
        output.append(outputLine)
        # for j in range(repeat_time):
        #     output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end - start))
    out_u, v = np.meshgrid(out_u, v)
    # fig = plt.figure()
    # axes = Axes3D(fig)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.view_init(100, 100)
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')
    # axes = fig.add_subplot(122, projection='3d')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary)
    pylab.show()

def fresfour(z, Nx, Nu, x_right, x_left, u_right, u_left, func, f, is2d):
    # x = np.linspace(x_left, x_right, Nx)
    F = []
    u = np.linspace(u_left, u_right, Nu)
    out_u = np.linspace(1.1 * u_left, 1.1 * u_right, int(0.5 * Nu))
    k = 2 * np.pi / (0.000633)
    # z = 800
    temp = (1j * k / (2 * z))
    start = time.clock()
    # Ffour = [0] * (Nu)
    new_u = []
    for i in range(Nu - 1):
        new_u.append((u[i + 1] + u[i]) / 2)
    Ffour = fourier.fourierArr(Nx, x_right, x_left, f, new_u, func)
    # for i in range(Nu - 1):
    #     Ffour[i] = fourierDot(Nx, x_right, x_left, f, (u[i + 1] + u[i]) / 2, func)
    print("done args")
    end = time.clock()
    print('time:' + str(end - start))
    # if(i < 6000):
    #     Fmid[i] = 0
    Fabs = list(map(lambda x: abs(x), Ffour))
    plt.plot(new_u, Fabs)
    plt.xlabel('x, mm')
    plt.ylabel('y, mm')
    plt.grid()
    plt.show()

    for n in range(len(out_u)):
        buff = 0
        for i in range(Nu - 1):
            buff += Ffour[i] * np.exp(temp * (new_u[i] - out_u[n]) ** 2)
        F.append(buff)
    F = list(map(lambda x: x * np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z)), F))
    Fabs = list(map(lambda x: abs(x), F))
    end = time.clock()

    print('time:' + str(end - start))
    if (is2d):
        F2dAbs = []
        for i in range(Nu):
            row = [j * Fabs[i] for j in Fabs]
            F2dAbs.append(row)

        arrayF2dAbs = np.array(F2dAbs)
        z, y = np.meshgrid(u, u)

        fig = pylab.figure()
        ax = fig.gca(projection='3d', proj_type='ortho')
        # axes = Axes3D(fig)
        ax.plot_surface(z, y, arrayF2dAbs, cmap=plt.cm.jet)
        # pylab.savefig("plots/pe3d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param( " + str(pe_param) + ").png")
        pylab.show()
        print("another done")
    else:
        plt.plot(out_u, Fabs)
        plt.xlabel('x, mm')
        plt.ylabel('y, mm')
        plt.grid()
        # plt.savefig(
        #     "plots/pe_odd_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(
        #         x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ").png")
        plt.show()



def accelerationVer2():
    Nx = 4001
    Nu = 1001
    u_left = -3
    u_right = 3
    x_left = -2
    x_right = 2
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    f = 300
    z_quan = 2
    z_0 = 200
    z_step = 200
    v = np.linspace(z_0, z_0 + z_step * (z_quan - 1), z_quan)
    output = []
    start = time.time()
    Fmid = [0] * (Nx - 1)
    func = beams.Ai
    new_x = []
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
        new_x.append((x[i + 1] + x[i]) / 2)
    for i in range(z_quan):
        outputLine = fresnel_four(f, i * z_step + z_0, Nx, Nu, new_x, u, Fmid)
        output.append(outputLine)
        # for j in range(repeat_time):
        #     output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end - start))
    u, v = np.meshgrid(u, v)
    # fig = plt.figure()
    # axes = Axes3D(fig)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.view_init(100, 100)
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')
    # axes = fig.add_subplot(122, projection='3d')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary)
    pylab.show()

def fresnel_four(f, z, Nx, Nu, x, u, Fmid):

    F = [0] * Nu

    k = 2 * np.pi / (0.000633)
    start = time.clock()
    # new_x = []
    # for i in range(Nx - 1):
    #     new_x.append((x[i + 1] + x[i]) / 2)
    step_x = x[1] - x[0]
    temp = (-1j * k / (z))
    dx = (x[-1] - x[0])/Nx
    secondtemp = 1j * (k / 2) * ((f - z) / (f * z))
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * np.exp(temp * (x[i]) * u[n] + secondtemp * (x[i]) ** 2)  # сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
        F[n] *= np.exp(-temp/2 * u[n]**2)*dx
        # F[n] *=
        # print(F[j])
    F = list(map(lambda x: x * np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z)), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print('time:' + str(end - start))

    plt.plot(u, Fabs)
    # plt.show()
    return Fabs
    # plt.plot(u, Fabs)
    # plt.grid()
    # plt.savefig(
    #     "plots/fresnel_four_ai_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(
    #         x_left) + ", " + str(
    #         x_right) + ", " + str(Nx) + ")_param(" + str() + ")_z(" + str(z) + ")_f(" + str(f) + ").png")
    # plt.show()


def normalizeArray(arr):
    max = np.max(arr)
    arr /= max
    return arr


def beamsForDifferentFocus(f, func):
    Nx = 2001
    Nu = 2001
    u_left = -1
    u_right = 10
    x_left = -2
    x_right = 2
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    for i in range(len(f)):
        a = fourier.fourierArr(Nx, x_right, x_left, f[i], u, func)
        # plt.plot(a, u)
        leg.append('f = ' + str(f[i]) + ', mm')
    plt.legend(leg)
    plt.show()

def beamsForDifferentParam(params, func):
    Nx = 3001
    Nu = 1501
    u_left = -1
    u_right = 10
    x_left = -1
    x_right = 1
    f = 200
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    print(func == beams.Ai)
    isAi = True if(func == beams.Ai) else False
    print(isAi)
    for i in range(len(params)):
        if(isAi):
            beams.setAiParam(params[i])
            print(beams.getAiParam())
        else:
            beams.setPeParam(params[i])
            print(beams.getPeParam())
        a = fourier.fourierArr(Nx, x_right, x_left, f, u, func)
        # plt.plot(a, u)
        leg.append('alpha = ' + str(params[i]))
    plt.legend(leg)
    plt.show()

def beamsForDifferentInputRange(rights, lefts, func):
    Nx = 2001
    Nu = 2501
    u_left = -1
    u_right = 10
    x_left = -2
    x_right = 2
    f = 200
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    for i in range(len(rights)):
        a = fourier.fourierArr(int((rights[i]-lefts[i])*900+1), rights[i], lefts[i], f, u, func)
        # plt.plot(a, u)
        leg.append('right = ' + str(rights[i]) + "; left = " + str(lefts[i]))
    plt.legend(leg)
    plt.show()

def accelerationCompare():
    Nx = 2001
    Nu = 1001
    u_left = -3
    u_right = 3
    x_left = -1
    x_right = 1
    f = 150
    z = 250
    param = 20
    beams.setPeParam(param)
    beams.setAiParam(param)
    x = np.linspace(x_left, x_right, Nx)
    u = np.linspace(u_left, u_right, Nu)
    Ai_input = [0] * (Nx - 1)
    Pe_input = [0] * (Nx - 1)
    new_x = []
    for i in range(Nx - 1):
        Ai_input[i] = beams.Ai((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
        Pe_input[i] = beams.Pe((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
        new_x.append((x[i + 1] + x[i]) / 2)
    Ai_output = fresnel_four(f, z, Nx, Nu, new_x, u, Ai_input)
    Pe_output = fresnel_four(f, z, Nx, Nu, new_x, u, Pe_input)
    plt.plot(u, Ai_output, u, Pe_output)
    plt.show()


def pearceyShit():
    Nx = 1001
    Nu = 1001
    u_left = -1
    u_right = 1
    x_left = -1
    x_right = 1
    param_start = 0
    param_end = 3
    N_param = 30
    start = time.clock()
    params = np.linspace(param_start, param_end, N_param)
    u = np.linspace(u_left, u_right, Nu)
    new_u = np.array(list(map(lambda x: x/0.01, u)))
    output = []
    for i in range(len(params)):
        beams.setPeParam(params[i]/0.1)
        temp = fourier.fourierArr(Nx, x_right, x_left, 200, new_u, beams.Pe)
        temp = list(map(abs, temp))
        output.append(temp)

    end = time.clock()
    output = np.array(output)
    print('time:' + str(end - start))
    u, params = np.meshgrid(u, params)
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(u, params, output, cmap=plt.cm.jet)
    pylab.show()

