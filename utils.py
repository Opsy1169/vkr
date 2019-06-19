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
    Nx = 4001
    Nu = 701
    u_left = -4
    u_right = 4
    x_left = -3
    x_right = 3
    F = [0] * Nu
    f = 200
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


def sinc():
    Nx = 1001
    Nu = 501
    u_left = -3
    u_right = 3
    x_left = -1
    x_right = 1
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, False, beams.Rect)

def gauss():
    Nx = 301
    Nu = 1001
    u_left = -3
    u_right = 3
    x_left = -3
    x_right = 3
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    # x = np.linspace(x_left, x_right, Nx)
    # Fmid = [0]*(Nx)
    # for i in range(Nx):
    #     Fmid[i] = beams.Gauss(x[i])
    # F = np.fft.fft(Fmid)
    # F = np.fft.fftshift(F)
    # plt.plot(x, F)
    # plt.show()
    fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, f, False, beams.Gauss1)

def gauss1():
    x = np.linspace(-5, 5, 200)
    a = [1, 2, 3, 4]
    print(np.fft.fftshift(a))
    f = list(map(lambda a: np.exp(-a**2), x))
    F = np.fft.fft(np.fft.fftshift(f))
    F = np.fft.fftshift(F)
    plt.plot(x, F)
    plt.grid()
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
    x_left = -1.5
    x_right = 1.5
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    f = 200
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

def fresnel_four_comparing(f, z, Nx, x_left, x_right, Nu, u, func):
    F = [0] * Nu
    k = 2 * np.pi / (0.000633)
    start = time.clock()
    x = np.linspace(x_left, x_right, Nx)
    step_x = x[1] - x[0]
    Fmid = [0] * (Nx - 1)
    new_x = []
    max_x = []
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
        new_x.append((x[i + 1] + x[i]) / 2)
    x = new_x
    temp = (-1j * k / (z))
    dx = (x[-1] - x[0])/Nx
    secondtemp = 1j * (k / 2) * ((f - z) / (f * z))
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * np.exp(temp * (x[i]) * u[n] + secondtemp * (x[i]) ** 2)  # сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
        F[n] *= np.exp(-temp/2 * u[n]**2)*dx
    F = list(map(lambda x: x * np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z)), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print('time:' + str(end - start))
    plt.plot(u, Fabs)
    return Fabs



def beamsForDifferentFocus(f, func):
    Nx = 2001
    Nu = 2001
    u_left = -3
    u_right = 5
    x_left = -1
    x_right = 1
    z = 400
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    max_val = []
    max_x = []
    Fmid = [0] * (Nx - 1)
    new_x = []
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
        new_x.append((x[i + 1] + x[i]) / 2)
    for i in range(len(f)):
        # a = fourier.fourierArr(Nx, x_right, x_left, f[i], u, func)
        a = fresnel_four(f[i], z, Nx, Nu, new_x, u, Fmid)
        a = np.array(a)
        # plt.plot(a, u)
        leg.append('f = ' + str(f[i]) + ', mm')
        # max_val.append(max(a))
        max_x.append(u[a.argmax()])
    print(max_x)
    plt.legend(leg)
    plt.show()

    plt.plot(f, max_x)
    plt.show()

def beamsForDifferentParam(params, func):
    Nx = 2001
    Nu = 1001
    u_left = -1
    u_right = 6
    x_left = -2
    x_right = 2
    f = 1000
    z = 400
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    Fmid = [0] * (Nx - 1)
    new_x = []
    max_x = []
    print(func == beams.Ai)
    isAi = True if(func == beams.Ai) else False
    print(isAi)
    for i in range(len(params)):
        if(isAi):
            beams.setAiParam(params[i])
        else:
            beams.setPeParam(params[i])
        a = fourier.fourierArr(Nx, x_right, x_left, f, u, func)
        # print(beams.getAiParam())
        # a = fresnel_four_comparing(f, z, Nx, x_left, x_right, Nu,  u, beams.Ai)
        a = np.array(a)
        max_x.append(u[a.argmax()])
        leg.append('alpha = ' + str(params[i]))
    plt.rcParams.update({'font.size': 12})
    plt.xlabel('х, мм')
    plt.ylabel('амплитуда')
    plt.legend(leg)
    plt.grid()
    print(max_x)
    plt.show()

def beamsForDifferentInputRange(rights, lefts, func):
    Nx = 1001
    Nu = 1501
    u_left = -1
    u_right = 6
    x_left = -2
    x_right = 2
    f = 1000
    z = 400
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    Fmid = [0] * (Nx - 1)
    new_x = []
    max_x = []
    for i in range(len(rights)):
        a = fourier.fourierArr(int((rights[i]-lefts[i])*900+1), rights[i], lefts[i], f, u, func)
        # b = int((rights[i]-lefts[i])*900+1)
        # a = fresnel_four_comparing(f, z,  b,  lefts[i], rights[i], Nu, u, beams.Ai)
        a = np.array(a)
        max_x.append(u[a.argmax()])
        # plt.plot(a, u)
        leg.append( "left = " + str(lefts[i])  + '; right = ' + str(rights[i]))
    plt.rcParams.update({'font.size': 12})
    plt.xlabel('х, мм')
    plt.ylabel('амплитуда')
    plt.legend(leg)
    plt.grid()
    print(max_x)
    plt.show()

def accelerationCompare():
    Nx = 2001
    Nu = 1001
    u_left = -3
    u_right = 3
    x_left = -1
    x_right = 1
    f = 300
    z = 500
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
        Pe_input[i] = beams.PeOdd((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
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

def plotshit():
    a = [0.02, 0.1, 0.15, 0.16, 0.13, 0.06, 0.02, -0.01]
    b = [-0.11, 0, -0.02, -0.17, -0.3, -0.45, -0.6, -1]
    c = [100, 150, 200, 250, 300, 350, 400, 500]
    c1 =[100, 150, 200, 250, 300]
    plt.rcParams.update({'font.size': 14})
    ##################################
    # f = [100, 150, 200, 250, 300, 350, 400, 450]
    # f_2_20 = [0.88, 1.32, 1.75, 2.2, 2.6, 3.12, 3.55, 4.01]
    # f_2_40 = [1.6, 2.43, 3.24, 4.1, 4.91, 5.72, 6.54, 7.4 ]
    # f_1_20 = [0.22, 0.36, 0.47, 0.58, 0.69, 0.81, 0.93, 1.05]
    # f_1_40 = [0.37, 0.58, 0.78, 0.99, 1.15, 1.37, 1.57, 1.77]
    # fig = plt.figure()
    # ax = plt.subplot(111)
    # ax.plot(f, f_2_20, label='alpha=20,x=2')
    # ax.plot(f, f_2_40, label='alpha=40,x=2')
    # ax.plot(f, f_1_20, label='alpha=20,x=1')
    # ax.plot(f, f_1_40, label='alpha=40,x=1')
    # ax.grid()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    #
    # # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #           fancybox=True, shadow=True, ncol=5)
    # ax.set_xlabel('Фокусное расстояние, мм')
    # ax.set_ylabel('Ширина, мм')
    ##############################################################

    # alpha = [10, 20, 30, 40, 50]
    # alpha_200_2 = [0.82, 1.5, 2.25, 2.9, 3.65]
    # alpha_400_2 = [1.67, 3.05, 4.5, 5.88, 7.32]
    # alpha_200_1 = [0.24, 0.28, 0.38, 0.5, 0.58]
    # alpha_400_1 = [0.3, 0.54, 0.78, 1.02, 1.25]
    #
    # ax = plt.subplot(111)
    # ax.plot(alpha, alpha_200_2, label='f=200,x=2')
    # ax.plot(alpha, alpha_400_2, label = 'f=400,x=2')
    # ax.plot(alpha, alpha_200_1, label='f=200,x=1')
    # ax.plot(alpha,  alpha_400_1, label='f=400,x=1')
    # ax.grid()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    #
    # # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #           fancybox=True, shadow=True, ncol=5)
    # ax.set_xlabel(u'Параметр \u03B1')
    # ax.set_ylabel('Ширина, мм')
    # plt.show()
    ##############################################################

    # input = [1, 2, 3, 4]
    # input_200_20 = [0.47, 1.58, 3.56, 6.28]
    # input_400_20 = [0.94, 3.15, 7.06, 12.8]
    # input_200_40 = [0.79, 3.06, 6.93, 12.45]
    # input_400_40 = [1.52, 6, 13.8, 24.8]
    #
    # ax = plt.subplot(111)
    # ax.plot(input, input_200_20, label='f=200,\u03B1=20')
    # ax.plot(input, input_400_20, label = 'f=400,\u03B1=20')
    # ax.plot(input, input_200_40, label='f=200,\u03B1=40')
    # ax.plot(input,  input_400_40, label='f=400,\u03B1=40')
    # ax.grid()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    #
    # # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #           fancybox=True, shadow=True, ncol=5)
    # ax.set_xlabel('Правая граница входного размера, мм')
    # ax.set_ylabel('Ширина, мм')
    # plt.show()
    ##############################################################

    #  ACCELERATION
    ###########################################################
    # f = [100, 150, 200, 250, 300, 350 ]
    # f_2_20 = [-2.66, -0.7399999999999998, -0.20399999999999974, -0.008000000000000007, 0.08400000000000007, 0.1120000000000001]
    # f_2_40 = [-1.272, -0.2919999999999998, -0.016000000000000014, 0.08400000000000007, 0.13200000000000012, 0.1280000000000001]
    # f_1_20 = [-1.64, -0.6320000000000001, -0.23599999999999977, -0.02400000000000002, 0.07600000000000007, 0.1080000000000001]
    # f_1_40 = [-1.22, -0.31999999999999984, -0.016000000000000014, 0.08400000000000007, 0.1120000000000001, 0.13600000000000012]
    #
    # ax = plt.subplot(111)
    # ax.plot(f, f_2_20, label='alpha=20,x=2')
    # ax.plot(f, f_2_40, label='alpha=40,x=2')
    # ax.plot(f, f_1_20, label='alpha=20,x=1')
    # ax.plot(f, f_1_40, label='alpha=40,x=1')
    # ax.grid()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    #
    # # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #           fancybox=True, shadow=True, ncol=5)
    # ax.set_xlabel('Фокусное расстояние, мм')
    # ax.set_ylabel('Отклонение, мм')
    ###################################################
    # alpha = [10, 20, 30, 40, 50]
    # alpha_150_1 =  [-0.19799999999999995, -0.32399999999999984, -0.48599999999999977, -0.6299999999999999, -0.8820000000000001]
    # alpha_300_1 = [0.15600000000000014, 0.11399999999999988, 0.1259999999999999, 0.07200000000000006, -0.01200000000000001]
    # alpha_150_2 = [-0.19199999999999973, -0.29400000000000004, -0.45599999999999996, -0.738, -1.554]
    # alpha_300_2 = [0.1499999999999999, 0.13200000000000012, 0.1080000000000001, 0.08400000000000007, 0.03000000000000025]
    # alpha_150_1.reverse()
    # alpha_300_1.reverse()
    # alpha_150_2.reverse()
    # alpha_300_2.reverse()
    #
    # ax = plt.subplot(111)
    # ax.plot(alpha, alpha_150_1, label='f=150,x=1')
    # ax.plot(alpha, alpha_300_1, label = 'f=300,x=1')
    # ax.plot(alpha, alpha_150_2, label='f=150,x=2')
    # ax.plot(alpha,  alpha_300_2, label='f=300,x=2')
    # ax.grid()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    #
    # # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #           fancybox=True, shadow=True, ncol=5)
    # ax.set_xlabel(u'Параметр \u03B1')
    # ax.set_ylabel('Отклонение, мм')
    # plt.show()
    ##################################################
    #
    input = [1, 2, 3, 4]
    input_150_20 = [-0.637, -0.736, -0.758, -0.747]
    input_300_20 = [0.07799999999999985, 0.07799999999999985, 0.07799999999999985, 0.07799999999999985]
    input_150_40 = [-0.31800000000000006, -0.29600000000000004, -0.29600000000000004, -0.29600000000000004]
    input_300_40 = [0.11099999999999999, 0.133, 0.12199999999999989, 0.12199999999999989]

    ax = plt.subplot(111)
    ax.plot(input, input_150_20, label='f=150,\u03B1=20')
    ax.plot(input, input_300_20, label = 'f=300,\u03B1=20')
    ax.plot(input, input_150_40, label='f=150,\u03B1=40')
    ax.plot(input,  input_300_40, label='f=300,\u03B1=40')
    ax.grid()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              fancybox=True, shadow=True, ncol=5)
    ax.set_xlabel('Правая граница входного размера, мм')
    ax.set_ylabel('Отклонение, мм')
    plt.show()

def scipyfft():
    x = np.linspace(-1, 1, 1000)
    Fmid = []
    f = 1000
    k = 2 * np.pi / (0.000633 * f)
    for i in range(len(x)):
        Fmid.append(beams.Ai(x[i]))
    out = sp.fft(x)
    out = np.abs(out)
    plt.plot(x, out)
    plt.show()

def norm(arr):
    norm = 0
    for i in range(len(arr)):
        norm += arr[i]**2
    return norm

def sko(diffArr, sourceArr):
    return norm(diffArr)/norm(sourceArr)


def precision():
    Nx = 500
    Nx1 = 4000
    Nu = 701
    u_left = -2
    u_right = 3
    x_left = -2
    x_right = 2
    f = 400
    z = 400
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    Fmid = [0] * (Nx - 1)
    new_x = []
    max_x = []
    out1 = fourier.fourierArr(Nx, x_right, x_left, f, u, beams.Ai)
    out2 = fourier.fourierArr(Nx1, x_right, x_left, f, u, beams.Ai)
    out1_abs = out1.__abs__()
    out2_abs = out2.__abs__()
    # print(np.linalg.norm(out1-out2))
    # print(np.linalg.norm(out1_abs - out2_abs))
    # print(sko(out1 - out2, out1))
    diff = out1_abs - out2_abs
    plt.plot(u, out1_abs, u, out2_abs)
    plt.show()
    print(norm(out2_abs))
    print(norm(out1_abs))
    print(norm((out1_abs - out2_abs)))
    print(sko((out1_abs - out2_abs), out1_abs))
        # print(beams.getAiParam())
        # a = fresnel_four_comparing(f, z, Nx, x_left, x_right, Nu,  u, beams.Ai)



