

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


pe_param = 20

output = 0

def Pe(x):
    return np.exp(1j * (x ** 4) + pe_param* 1j * (x ** 2))


def PeOdd(x):
    if(x >= 0):
        return np.exp(1j * (x ** 4) + pe_param * 1j * (x ** 2))
    else:
        return np.exp(1j * -(x ** 4) + pe_param * 1j * -(x ** 2))

def PeOdd2d(x, y):
    if((x >= 0) & (y >= 0)):
        return np.exp(1j * (x ** 4) + pe_param * 1j * (x ** 2) + 1j * (y ** 4) + pe_param * 1j * (y ** 2))
    elif((x >= 0) & (y < 0)):
        return np.exp(1j * (x ** 4) + pe_param * 1j * (x ** 2) + 1j * -(y ** 4) + pe_param * 1j * -(y ** 2))
    elif ((x < 0) & (y > 0)):
        return np.exp(1j * -(x ** 4) + pe_param * 1j * -(x ** 2) + 1j * (y ** 4) + pe_param * 1j * (y ** 2))
    else:
        return np.exp(1j * -(x ** 4) + pe_param * 1j * -(x ** 2) + 1j * -(y ** 4) + pe_param * 1j * -(y ** 2))

Aiparam = 30

sigma = 0.1
def Gauss(x):
    return np.exp(-x**2/sigma**2)

def Ai(x):
    return np.exp(-1j * Aiparam * (x ** 3)/3)

def AiEven(x):
    return np.exp(-1j * Aiparam * (abs(x) ** 3)/3)

def AiEven2d(x, y):
    return np.exp(1j * Aiparam * (abs(x) ** 3)/3 + 1j * Aiparam * (abs(y) ** 3)/3)

def AiI(x):
    return cmath.exp(I * Aiparam * (x ** 3) + I*x)

def Pe2d(x, y):
    return np.exp(1j*(x**4 + y**4) + 1j*(x**2 + y**2))
def Ai2d(x, y):
    return np.exp(1j*Aiparam*x**3 + 1j*Aiparam*y**3)

def getInitPe2d():
    xright = 2
    xleft = -2
    yright = 2
    yleft = -2
    nx = 151
    ny = 151
    x = np.linspace(xleft, xright, nx)
    y = np.linspace(yleft, yright, ny)
    xstep = abs(x[1] - x[0])
    ystep = abs(y[1] - y[0])
    xmid = []
    ymid = []
    for i in range(nx-1):
        xmid.append((x[i+1] + x[i])/2)
        ymid.append((y[i+1] + y[i])/2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
             xarr.append(Pe2d(xmid[i], ymid[j])*xstep*ystep)

        Fphase.append(xarr)
    # x, y = np.meshgrid(u, u)

    return Fphase, xmid, ymid
    #
    # x, y = np.meshgrid(x[0:nx-1:1], y[0:ny-1:1])
    # a= 1
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, Fphase, cmap=plt.cm.jet)
    # pylab.show()


def plotPePhase():
    a = np.linspace(-2, 2, 100)
    Fabs = list(map(lambda x: np.angle(PeOdd(x)).real, a))
    plt.plot(a, Fabs)
    plt.grid()
    plt.show()


def getInitAi2d():
    xright = 2
    xleft = -2
    yright = 2
    yleft = -2
    nx = 151
    ny = 151
    x = np.linspace(xleft, xright, nx)
    y = np.linspace(yleft, yright, ny)
    xstep = abs(x[1] - x[0])
    ystep = abs(y[1] - y[0])
    xmid = []
    ymid = []
    for i in range(nx-1):
        xmid.append((x[i+1] + x[i])/2)
        ymid.append((y[i+1] + y[i])/2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
            xarr = np.append(xarr, Ai2d(xmid[i], ymid[j])*xstep*ystep)
        Fphase.append(xarr)

    # Fphase = np.array(Fphase)
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, Fphase, cmap=plt.cm.jet)
    # pylab.show()
    # x, y = np.meshgrid(u, u)

    Fphase = np.array(Fphase)
    return Fphase, xmid, ymid

def getInitPeOdd2d():
    xright = 2
    xleft = -2
    yright = 2
    yleft = -2
    nx = 501
    ny = 501
    x = np.linspace(xleft, xright, nx)
    y = np.linspace(yleft, yright, ny)
    xstep = abs(x[1] - x[0])
    ystep = abs(y[1] - y[0])
    xmid = []
    ymid = []
    for i in range(nx - 1):
        xmid.append((x[i + 1] + x[i]) / 2)
        ymid.append((y[i + 1] + y[i]) / 2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
            xarr.append(PeOdd2d(x[i], y[j]) * xstep * ystep)
        Fphase.append(xarr)

    for i in range(nx - 1):
        for j in range(ny - 1):
            a = Fphase[i][j]
            a = np.angle(a).real
            Fphase[i][j] = a
    Fphase = np.array(Fphase)

    x, y = np.meshgrid(x[0:nx - 1:1], y[0:ny - 1:1])
    a = 1
    Fphase = np.array(Fphase)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    # ax.set_xlabel('X')
    ax.plot_surface(x, y, Fphase, cmap=plt.cm.binary, antialiased=False)
    pylab.show()

def getInitAiEven2d():
    xright = 2
    xleft = -2
    yright = 2
    yleft = -2
    nx = 501
    ny = 501
    x = np.linspace(xleft, xright, nx)
    y = np.linspace(yleft, yright, ny)
    xstep = abs(x[1] - x[0])
    ystep = abs(y[1] - y[0])
    xmid = []
    ymid = []
    for i in range(nx - 1):
        xmid.append((x[i + 1] + x[i]) / 2)
        ymid.append((y[i + 1] + y[i]) / 2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
            xarr.append(AiEven2d(x[i], y[j]) * xstep * ystep)
        Fphase.append(xarr)

    for i in range(nx - 1):
        for j in range(ny - 1):
            a = Fphase[i][j]
            a = np.angle(a).real
            Fphase[i][j] = a
    Fphase = np.array(Fphase)

    x, y = np.meshgrid(x[0:nx - 1:1], y[0:ny - 1:1])
    a = 1
    Fphase = np.array(Fphase)
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(x, y, Fphase, cmap=plt.cm.binary)
    # pylab.show()
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.plot_surface(x, y, Fphase, cmap=plt.cm.binary, antialiased=False)
    pylab.show()


def plotPhasetAi2d():
    xright = 2
    xleft = -2
    yright = 2
    yleft = -2
    nx = 501
    ny = 501
    x = np.linspace(xleft, xright, nx)
    y = np.linspace(yleft, yright, ny)
    xstep = abs(x[1] - x[0])
    ystep = abs(y[1] - y[0])
    xmid = []
    ymid = []
    for i in range(nx-1):
        xmid.append((x[i+1] + x[i])/2)
        ymid.append((y[i+1] + y[i])/2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
            xarr.append( Pe2d(x[i], y[j])*xstep*ystep)
        Fphase.append(xarr)

    for i in range(nx - 1):
        for j in range(ny - 1):
            a = Fphase[i][j]
            a = np.angle(a).real
            Fphase[i][j] = a
    Fphase = np.array(Fphase)

    x, y = np.meshgrid(x[0:nx-1:1], y[0:ny-1:1])
    a = 1
    Fphase = np.array(Fphase)
    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(x, y, Fphase, cmap=plt.cm.binary)
    pylab.show()
    # x, y = np.meshgrid(u, u)

    # Fphase = np.array(Fphase)
    # return Fphase, xmid, ymid

def Fourier2d():

    init, x, y = getInitPe2d()
    start = time.clock()
    # init, x, y = getInitAi2d()
    uleft = -2
    uright = 2
    vleft = -2
    vright = 2

    nv = 101
    nu = 101
    koef = 2 * np.pi / (0.000633 * 1000)
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
            for k in range(len_x):#x
                x_k = x[k]
                x_u = x_k*u_i
                for p in range(len_y):#y
                    integrateelem += init[k][p] * np.exp( koef_mul_i *(x_u + y[p]*v_j))
            outputline.append( integrateelem)
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
    Nx = 10001
    Nu = 2001
    u_left = -1
    u_right = 10
    x_left = -7
    x_right = 7
    # zeros = np.zeros(2000)
    x = np.linspace(x_left, x_right, Nx)
    # x = np.concatenate((zeros, x))
    F = [0] * Nu
    f = 1000
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633 * f)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = PeOdd((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Fmid[i] * np.exp((-1j * k * (x[i] + x[i + 1]) * u[j]) / 2)
        # F[j] = F[j] * (x_right-x_left) / (Nx)
    print("done")
    F = list(map(lambda x: x*(x_right-x_left) / (Nx), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print("time:" + str(end-start))
    # another = list(map(np.real, F))
    # anotheranother =list(map(np.imag, F))

    F2dAbs = []

    # for i in range(Nu):
    #     row = [j * Fabs[i] for j in Fabs]
    #     F2dAbs.append(row)
    #
    # arrayF2dAbs = np.array(F2dAbs)
    # z, y = np.meshgrid(u, u)
    #
    # fig = pylab.figure()
    # axes = Axes3D(fig)
    # axes.plot_surface(z, y, arrayF2dAbs, cmap=plt.cm.jet)
    # pylab.savefig("plots/pe3d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param( " + str(pe_param) + ").png")
    # pylab.show()
    # print("another done")
    plt.plot(u, Fabs)
    plt.grid()
    plt.legend(['Pearcey output abs: f = ' + str(f) + "; param = " + str(pe_param)])
    plt.savefig("plots/pe_odd_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ").png")
    plt.show()



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
    k = 2*np.pi/( l )
    #
    x0 = l
    zc = (2 * k * l ** 2)
    # z = 75*l/8
    z = (2 * k * l ** 2)*0.9

    # for j in range(Nu):
    #     F[j] = 0
    #     for i in range(Nx - 1):
    #         F[j] +=  np.exp((1j * ((x[i] + x[i + 1])/2)**2 * u[j]/x0) + 1j*(1 - z/zc)*((x[i] + x[i + 1])/2)**4)*((x[i+1]) - x[i])

    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] +=  np.exp((1j * ((x[i] + x[i + 1])/2)**2 * u[j]/(x0*(1- z/zc)**(1/2))) + 1j*((x[i] + x[i + 1])/2)**4)*((x[i+1]) - x[i])

    # F[j] = F[j] * (x_right-x_left) / (Nx)
    print("done")
    F = list(map(lambda x: x /((1 - z/zc)**(1/4)), F))
    Fabs = list(map(abs, F))
    plt.plot(u, Fabs)
    plt.grid()
    # plt.savefig("plots/pe_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(
    #     x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ").png")
    plt.show()

def Airy():
    a = 3
    b = 1
    Nx = 2501
    Nu = 2001
    u_left = -.1
    u_right = 3
    x_left = -2
    x_right = 2
    f = 200
    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2*np.pi/(0.000633*f)

    # print(x)

    # print(u)
    start = time.clock()
    for j in range(Nu):
        F[j] = 0
        for i in range(Nx - 1):
            F[j] += Ai((x[i] + x[i + 1]) / 2) * np.exp((-1j * k * ((x[i] + x[i + 1])) * u[j])/2) * (x[i + 1] - x[i])#сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье


    # F = list(map(lambda x: np.sqrt(k)*x, F))
    G = [0]*Nu
    # mp.dps = 50
    # for i in range(Nu):
    #     G[i] = Float(abs(mp.quad(lambda x: Ai(x) * cmath.exp(-I * k * x * u[i]), [-1, 3])), 35)
    #     print(G[i])
    # print(F)
    # print(F)
    # print("airy")
    # print(k, F[0], u[0])
    end = time.clock()

    print("TIME: ", str(end-start))
    Fabs = list(map(abs, F))
    # print(Fabs)

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
    # pylab.savefig("plots/ai2d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(Aiparam) + ").png")
    # pylab.show()
    # print("another done")
    plt.plot(u, Fabs)
    plt.grid()
    plt.legend(['Airy output abs: f = ' + str(f) + "; param = " + str(Aiparam)])
    plt.savefig("plots/ai_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + "param(" + str(Aiparam) + ")_f(" + str(f) + ").png")
    plt.show()

def InitialAiryPhaze():
    x = np.linspace(-2, 2, 400)
    airyinit = list(map(Ai, x))
    airyinitphase = list(map(np.real, airyinit))
    plt.plot(x, airyinitphase)
    plt.grid()
    plt.show()


def fresnel(z, Nx, Nu, x_right, x_left, u_right, u_left):
    a = 3
    b = 1
    # Nx = 2001
    # Nu = 1001
    # u_left = -5
    # u_right = 5
    # x_left = -4
    # x_right = 4

    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633)
    # z = 800
    start = time.clock()
    Fmid = [0]*(Nx-1)
    for i in range(Nx - 1):
        Fmid[i] = Ai((x[i] + x[i + 1]) / 2)* (x[i + 1] - x[i])
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * np.exp((1j * k /(2*z))*
                                     ((x[i] + x[i + 1])/2 - u[n])**2)   # сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
    F = list(map(lambda x: x*np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z)), F))
    print(np.sqrt(-1j * k / (2 * np.pi * z)))
    print('a')

    Fabs = list(map(lambda x: abs(x)**2, F))
    end = time.clock()

    print('time:' + str(end-start))


    # plt.plot(u, Fabs)
    # plt.grid()
    # plt.savefig("plots/fresnel_pe_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ")_z(" + str(z) + ").png")
    # plt.show()
    return Fabs


def acceleration():
    Nx = 2001
    Nu = 401
    u_left = -4
    u_right = 4
    x_left = -4
    x_right = 4
    u = np.linspace(u_left, u_right, Nu)

    z_quan = 30
    repeat_time = 40
    z_0 = 400
    z_step = 40
    v = np.linspace(z_0, z_0+z_step*z_quan-1, z_quan*repeat_time)
    output = []
    start = time.time()
    for i in range(z_quan):
        outputLine = fresnel(i*z_step+z_0, Nx, Nu, x_right, x_left, u_right, u_left)
        for j in range(repeat_time):
            output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end-start))
    u, v = np.meshgrid(u, v)
    # fig = plt.figure()
    # axes = Axes3D(fig)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type = 'ortho')
    ax.view_init(100, 100)
    # axes = fig.add_subplot(122, projection='3d')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary,antialiased=False)
    pylab.show()



def fresnel_four():
    a = 3
    b = 1
    Nx = 3001
    Nu = 1501
    u_left = -.1
    u_right = 8
    x_left = -4
    x_right = 4

    x = np.linspace(x_left, x_right, Nx)
    F = [0] * Nu
    u = np.linspace(u_left, u_right, Nu)
    k = 2 * np.pi / (0.000633)
    z = 1200
    f = 650
    # print(x)

    # print(u)
    a = (1.777514461805397 - 1.777514461805397j)
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = Ai((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * np.exp((-1j * k / (z)) * ((x[i] + x[i + 1]) / 2) * u[n] - 1j * (k / 2) * ((f-z)/(f*z)) * ((x[i] + x[i + 1]) / 2)**2) # сейчас в формуле в экспоненте стоит -I, чтобы было похоже на результат Фурье
        F[n] *= np.exp(-1j * (k / 2*z) * u[n]**2)
        # F[n] *=
        # print(F[j])
    F = list(map(lambda x: x * np.exp(1j * k * z) * np.sqrt(1j * k / (2 * np.pi * z)), F))
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
    # F2dAbs = []
    Fabs = list(map(abs, F))
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
    # pylab.savefig("plots/ai2d_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(x_right) + ", " + str(Nx) + ")_param(" + str(Aiparam) + ").png")
    # pylab.show()


    # end = time.clock()
    # print('time:' + str(end - start))
    plt.plot(u, Fabs)
    plt.grid()
    plt.savefig(
        "plots/fresnel_four_ai_u(" + str(u_left) + ", " + str(u_right) + ", " + str(Nu) + ")_x(" + str(x_left) + ", " + str(
            x_right) + ", " + str(Nx) + ")_param(" + str(pe_param) + ")_z(" + str(z) + ")_f(" + str(f) + ").png")
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


def fastfour():
    x = np.linspace(-5, 5, 200)
    a = [1, 2, 3, 4]
    print(np.fft.fftshift(a))
    f = list(map(lambda a: np.exp(-a**2), x))
    F = np.fft.fft(f)
    F = np.fft.fftshift(F)
    plt.plot(x, F)
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
# fresnel()
# scipyfresnel()
# getInitPe2d()
# getInitAi2d()
# Fourier2d()
# halfPearcey()
# fresnel_four()
# plotPhasetAi2d()
# fastfour()
# acceleration()
getInitPeOdd2d()
# plotPePhase()
# getInitAiEven2d()
# Nx = 2001
# Nu = 401
# u_left = -4
# u_right = 4
# x_left = -4
# x_right = 4
# fresnel(400, Nx, Nu, x_right, x_left, u_right, u_left)