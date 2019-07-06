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


class Params:
    pe_param = 1
    Aiparam = 5
    sigma = 2

def setPeParam( param):
    Params.pe_param = param

def setAiParam( param):
    Params.Aiparam = param

def getAiParam():
    return Params.Aiparam

def getPeParam():
    return Params.pe_param

def Pe(x):
    return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2))

def correctPe(x):
    return np.exp(1j * (x ** 4) + 1j * ((getPeParam() * x) ** 2))

def correctAi(x):
    return np.exp(1j * ((getAiParam() * x) ** 3) / 3)

def PeGauss(x):
    return np.exp(x**2/Params.sigma)*np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2))

def Pe_with_param(x, t):
    return np.exp(1j * (x ** 4) + t * 1j * (x ** 2))


def PeOdd(x):
    if (x >= 0):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2))
    else:
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2))


def PeOdd2d(x, y):
    if ((x >= 0) & (y >= 0)):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2) + 1j * (y ** 4) + getPeParam() * 1j * (y ** 2))
    elif ((x >= 0) & (y < 0)):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2) + 1j * -(y ** 4) + getPeParam() * 1j * -(y ** 2))
    elif ((x < 0) & (y > 0)):
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2) + 1j * (y ** 4) + getPeParam() * 1j * (y ** 2))
    else:
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2) + 1j * -(y ** 4) + getPeParam() * 1j * -(y ** 2))



def Gauss(x):
    return np.sqrt(2*np.pi)*Params.sigma*np.exp(-((x ** 2) / (2*Params.sigma ** 2)))

def Gauss1(x):
    return np.exp(-(x ** 2))

def Ai(x):
    return np.exp(1j * getAiParam() * (x ** 3) / 3)


def AiEven(x):
    return np.exp(-1j * getAiParam() * (abs(x) ** 3) / 3)


def AiEven2d(x, y):
    return np.exp(1j * getAiParam() * (abs(x) ** 3) / 3 + 1j * getAiParam() * (abs(y) ** 3) / 3)


def AiI(x):
    return cmath.exp(I * getAiParam() * (x ** 3) + I * x)


def Pe2d(x, y):
    return np.exp(1j * (x ** 4 + y ** 4) + 1j * (x ** 2 + y ** 2))

def Pe2d_right(q, x):
    return np.exp(1j * (q ** 4 ) + 1j * (q ** 2 * x))


def Ai2d(x, y):
    return np.exp(1j * getAiParam() * x ** 3 + 1j * getAiParam() * y ** 3)

def getInitPe2d():
    xright = 1
    xleft = -1
    yright = 1
    yleft = -1
    nx = 151
    ny = 151
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
            xarr.append(Pe2d(xmid[i], ymid[j]) * xstep * ystep)

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



def plotAiPhase():
    a = np.linspace(-1, 1, 100)
    Fabs = list(map(lambda x: np.angle(Ai(x)), a))
    plt.plot(a, Fabs)
    plt.grid()
    plt.show()

def plotPePhase():
    a = np.linspace(-2, 2, 100)
    Fabs = list(map(lambda x: np.angle(PeOdd(x)).real, a))
    plt.plot(a, Fabs)
    plt.grid()
    plt.show()

def Rect(x):
    return 1


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
    for i in range(nx - 1):
        xmid.append((x[i + 1] + x[i]) / 2)
        ymid.append((y[i + 1] + y[i]) / 2)
    xmid = np.array(xmid)
    ymid = np.array(ymid)
    Fphase = []
    for i in range(len(xmid)):
        xarr = []
        for j in range(len(ymid)):
            xarr = np.append(xarr, Ai2d(xmid[i], ymid[j]) * xstep * ystep)
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
    xleft = 0
    yright = 0
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
            xarr.append(Ai2d(x[i], y[j]) * xstep * ystep)
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
    fig = pylab.figure()
    # axes = Axes3D(fig)
    axes = fig.gca(projection='3d', proj_type='ortho')
    axes.set_xlabel('x, mm')
    axes.set_ylabel('y, mm')
    axes.plot_surface(x, y, Fphase, cmap=plt.cm.binary)
    pylab.show()
    # x, y = np.meshgrid(u, u)

    # Fphase = np.array(Fphase)
    # return Fphase, xmid, ymid

