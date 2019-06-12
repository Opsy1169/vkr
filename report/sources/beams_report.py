class Params:
    pe_param = 10
    Aiparam = 20
    sigma = 1000

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

def PeGauss(x):
    return np.exp(x**2/Params.sigma)*np.exp(1j * (x ** 4)
                                            + getPeParam() * 1j * (x ** 2))

def Pe_with_param(x, t):
    return np.exp(1j * (x ** 4) + t * 1j * (x ** 2))

def PeOdd(x):
    if (x >= 0):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2))
    else:
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2))

def PeOdd2d(x, y):
    if ((x >= 0) & (y >= 0)):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2) + 1j *
                      (y ** 4) + getPeParam() * 1j * (y ** 2))
    elif ((x >= 0) & (y < 0)):
        return np.exp(1j * (x ** 4) + getPeParam() * 1j * (x ** 2) + 1j
                      * -(y ** 4) + getPeParam() * 1j * -(y ** 2))
    elif ((x < 0) & (y > 0)):
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2) + 1j
                      * (y ** 4) + getPeParam() * 1j * (y ** 2))
    else:
        return np.exp(1j * -(x ** 4) + getPeParam() * 1j * -(x ** 2) + 1j
                      * -(y ** 4) + getPeParam() * 1j * -(y ** 2))

def Gauss(x):
    return np.exp(-x ** 2 / Params.sigma ** 2)

def Ai(x):
    return np.exp(1j * getAiParam() * (x ** 3) / 3)


def AiEven(x):
    return np.exp(-1j * getAiParam() * (abs(x) ** 3) / 3)


def AiEven2d(x, y):
    return np.exp(1j * getAiParam() * (abs(x) ** 3) / 3 + 1j
                  * getAiParam() * (abs(y) ** 3) / 3)


def AiI(x):
    return cmath.exp(I * getAiParam() * (x ** 3) + I * x)

def Pe2d(x, y):
    return np.exp(1j * (x ** 4 + y ** 4) + 1j * (x ** 2 + y ** 2))


def Ai2d(x, y):
    return np.exp(1j * getAiParam() * x ** 3 + 1j
                  * getAiParam() * y ** 3)


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


