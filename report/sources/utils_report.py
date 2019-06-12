def accelerationVer2(Nx, x_right, x_left, Nu, u_left, u_right):
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    z_quan = 45
    z_0 = 100
    z_step = 20
    v = np.linspace(z_0, z_0 + z_step * (z_quan - 1), z_quan)
    output = []
    start = time.time()
    Fmid = [0] * (Nx - 1)
    func = beams.PeOdd
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    for i in range(z_quan):
        outputLine = fresnel_four(i * z_step + z_0, Nx, Nu, x, u, Fmid)
        output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end - start))
    u, v = np.meshgrid(u, v)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.view_init(100, 100)
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary)
    pylab.show()

def fresnel_four(z, Nx, Nu, x, u, Fmid):

    F = [0] * Nu

    k = 2 * np.pi / (0.000633)
    f = 150
    start = time.clock()
    new_x = []
    for i in range(Nx - 1):
        new_x.append((x[i + 1] + x[i]) / 2)
    step_x = new_x[1] - new_x[0]
    temp = (-1j * k / (z))
    secondtemp = 1j * (k / 2) * ((f - z) / (f * z))
    for n in range(Nu):
        F[n] = 0
        for i in range(Nx - 1):
            F[n] += Fmid[i] * np.exp(temp * (new_x[i]) * u[n] + secondtemp * (new_x[i]) ** 2        F[n] *= np.exp(-temp/2 * u[n]**2)
        F = list(map(lambda x: x * np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z)), F))
    Fabs = list(map(abs, F))
    end = time.clock()
    print('time:' + str(end - start))

    plt.plot(u, Fabs)
    return Fabs

def normalizeArray(arr):
    max = np.max(arr)
    arr /= max
    return arr


def beamsForDifferentFocus(Nx, x_right, x_left, Nu, u_left, u_right, f, func):
    u = np.linspace(u_left, u_right, Nu)
    x = np.linspace(x_left, x_right, Nx)
    leg = []
    for i in range(len(f)):
        a = fourier.fourierArr(Nx, x_right, x_left, f[i], u, func)
        leg.append('f = ' + str(f[i]) + ', mm')
    plt.legend(leg)
    plt.show()

def beamsForDifferentParam(Nx, x_right, x_left, f, Nu, u_left, u_right, params, func):
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
        leg.append('alpha = ' + str(params[i]))
    plt.legend(leg)
    plt.show()

def beamsForDifferentInputRange(f, Nu, u_left, u_right, rights, lefts, func):
    f = 1000
    u = np.linspace(u_left, u_right, Nu)
    leg = []
    for i in range(len(rights)):
        a = fourier.fourierArr((rights[i]-lefts[i])*800+1, rights[i], lefts[i], f, u, func)
        leg.append('right = ' + str(rights[i]) + "; left = " + str(lefts[i]))
    plt.legend(leg)
    plt.show()
