def acceleration(Nx, x_right, x_left, Nu, u_left, u_right):
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



def accelerationVer3(Nx, x_right, x_left, Nu, u_left, u_right):
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
        outputLine = fresnel.fresnelArr(i * z_step + z_0,new_u, out_u,  Ffour)
        output.append(outputLine)

    output = np.array(output)
    end = time.time()
    print('acc time: ' + str(end - start))
    out_u, v = np.meshgrid(out_u, v)
    fig = plt.figure()
    ax = fig.gca(projection='3d', proj_type='ortho')
    ax.view_init(100, 100)
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')
    ax.plot_surface(u, v, output, cmap=plt.cm.binary)
    pylab.show()

def fresfour(z, Nx, Nu, x_right, x_left, u_right, u_left, func, f, is2d):
    F = []
    u = np.linspace(u_left, u_right, Nu)
    out_u = np.linspace(1.1 * u_left, 1.1 * u_right, int(0.5 * Nu))
    k = 2 * np.pi / (0.000633)
    temp = (1j * k / (2 * z))
    start = time.clock()
    new_u = []
    for i in range(Nu - 1):
        new_u.append((u[i + 1] + u[i]) / 2)
    Ffour = fourier.fourierArr(Nx, x_right, x_left, f, new_u, func)
    print("done args")
    end = time.clock()
    print('time:' + str(end - start))
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
        ax = fig.gca(projection='3d', proj_type='ortho')        ax.plot_surface(z, y, arrayF2dAbs, cmap=plt.cm.jet)        		pylab.show()
        print("another done")
    else:
        plt.plot(out_u, Fabs)
        plt.xlabel('x, mm')
        plt.ylabel('y, mm')
        plt.grid()
        plt.show()



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
