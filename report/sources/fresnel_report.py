def fresnelDot(z, Nx, x_right, x_left, u, func):
    x = np.linspace(x_left, x_right, Nx)
    k = 2 * np.pi / (0.000633)
    temp = (1j * k / (2 * z))
    start = time.clock()
    Fmid = [0] * (Nx - 1)
    for i in range(Nx - 1):
        Fmid[i] = func((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i])
    out = 0
    for i in range(Nx - 1):
        out += Fmid[i] * np.exp(temp * ((x[i] + x[i + 1]) / 2 - u) ** 2)
    out *= np.exp(1j * k * z) * np.sqrt(-1j * k / (2 * np.pi * z))
    return out

def fresnelArr(z, new_u, out_u, Ffour):
    F = [0] * len(out_u)
    wavelength = 0.000633
    k = 2 * np.pi / (wavelength)
    temp = (1j * k / (2 * z))
    step_u = new_u[1] - new_u[0]
    for n in range(len(out_u)):
        F[n] = 0
        out_u_n = out_u[n]
        for i in range( len(new_u)):
            F[n] += Ffour[i] * np.exp(temp * (new_u[i] - out_u_n) ** 2)
    F = list(map(lambda x: x * np.exp(1j*k*z) * np.sqrt(-1j*k/(2*np.pi * z)), F))

    Fabs = list(map(lambda x: abs(x) ** 2, F))
    plt.plot(out_u, Fabs)
    plt.show()
    end = time.clock()
    return Fabs
