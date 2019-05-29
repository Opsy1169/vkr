for i in range(Nx - 1):
    F[n] += Fmid[i] * np.exp((1j * k /(2*z))*
                             ((x[i] + x[i + 1])/2 - u[n])**2
