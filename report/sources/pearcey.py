def Pearcey(x):
    		return cmath.exp(I * (x ** 4) + ( Peparam * I * (x ** 2)))

def PlotInitialAiry(a, b, n):
    		x = np.linspace(a, b, n)
    		peinit = list(map(Pearcey, x))
    		peinitphase = list(map(cmath.phase, airyinit))
		peinitabs = list(map(abs, airyinit))
		peinitreal = list(map(numpy.real, airyinit))
		peinitimag = list(map(numpy.imag, airyinit))
    		plt.plot(x, peinitphase)
    		plt.grid()
    		plt.show()
		plt.plot(x, peinitabs)
    		plt.grid()
    		plt.show()
		plt.plot(x, peinitreal)
    		plt.grid()
    		plt.show()
		plt.plot(x, peinitimag)
    		plt.grid()
    		plt.show()
	}