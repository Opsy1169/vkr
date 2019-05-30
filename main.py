import beams
import fresnel
import fourier
import utils

if __name__ == '__main__':
    # utils.plotSmthn()
    # utils.Pearcey()
    # utils.PearceyOdd()
    # utils.Airy()

    # Fourier2d()
    # halfPearcey()
    # fresnel_four()
    # plotPhasetAi2d()
    # fastfour()
    # acceleration()
    # getInitPeOdd2d()
    # plotPePhase()
    # getInitAiEven2d()
    Nx = 5001
    Nu = 2001
    u_left = -1
    u_right = 4
    x_left = -1
    x_right = 1
    # fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, 1000, False, beams.Pe)
    # fresnel.fresnel(5000, Nx, Nu, x_right, x_left, u_right, u_left, beams.Pe)
    # fresnel.fresnel(5000, Nx, Nu, x_right, x_left, u_right, u_left, beams.Pe_with_param)
    # utils.fresfour(400, Nx, Nu, x_right, x_left, u_right, u_left, beams.Ai, 500, False)
    # beams.plotPhasetAi2d()
    # utils.beamsForDifferentInputRange([2, 5], [-2, -5], beams.Pe)
    utils.accelerationVer2()
    # utils.acceleration()
    # beams.plotAiPhase()
    # fourier.Fourier2d()