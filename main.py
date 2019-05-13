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
    Nx = 4001
    Nu = 8001
    u_left = -5
    u_right = 5
    x_left = -4
    x_right = 4
   # fresnel.fresnel(1500, Nx, Nu, x_right, x_left, u_right, u_left, beams.Ai)
    # fresnel.fresnel(5000, Nx, Nu, x_right, x_left, u_right, u_left, beams.Pe_with_param)
    # utils.fresfour(500, Nx, Nu, x_right, x_left, u_right, u_left, beams.Ai, 500, False)
    utils.accelerationVer2()
    # utils.acceleration()
    # beams.plotAiPhase()