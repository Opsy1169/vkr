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
    Nx = 1501
    Nu = 1001
    u_left = -2
    u_right = 2
    x_left = -1
    x_right = 1
    # utils. pearceyShit()
    # fourier.Fourier(Nx, Nu, u_left, u_right, x_left, x_right, 200, False, beams.Pe)
    # fresnel.fresnel(5000, Nx, Nu, x_right, x_left, u_right, u_left, beams.Pe)
    # fresnel.fresnel(5000, Nx, Nu, x_right, x_left, u_right, u_left, beams.Pe_with_param)
    # utils.fresfour(400, Nx, Nu, x_right, x_left, u_right, u_left, beams.Ai, 500, False)
    # beams.plotPhasetAi2d()
    # utils.beamsForDifferentParam([5, 20], beams.Ai)
    # utils.precision()
    utils.beamsForDifferentInputRange([5, 2], [-5, -2], beams.Ai)
    # utils.beamsForDifferentFocus([100, 150, 200, 250, 300, 350], beams.Ai)
    # 250, 200, 150, 100
    # 450, 400, 350, 300
    # utils.accelerationVer2()
    # utils.accelerationCompare()
    # utils.plotshit()
    # utils.sinc()
    # utils.Pearcey()
    # fourier.Fourier2d()
    # utils.acceleration()
    # beams.plotAiPhase()
    # fourier.Fourier2d()