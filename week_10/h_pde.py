import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.animation import FuncAnimation

# Globals
N = 450
u = 1
dx = 1 / N
h = dx / abs(u)
C = u * h / dx * 0.5

# 1) Using simple Python loop:
# def timeStep1(rho):
#     rhoNew = np.zeros_like(rho)
#     N = len(rho)
#     for j in range(N):
#         jp1 = (j + 1) % N
#         jm1 = (j - 1) % N
#         rhoNew(j) = ... schema ...
#     return rhoNew


# 2) Do not use a for loop!
def lax_step(rho):
    # np.roll(rho, 1):  Roll the array by one position to the right [1,2,3] -> [3,1,2]
    # np.roll(rho, -1):  Roll the array by one position to the left  [1,2,3] ->j [2,3,1]
    # rhoJplus1 = np.roll(rho, ????) Think twice!
    # rhoJminus1 = ...
    # rhoNew = 0.5*rhoJplus1 - 0.5*rhoJminus1

    left = np.roll(rho, 1)
    right = np.roll(rho, -1)
    rho_new = left * 0.5 * (1 + C) + right * 0.5 * (1 - C)
    return rho_new

def upwind_step(rho):
    if u > 0:
        return ndimage.convolve(rho, np.array([0, 1-C, C]), mode="wrap")
    else:
        return ndimage.convolve(rho, np.array([-C, 1 + C, 0]), mode="wrap")


# 3) Using convolve (see SOR)
def lax_wendroff_step(rho):
    #     a_j = c * a_{j-1} + b * a_{j+1}
    #     W = np.array([b 0 c])
    #     aNew = convolve(a, W, mode="periodic")
    return ndimage.convolve(
        rho, np.array([-0.5 * C * (1 - C), 1 - C**2, 0.5 * C * (1 + C)]), mode="wrap"
    )


if __name__ == "__main__":
    # 4) Initialziation:
    fig = plt.figure(constrained_layout=True)
    rho = np.zeros(N)

    # for j in range(int(N/2)-10, int(N/2)+10):
    #     rho[j] = 1
    # rho[int(N / 2) - 10 : int(N / 2) + 10] = 1
    rho[int(N / 4) : int(2 * N / 4)] = 0.4

    x = np.linspace(0, 1, N)

    # for _ in range(5):
    #     rho = lax_step(rho)
    plt.plot(x, rho)

    def update(i):
        global rho
        for i in range(5):
            # rho = lax_step(rho)
            # rho = lax_wendroff_step(rho)
            rho = lax_wendroff_step(rho)
        plt.plot(x, rho)

    ani = FuncAnimation(fig, update, 1000, interval=60, repeat=False)
    # ani.save("lax.mp4")
    plt.show()