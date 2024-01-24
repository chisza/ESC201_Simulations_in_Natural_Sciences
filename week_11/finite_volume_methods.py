# What we want to do?
# We have a grid with a density -> Gaussian
# Implement 2 methods: CIR and CTU

# Finite Difference Methods

# libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, signal
from matplotlib.animation import FuncAnimation

# Initial Condition
Nx = 200
Ny = 200
dx = 1.0 / Nx
dy = 1.0 / Ny
sigma_x = 10 * dx
sigma_y = 10 * dy
x0 = 2 * sigma_x
y0 = 2 * sigma_y
ca = 0.48
cb = 0.48


def ctu(rho):
	"""Corner Transfer"""
	rho_1 = np.roll(rho, 1, axis = 0)
	rho_1[0, :] = 0
	rho_temp = (1 - ca) * rho + ca * rho_1
	print(rho_1)

	rho_2 = np.roll(rho_temp, 1, axis = 1)
	rho_2[:, 0] = 0
	rho_new = (1 - cb) * rho_temp + cb * rho_2
	return rho_new

def cir(rho):
	"""CIR method"""
	rho_1 = np.roll(rho, 1, axis = 0)
	rho_1[0, :] = 0
	rho_2 = np.roll(rho, 1, axis = 1)
	rho_2[:, 0] = 0

	rho_new = rho - ca * (rho - rho_1) - cb * (rho - rho_2)
	return rho_new


def ctu_con(rho):
	"""Corner Transfer with Convolution"""
	rho_new = signal.convolve2d(rho, np.array([[0, 0, 0],
												[0, (1-ca) * (1-cb), ca*(1-cb)],
												[0, cb * (1-ca), ca  * cb]]),
								boundary="fill", mode="same")

	return rho_new

def cir_con(rho):
	"""CIR method with convolution"""
	rho_new = signal.convolve2d(rho, np.array([[0, 0, 0],
												[0, 1-ca-cb, ca],
												[0, cb, 0]]),
								boundary="fill", mode="same")

	return rho_new



def initialCondition():
	rho = np.zeros((Nx, Ny))
	for j in range(Nx):
		for l in range(Ny):
			x = j * dx
			y = l * dy
			rho[l, j] = np.exp(-(x - x0)**2 / (2 * sigma_x**2)
							   - (y - y0)**2 / (2 * sigma_y**2))
	return rho


if __name__ == "__main__":
	fig = plt.figure(constrained_layout=True)
	rho = initialCondition()

	def update(i):
		global rho
		for i in range(20):
			#rho = ctu(rho)
			#rho = cir(rho)
			rho = cir_con(rho)
			#rho = ctu_con(rho)
		plt.imshow(rho)

	ani = FuncAnimation(fig, update, 200, interval=60, repeat=False)
	ani.save("finite_volume_methods_cir.mp4")
	print("Finished")
	plt.show()


