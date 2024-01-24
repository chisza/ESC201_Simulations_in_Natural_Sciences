# exercise on hyperbolic pdes
# linear advection

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.animation import FuncAnimation

# 0. the method that is not working
# 1. LAX
# 2. 1. Order Upwind
# 3. LAX-Wendoff

# solve pdes with convolve

# initial conditions
N = 400 # grid size
u = 1 # constant velocity
dx = 1 / N # grid velocity
h = dx / abs(u) * 0.5 # time step (abs() so that it would work with negative velocities
C = (u * h) / dx # constant

 # Rho is a matrix of size N, initialized to =

# with for loop
def lax_time_step(rho):
	rho = np.zeros(N)
	rho_new = np.zeros_like(rho) # initialize a new rho, the 0s are to be overwritten
	n = len(rho)
	for j in range(n):
		rho_j_plus_1 = (j + 1) % n
		rho_j_minus_1 = (j - 1) % n
		rho_new[j] = 0.5 * (rho_j_plus_1 + rho_j_minus_1) - 0.5 * C * (rho_j_plus_1 - rho_j_minus_1)
	return rho_new

# without for loop
def lax(rho):
	"""Linear advection with LAX"""
	right = np.roll(rho, -1) #N Roll the array by one position to the right [1,2,3] -> [3,1,2]
	left = np.roll(rho, 1) # Roll the array by one position to the left  [1,2,3] ->j [2,3,1]
	rho_new = 0.5 * (right + left) - 0.5 * C * (right - left) # apply the equation
	return rho_new


def first_oder_upwind(rho):
	right = np.roll(rho, -1)  # N Roll the array by one position to the right [1,2,3] -> [3,1,2]
	left = np.roll(rho, 1)  # Roll the array by one position to the left  [1,2,3] ->j [2,3,1]

	if u > 0:
		return rho - C * (rho - left)
	else:
		return rho - C * (right - rho)

def lax_wendoff(rho):
	right = np.roll(rho, -1)  # N Roll the array by one position to the right [1,2,3] -> [3,1,2]
	left = np.roll(rho, 1)  # Roll the array by one position to the left  [1,2,3] ->j [2,3,1]

	if u > 0:
		return 0.5 * C * (1 + C) * left + (1 - C**2) * rho - 0.5 * C * (1 - C) * right

# with SOR
def first_order_upwind_sor(rho):
	if u > 0:
		# rho should be convolved
		# make a stencil from the formula
		return ndimage.convolve(rho, np.array([0, 1-C, C]), mode="wrap")
	else:
		return ndimage.convolve(rho, np.array([-C, 1 + C, 0]), mode="wrap")

def lax_wendroff_sor(rho):
	return ndimage.convolve(
		rho, np.array([-0.5 * C * (1 - C), 1 - C ** 2, 0.5 * C * (1 + C)]), mode="wrap"
	)

def update_parallel(start, end, rho_shared):
	global rho
	for i in range(5):
		rho_shared[start:end] = lax_wendroff_sor(rho[start:end])


if __name__ == "__main__":

	fig = plt.figure(constrained_layout=True)
	rho = np.zeros(N)

	# for j in range(int(N/2)-10, int(N/2)+10):
	#     rho[j] = 1
	rho[int(N / 2) - 10 : int(N / 2) + 10] = 1

	x = np.linspace(0, 1, N)

	# for _ in range(5):
	#     rho = lax_step(rho)
	plt.plot(x, rho)

	def update(i):
		global rho
		print(rho)
		for i in range(20):
			# rho = lax_step(rho)
			# rho = lax_wendroff_step(rho)
			#rho = lax_wendroff_sor(rho)
			rho = first_oder_upwind(rho)
		plt.plot(x, rho)

	ani = FuncAnimation(fig, update, 10, interval=60, repeat=False)
	#ani.save("lax_wendroff_sor.mp4")
	plt.show()

