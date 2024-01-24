import numpy as np
from scipy import ndimage  # , signal
import matplotlib.pyplot as plt

N = 100
omega = 2 / (1 + (np.pi / N))
BOX_HEIGHT = 0.01
DELTA_X = BOX_HEIGHT / (N - 1)
print(DELTA_X)
DELTA_Y = DELTA_X
h = 1e-12  # s, time step
tEnd = 4000 * h  # s
VE = 1e6
e = -1.6e-19  # Coulomb
me = 9.11e-31  # kg


def placePlate(p1, p2, R, used_U, potential):
	"""
	Parameters
	----------
	p1 : tuple (x,y)
		Starting poing of line segment
	p2 : tuple (x,y)
		End point of line segment
	R : (N,N) matrix
		Defines the boundary conditions, zero entries at the boundary
	U : (N,N) matrix
		Potential matrix.
	potential : float
		Potential value of the plate.

	Returns
	-------
	None.

	"""
	x1, y1 = p1
	x2, y2 = p2
	if x2 < x1:
		t = x2
		x2 = x1
		x1 = t
		t = y2
		y2 = y1
		y1 = t
	print(p1, p2)
	a = (y2 - y1) / (x2 - x1)
	# a =0 # NOTE: use this for points, instead of bars
	j1 = int(x1 / DELTA_X)
	j2 = int(x2 / DELTA_X)
	l1 = int(y1 / DELTA_Y)
	l2 = int(y2 / DELTA_Y)

	n = max(j2 - j1 + 1, l2 - l1 + 1)
	for i in range(n + 1):
		x = x1 + i * (x2 - x1) / n
		y = y1 + a * (x - x1)
		j = int(x / DELTA_X)
		l = int(y / DELTA_Y)
		R[l, j] = 0
		used_U[l, j] = potential
	plt.plot([x1, x2], [y1, y2])


# grid to be modified
def SOR():
	U = np.zeros((N, N))
	weight = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])
	# setting correction grid
	R = np.copy(U)
	R[:, :] = omega / 4
	R[0, :] = 0
	R[N - 1, :] = 0
	R[:, 0] = 0
	R[:, N - 1] = 0

	starting_point1 = (0.0025, 0.0025)
	end_point2 = (0.0075, 0.0075)
	placePlate(starting_point1, end_point2, R, U, 1000)
	# starting_point1 = (0.0035, 0.0025)
	# end_point2 = (0.0085, 0.0060)
	# placePlate(starting_point1, end_point2, R, U, -1000)

	B = np.zeros((N, N), dtype=bool)
	B[::2, ::2] = True
	B[1::2, 1::2] = True

	Out = np.zeros_like(U)
	correction = np.ones_like(U)
	step = 0

	while np.max(np.abs(correction)) > 0.1:
		correction = R * ndimage.convolve(U, weight, output=Out, mode="constant")
		U[B] += correction[B]
		step += 1
		correction = R * ndimage.convolve(U, weight, output=Out, mode="constant")
		U[~B] += correction[~B]
		step += 1

	return U


def generate_electrons(n):
	y = np.linspace(0.6 * BOX_HEIGHT, 0.9 * BOX_HEIGHT, n)
	x = np.zeros_like(y)
	# Take a random angle phi in the range -pi/2 to pi/2,
	# either uniform or normal distributed
	# vx = VE * np.cos(phi)
	# vy = VE * np.sin(phi)
	vx = VE * np.ones_like(x)
	vy = np.zeros_like(y)
	return (x, y, vx, vy)


def coordToIndex(x, y):
	j = np.array(x / DELTA_X, dtype="int")
	l = np.array(y / DELTA_Y, dtype="int")
	return (j, l)


# How to calculate t and u and calculate the acceleration
# -------------------------------------------------------
def accel(x, y):
	global U
	j, l = coordToIndex(x, y)

	# Make sure j,l are inside grid
	J, L = U.shape
	j = np.maximum(j, np.zeros_like(j))
	j = np.minimum(j, (J - 2) * np.ones_like(j))
	l = np.maximum(l, np.zeros_like(l))
	l = np.minimum(l, (L - 2) * np.ones_like(l))

	t = (x - j * DELTA_X) / DELTA_X
	u = (y - l * DELTA_Y) / DELTA_Y
	U1 = U[l, j]
	U2 = U[l, j + 1]
	U3 = U[l + 1, j]
	U4 = U[l + 1, j + 1]
	# U1 = U[j, l]
	# U2 = U[j + 1, l]
	# U3 = U[j, l + 1]
	# U4 = U[j + 1, l + 1]

	ax = -(e / me) * (1 / DELTA_X) * ((1 - u) * (U2 - U1) + u * (U4 - U3))
	ay = -(e / me) * (1 / DELTA_Y) * ((1 - t) * (U3 - U1) + t * (U4 - U2))
	return (ax, ay)


def drift(x, y, vx, vy):
	xi = x + h * 0.5 * vx
	yi = y + h * 0.5 * vy
	return (xi, yi)


def kick(x, y, vx, vy):
	# print("1vx, vy", vx, vy)
	ax, ay = accel(x, y)
	vx = vx + h * ax
	vy = vy + h * ay
	# print("vx, vy", vx, vy)
	return (vx, vy)


def leap_frog(s):
	x, y, vx, vy = s
	x, y = drift(x, y, vx, vy)
	vx, vy = kick(x, y, vx, vy)
	x, y = drift(x, y, vx, vy)
	return (x, y, vx, vy)


if __name__ == "__main__":
	U = SOR()
	steps = int(tEnd / h)
	num_electrons = 3

	x_list = []
	y_list = []
	s = generate_electrons(num_electrons)

	for step in range(steps):
		s = leap_frog(s)
		x, y, vx, vy = s
		x_list.append(x)
		y_list.append(y)

	plt.imsave("sweep.png", U)
	plt.plot(x_list, y_list)
	plt.savefig(("electron-paths.png"))
	#plt.show()