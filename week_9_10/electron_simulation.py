# libraries
import numpy as np
from scipy import ndimage
import math
import matplotlib.pyplot as plt

# define the conditions
height, width = 1000, 1000
box_height = 0.01

# omega definition
w = 2 / (1 + math.pi / height)
print(w)

# definition of empty grid
U = np.zeros((height, width))
R = np.zeros((height, width))

# define deltax and deltay
DELTA_X = box_height / (height - 1)
DELTA_Y = DELTA_X

# Electron charge, mass and time units
e = -1.6e-19 # Coulomb
me = 9.11e-31 # kg
h = 10e-12 # s, time step
tEnd = 6e-9 # s


def placeplate(p1, p2, r_grid, u_grid, potential):
	"""places the plate in the box,
	write a function, because more than one plate could be necessary"""

	# the coordinates of the grid
	x1, y1 = p1
	x2, y2 = p2

	# What is this?
	if x2 < x1:
		t = x2
		x2 = x1
		x1 = t
		t = y2
		y2 = y1
		y1 = t

	a = (y2 - y1) / (x2 - x1)
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
		r_grid[l, j] = 0
		u_grid[l, j] = potential
	plt.plot([x1, x2], [y1, y2])


def thegrid(height, width, w, U, R):
	"""define the grid that should be used"""

	# definition stencil
	stencil = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])

	# initialize boundary conditions -> R -> outside 0, inside w/4
	R.fill(w / 4)
	R[1, :] = 0
	R[-1, :] = 0
	R[:, 1] = 0
	R[:, -1] = 0

	# put center piece it in U / R entries at the same place are 0
	#U[int(height / 2), int(width / 3):-int(width / 3)] = 1000
	#R[int(height / 2), int(width / 3):-int(width / 3)] = 0
	starting_point1 = (0.0025, 0.0025)
	end_point1 = (0.0075, 0.0075)
	placeplate(starting_point1, end_point1, R, U, 1000)

	# chess board
	C = np.zeros((height, width), dtype=bool)
	C[::2, ::2] = True
	C[1::2, 1::2] = True

	# loop
	correction_value = 2

	# steps
	steps = 0
	Out = np.zeros_like(U)

	while abs(correction_value) > 1:
		#Out = np.zeros_like(U)
		ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0)  # stencil angewwendet azu fu
		U[C] = U[C] + R[C] * Out[C]

		correction = R[C] * Out[C]  # this is the correction that is calculated and added every time, it gets so small over time, that is shows that it converged
		correction_value = np.max(correction)
		print(correction_value)

		ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0)  # stencil angewwendet azu fu
		U[~C] = U[~C] + R[~C] * Out[~C]

	return U


def generateelectrons(n):

	y = np.linspace(0.7*box_height, 0.9*box_height, n)
	x = np.zeros_like(y)
	# Take a random angle phi in the range -pi/2 to pi/2,
	# either uniform or normal distributed
	# vx = 1e6 * cos(phi)
	# vy = 1e6 * sin(phi)
	vx = 1e6*np.ones_like(x)
	vy = np.zeros_like(y)
	# can be vectors -> can work with an array of electrons
	return (x, y, vx, vy)

def coordToIndex(x, y):
	# delta_x and delta_y -> grid spacing 1 cm /100 grid cells
	j = np.array(x / DELTA_X, dtype='int')
	l = np.array(y / DELTA_Y, dtype='int')
	return (j, l)

# calculate t and u and calculate the acceleration
def accel (x, y):
	# get the global U, search the error here first!
	global U

	j, l = coordToIndex(x, y)

	# Make sure j, l are inside the grid
	L, J = U.shape
	j = np.maximum(j, np.zeros_like(j))
	j = np.minimum(j, (J-2)*np.ones_like(j))
	l = np.maximum(l, np.zeros_like(l))
	l = np.minimum(l, (L-2)*np.ones_like(l))

	t = (x - j * DELTA_X) / DELTA_X
	u = (y - l * DELTA_Y) / DELTA_Y

	# the phi values
	PHI1 = U[l, j]
	PHI2 = U[l, j + 1]
	PHI3 = U[l + 1, j]
	PHI4 = U[l + 1, j + 1]

	ax = -(e/me) * 1/DELTA_X * ((1-u)*(PHI2-PHI1) + u*(PHI4-PHI3))
	ay = -(e/me) * 1/DELTA_Y * ((1-t)*(PHI3-PHI1) + t*(PHI4-PHI2))

	return (ax, ay)


# insert kick and drift

def drift(x, y, vx, vy):
	xi = x + 1/2 * h * vx
	yi = y + 1/2 * h * vy

	return (xi, yi)


def kick(x, y, vx, vy):
	# calculate the acceleration
	ax, ay = accel(x, y)
	vx = vx + h * ax
	vy = vy + h * ay

	return (vx, vy)


def leapfrog(s):
	x, y, vx, vy = s
	x, y = drift(x, y, vx, vy)
	vx, vy = kick(x, y, vx, vy)
	x, y = drift(x, y, vx, vy)
	return (x, y, vx, vy)


# Try
U_new = thegrid(height, width, w, U, R)
print(U_new)
steps = int(tEnd / h)
num_electrons = 3

x_list = []
y_list = []
s = generateelectrons(num_electrons)

for step in range(steps):
	s = leapfrog(s)
	x, y, vx, vy = s
	x_list.append(x)
	y_list.append(y)

plt.imsave("sweep.png", U_new)
plt.plot(x_list, y_list)
plt.savefig(("electron-paths.png"))

