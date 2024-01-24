# libraries
import numpy as np
from scipy import ndimage
import math
import matplotlib.pyplot as plt

# the overall conditions
# height and width of the
height, width = 100, 100

# omega definition
w = 2 / (1 + math.pi / height)

# the height of a single cell
BOX_HEIGHT = 0.01 # in meters is the whole grid!

# calculation of the delta for both x and y axis
# is a subdivision for the grid cells
# size of the subdivisions of the grid
DELTA_X = BOX_HEIGHT / (height - 1)
DELTA_Y = DELTA_X

# Electron charge, mass and time units
e = -1.6e-19 # Coulomb
me = 9.11e-31 # kg
h = 10e-12 # s, time step
tEnd = 6e-9 # s

# a function to place the plate in the grid

def bresenham_line(x0, y0, x1, y1):
	points = []
	dx = abs(x1 - x0)
	dy = abs(y1 - y0)
	sx = 1 if x0 < x1 else -1
	sy = 1 if y0 < y1 else -1
	err = dx - dy

	while True:
		# Bresenham algorithm
		points.append((x0, y0))
		print(points)

		if x0 == x1 and y0 == y1:
			break

		e2 = 2 * err
		if e2 > -dy:
			err -= dy
			x0 += sx
		if e2 < dx:
			err += dx
			y0 += sy

	return points

def place_plate(p1, p2, R, used_U, potential):
	x1, y1 = p1
	x2, y2 = p2

	# Ensure x1, y1 is the leftmost point
	if x2 < x1:
		x1, y1, x2, y2 = x2, y2, x1, y1

	# Calculate grid points using Bresenham's algorithm
	points = list(zip(bresenham_line(int(x1 / DELTA_X), int(y1 / DELTA_Y), int(x2 / DELTA_X), int(y2 / DELTA_Y))))

	# Set potential for each grid point
	for point in points:
		R[point[0][1], point[0][0]] = 0
		used_U[point[0][1], point[0][0]] = potential

# a function to generate the grid with the potential

def grid_to_modify(height, width, w):
	"""Function for the grid to be modified and where the
	plates should be placed in: calculate the surrounding potential
	Does the over relaxation"""

	# definition of empty grid
	U = np.zeros((height, width))

	# definition stencil
	stencil = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])

	# initialize boundary conditions -> R -> outside 0, inside w/4
	R = np.zeros((height, width))
	R.fill(w / 4)
	R[0, :] = 0
	R[-1, :] = 0
	R[:, 0] = 0
	R[:, -1] = 0

	# placing the plate
	# coordinates within 1 cm (is box size)

	start_1 = (0.09 * BOX_HEIGHT, 0.72 * BOX_HEIGHT)
	end_1 = (0.25 * BOX_HEIGHT, 0.58 * BOX_HEIGHT)
	place_plate(start_1, end_1, R, U, 1000)

	start_2 = (0.93 * BOX_HEIGHT, 0.31 * BOX_HEIGHT)
	end_2 = (0.82 * BOX_HEIGHT, 0.31 * BOX_HEIGHT)
	place_plate(start_2, end_2, R, U, 1000)

	start_3 = (0.93 * BOX_HEIGHT, 0.08 * BOX_HEIGHT)
	end_3 = (0.61 * BOX_HEIGHT, 0.08 * BOX_HEIGHT)
	place_plate(start_3, end_3, R, U, 1000)


	# chess board
	C = np.zeros((height, width), dtype=bool)
	C[::2, ::2] = True
	C[1::2, 1::2] = True

	# loop
	correction_value = 2
	print("The calculation")
	while abs(correction_value) > 1:
		Out = np.zeros_like(U)
		ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0)  # stencil angewwendet azu fu
		U[C] = U[C] + R[C] * Out[C]

		correction = R[C] * Out[C]  # this is the correction that is calculated and added every time, it gets so small over time, that is shows that it converged
		correction_value = np.max(correction)
		print(correction_value)

		ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0)  # stencil angewwendet azu fu
		U[~C] = U[~C] + R[~C] * Out[~C]

	# returns a full grid with placed plates and the calculated potential
	# everywhere
	return U

# generate the electrons
def generate_electrons(n):
	y = np.linspace(0.7*BOX_HEIGHT, 0.9*BOX_HEIGHT, n)
	x = np.zeros_like(y)
	#x, y = 0.001, 0.005
	# Take a random angle phi in the range -pi/2 to pi/2,
	# either uniform or normal distributed
	phi = np.random.uniform(-np.pi/2, np.pi/2)
	vx = 1e6 * math.cos(phi)
	vy = 1e6 * math.sin(phi)
	#vx = 1e6*np.ones_like(x)
	#vy = np.zeros_like(y)
	# can be vectors -> can work with an array of electrons
	return (x, y, vx, vy)


# coordinates that should be indexed
def coordToIndex(x, y): # x and y are in meters -> search in which häuschen -> j, l = index
	# delta_x and delta_y -> grid spacing 1 cm /100 grid cells
	# calculates from the coordinates the next eckpunkt
	j = np.array(x / DELTA_X, dtype='int')
	l = np.array(y / DELTA_Y, dtype='int')
	return (j, l)


# calculate t and u and calculate the acceleration
def accel (x, y):
	# get the global U, search the error here first!
	global my_grid
	U = my_grid
	j, l = coordToIndex(x, y) # find the potential from the häuschen

	# Make sure j, l are inside the grid
	L, J = U.shape
	j = np.maximum(j, np.zeros_like(j))
	j = np.minimum(j, (J-2)*np.ones_like(j))
	l = np.maximum(l, np.zeros_like(l))
	l = np.minimum(l, (L-2)*np.ones_like(l))

	# calculation relative to eckpunkt from häuschen I'm currently in
	t = (x - j * DELTA_X) / DELTA_X
	u = (y - l * DELTA_Y) / DELTA_Y

	# the phi values
	PHI1 = U[l, j]
	PHI2 = U[l, j + 1]
	PHI3 = U[l + 1, j]
	PHI4 = U[l + 1, j + 1]
	ax = -(e/me) * 1.0/DELTA_X * ((1-u)*(PHI2-PHI1) + u*(PHI4-PHI3))
	ay = -(e/me) * 1.0/DELTA_Y * ((1-t)*(PHI3-PHI1) + t*(PHI4-PHI2))

	return (ax, ay)


# insert kick and drift
def drift(x, y, vx, vy):
	xi = x + 0.5 * h * vx
	yi = y + 0.5 * h * vy

	return (xi, yi)


def kick(x, y, vx, vy):
	# calculate the acceleration
	ax, ay = accel(x, y)
	vx_i = vx + h * ax
	vy_i = vy + h * ay

	return (vx_i, vy_i)


def leapfrog(s):
	x, y, vx, vy = s
	x, y = drift(x, y, vx, vy)
	vx, vy = kick(x, y, vx, vy)
	x, y = drift(x, y, vx, vy)
	return (x, y, vx, vy)


# a try to make this work
my_grid = grid_to_modify(height, width, w) # SOR

# places to store the x and y values of the electrons
x_list = []
y_list = []

# generate a bunch of electrons
num_electrons = 20
s = generate_electrons(num_electrons)
ind = []


# calculate the positions for all the electrons
steps = int(tEnd / h)
times = []
y_hit = []
for step in range(steps):
	print(f"step : {step}")
	s = leapfrog(s)
	x, y, vx, vy = s
	# check if electron still in box
	x_valid = []
	y_valid = []
	vx_valid = []
	vy_valid = []

	#ind = []
	for i in range(len(x)):
		if i not in ind:
			# if it not outside the box
			if 0 < x[i] < BOX_HEIGHT and 0 < y[i] < BOX_HEIGHT:
				x_valid.append(x[i] / DELTA_X)
				y_valid.append(y[i] / DELTA_Y)
				print(i)
			# if it is outside the box
			else:
				y_hit.append(y[i] / DELTA_Y)
				#print(y[i])
				if 10 <= y[i] / DELTA_Y <= 30:
					ti = step * h
					times.append(ti)
				x_valid.append(None)
				y_valid.append(None)
				ind.append(i)
		else:
			x_valid.append(None)
			y_valid.append(None)

	x_list.append(np.array(x_valid))
	y_list.append(np.array(y_valid))



#plt.imshow(my_grid, cmap = "plasma", origin = "lower")
#plt.show()

plt.imshow(my_grid, cmap = "plasma", origin = "lower")
plt.plot(x_list, y_list, "r")
plt.show()

# histogram of the times
plt.hist(times)
plt.xlabel("Time till detector (s)")
plt.ylabel("Number of electrons")
plt.show()

# plot y-axis of the hit on the wall
plt.hist(y_hit)
plt.xlabel("Place the electron hits the box [micro-m]")
plt.ylabel("Number of electrons")
plt.show()


