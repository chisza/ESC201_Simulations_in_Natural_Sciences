import numpy as np
import matplotlib.pyplot as plt

# define the resolution
real = 2000
imaginary = 2000

# define arrays
x_real = np.linspace(-3, 2, real)
y_imaginary = np.linspace(-2, 2, imaginary)

# get the arrays for the axes of the grid
X_real, Y_imaginary = np.meshgrid(x_real, y_imaginary)
print(X_real, Y_imaginary)

# get a complex array that represents the coordinates we want to map
C = X_real + 1j * Y_imaginary
print(C)

# array to keep track of which grid points are bounded
counter = np.zeros((real, imaginary))

# maximal number of iterations
num_iter = 200

# the complex part that should be looked at, transformed to an array
#a = 0.285 + 0.01 * 1j
a = -0.744 + 0.148 * 1j

a = np.full((real, imaginary), a)

for j in range(num_iter):
	# gives a matrix with true and false, value does not get overwritten
	boundary = abs(C) < np.maximum(2*np.ones_like(C), np.abs(a))

	# the calculation
	C[boundary] = C[boundary] ** 2 + a[boundary]

	# count the calculations that were done with the considered
	# values
	counter[boundary] = counter[boundary] + 1

plt.imshow(counter)
plt.xlim(500, 1750)
plt.ylim(500, 1750)
plt.show()