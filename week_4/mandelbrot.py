# import libraries
import numpy as np
import matplotlib.pyplot as plt

# Mandelbrot

# define the resolution
real = 2000
imaginary = 2000

# define arrays
x_real = np.linspace(-2, 2, real)
y_imaginary = np.linspace(-2, 2, imaginary)

# get the arrays for the axes of the grid
X_real, Y_imaginary = np.meshgrid(x_real, y_imaginary)
print(X_real, Y_imaginary)

# get a complex array that represents the coordinates we want to map
C = X_real + 1j * Y_imaginary
print(C)

# get an array with the values that are to be updated and checked
Z = np.zeros((real, imaginary), dtype = complex)

# counter
counter = np.zeros((real, imaginary))

# make the calculations
maxiter =100 # maximum number of iterations

k = 0
while k < maxiter:
	# set the boundary
	boundary = abs(Z) < np.maximum(2*np.ones_like(C), abs(C)) # gives a matrix with true and false, value does not get overwritten
	# do the calculation
	Z[boundary] = Z[boundary] * Z[boundary] + C[boundary] # only use the ones that fullfill the condition
	counter[boundary] += 1 # add one to each which are still considered in the iteration
	# add to the iteration
	k += 1
print(Z)
print(counter)

plt.imshow(counter)
plt.show()