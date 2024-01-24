# ESC201 Exercise Week 2, Chiara Zanatta
# Newton's Method on Kepler's Equation

# import libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# define the eccentricity -> How elliptical is the ellipse?
e = 0.5
# define the semi-major axis
a = 1


def fKepler(E, M, e):
	"""Define the Kepler Equation, so that is equal 0"""
	return E - e * math.sin(E) - M


def fKeplerPrime(E, e):
	"""Define the prime of the Kepler's Equation"""
	return 1 - e * math.cos(E)


def newton(fKepler, fKeplerPrime, M, e, Estart):
	"""Solve Kepler's Equation with Newton's Method"""
	E = Estart  # set an initial E
	tol = 0.001  # set a tolerance
	dx = tol + 1  # set dx over tol to start the while loop
	while abs(dx) > tol:
		dx = -fKepler(E, M, e) / fKeplerPrime(E, e)  # calculate the distance from the sun
		Enew = E + dx  # use the calculated E to define the new E and continue till stop factor is reached
		E = Enew
	return E


def init():
	"""Initialize the Plot"""
	# plt.cla()
	ax.set_xlim(-2.1, 1.1)
	ax.set_ylim(-1.1, 1.1)
	return ln, # with a comma to unpack the tuple that is otherwise returned


def update(frame):
	"""Define the different frames for the animation"""
	global xdata, ydata
	M = frame
	E = newton(fKepler, fKeplerPrime, M, e, M)
	x = a * math.cos(E) - a * e
	y = a * math.sqrt(1 - e * e) * math.sin(E)

	xdata.append(x)
	ydata.append(y)

	ln.set_data(xdata, ydata)
	return ln,


fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro', animated=True)

ani = FuncAnimation(fig, update, frames=np.linspace(0, 2 * np.pi, 128),  # define start as 0 and stop as 2*pi to get a full 360° of the circle
					init_func=init, blit=True, interval=50, repeat=False)

ani.save("kepler.mp4", writer="ffmpeg", dpi=250)
plt.show()
