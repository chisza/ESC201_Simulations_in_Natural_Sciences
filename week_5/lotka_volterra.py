# Exercises week 5
# Solve the Lotka-Volterra Model with different ODE Methods
import numpy as np
import matplotlib.pyplot as plt


def odeSolver(t0, y0, dfFunc, h, nSteps, solverStepFunc):

	"""This is a general ODE solver that takes the
	derivative df/dt (dfFunc) and the algorithm for one time
	step (solverStepFunc) as function arguments.
	t0 = Initial time
	y0 = Initial function values (array)
	nSteps = Total number of integration steps
	solverStepFunc = Euler Method, Midpoint RK or RK4
	"""
	yn = y0
	tn = t0
	tlist = [t0]
	ylist = [y0]

	# Solve the specified method for the amount of times and then
	# and return the time steps in one list, and the calculated
	# yn in another list
	for n in range(nSteps):
		yn1 = solverStepFunc(tn, yn, h, dfFunc)
		tn1 = tn + h
		tlist.append(tn1)
		ylist.append(yn1)
		tn = tn1
		yn = yn1
	return (np.array(tlist), np.array(ylist))


def eulerStep(tn, yn, h, dfdt):
	"""Solves the Forward Euler Method"""
	yn1 = yn + h * dfdt(tn, yn)
	return yn1


def MidPointRK2Step(tn, yn, h, dfdt):
	"""Solves the Runge Kutta 2 Method"""
	yn1 = yn + h * dfdt(tn + h/2, yn + h/2 * dfdt(tn, yn))
	return yn1


def RK4Step(tn, yn, h, dfdt):
	"""Solves the Runge Kutta 4 Method"""
	k1 = h * dfdt(tn, yn)
	k2 = h * dfdt(tn + h/2, yn + k1/2)
	k3 = h * dfdt(tn + h/2, yn + k2/2)
	k4 = h * dfdt(tn + h, yn + k3)
	yn1 = yn + k1/6 + k2/3 + k3/3 + k4/6
	return yn1


def LotkaVolterra(t, y):
	"""Implements the Lotka Volterra System dm/dt and dfox/dt
	where y=(m,fox), dy/dt=(dm/dt, dfox/dt)"""
	prey = int(y[0])
	predator = int(y[1])
	K_m = 2
	K_mf = 0.02
	K_f = 1.06
	K_fm = 0.01
	dmdt = (K_m * prey - K_mf * prey * predator)
	dfoxdt = ((-K_f * predator) + K_fm * predator * prey)
	return np.array([dmdt, dfoxdt]) # = dy/dt

# with these values also euler oscillates nicely
h = 0.001
nSteps = 60000
t0 = 0
#y0 = np.array([100, 15])

# define different starting conditions
starting_conditions_prey = [1, 15, 100, 600, 100]
starting_conditions_predator = [16, 100, 44, 250, 15]

# the methods in a list
method_list = [eulerStep, MidPointRK2Step, RK4Step]

# do the for loop twice, once to get the plots for the population
# versus time, once for the phase plots (at some point, figure out
# how to do both plots in the same loop)

# loop for the population vs time
for i in range(len(starting_conditions_prey)):
	y0 = np.array([starting_conditions_predator[i], starting_conditions_predator[i]])

	# create figure
	oscillations = plt.figure(figsize=(10, 7))

	# setting values to rows and column variables
	rows = 2
	columns = 2

	for j in range(len(method_list)):
		t, y = odeSolver(t0, y0, LotkaVolterra, h, nSteps, method_list[j])

		# extract the function name of the method
		title = str(method_list[j])
		title = title.split(" ")
		title = title[1]

		# fill the plots
		oscillations.add_subplot(rows, columns, j+1)
		plt.title(title) # get the function name as title for the plot
		plt.plot(t, y[:, 0], label = "Prey")  # population versus time
		plt.plot(t, y[:, 1], label = "Predator")
		plt.legend()

	oscillations.savefig("Lotka_Volterra_oscillations" + str(starting_conditions_prey[i]) + "_" +
						str(starting_conditions_predator[i]) + ".png")

# loop for the phase plots
for i in range(len(starting_conditions_prey)):
	y0 = np.array([starting_conditions_predator[i], starting_conditions_predator[i]])

	# create figure
	phases = plt.figure(figsize=(10, 7))

	# setting values to rows and column variables
	rows = 2
	columns = 2

	for j in range(len(method_list)):
		t, y = odeSolver(t0, y0, LotkaVolterra, h, nSteps, method_list[j])

		title = str(method_list[j])
		title = title.split(" ")
		title = title[1]

		phases.add_subplot(rows, columns, j+1)
		plt.title(title)
		plt.plot(y[:, 0], y[:, 1])
		plt.xlabel("Prey")
		plt.ylabel("Predator")

	phases.savefig("Lotka_Volterra_phases" + str(starting_conditions_prey[i]) + "_" +
						str(starting_conditions_predator[i]) + ".png")
