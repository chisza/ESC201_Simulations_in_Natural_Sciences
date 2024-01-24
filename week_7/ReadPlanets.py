import pandas as pd
import numpy as np
import astropy.time
import dateutil.parser
import matplotlib.pyplot as plt
import itertools


def readPlanets(filename, N=-1):
	# Loading text files with Pandas
	df = pd.read_csv(filename, sep=',', header=None,
					 names=['name', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'])

	# Data is now in a Pandas dataframe
	# print(df)

	name = np.array(df.loc[:, 'name'])
	m = np.array(df.loc[:, 'm'])
	r = np.array(df.loc[:, 'x':'z'])
	v = np.array(df.loc[:, 'vx':'vz'])

	if N > 0:
		name = name[0:N - 1]
		m = m[0:N - 1]
		r = r[0:N - 1]
		v = v[0:N - 1]

	return (name, r, v, m)


# Loading text files with NumPy
def readPlanetsNp(filename):
	data = np.loadtxt(filename, delimiter=',', unpack=False,
					  dtype={'names': ('planet', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'),
							 'formats': ('S8', float, float, float, float,
										 float, float, float)})

	# print(data)
	name = np.array([d[0] for d in data])
	m = np.array([d[1] for d in data])
	r = np.array([[d[2], d[3], d[4]] for d in data])
	v = np.array([[d[5], d[6], d[7]] for d in data])
	return (name, r, v, m)


name, r, v, m = readPlanets("SolSystData.dat")


# name1, r1, v1, m1 = readPlanetsNp("SolSystData.dat")

# Date for the SolSystData.dat coordinates is 18.03.2008
# which corresponds to 2454543.500000000 julian date.
dt1 = dateutil.parser.parse('18.03.2008')
time1 = astropy.time.Time(dt1)
print(time1.jd)

# Change this date to "today"
dt2 = dateutil.parser.parse('30.10.2023')
time2 = astropy.time.Time(dt2)
print(time2.jd)

# Integrate the solar system from 0 to days
days = time2.jd - time1.jd
print(days)

# G * Msun = k**2
k = 0.01720209895


def accel(k, m, r):
	"""Get the accelerations for each planet"""

	# start an empty array for storing the acceleration
	# same length as planets, as each planet has one acceleration
	a = np.zeros_like(v)
	# loop through the positions of a and assign each an acceleration
	# corresponding to a planet
	for i in range(0, len(a)):
		# loop through the planets but the one currently looked at
		# to get the interactions between all the planet pairs
		for j in range(i + 1, len(a)):
			# calculate the force acting on the planets
			# direction matters!!!!!!!!!!!!!!!!!!!!!!
			# r[j] - r[i] when the force on planet i due to planet j
			dr = r[j] - r[i]
			force_i = ((k ** 2 * m[i] * m[j]) / (np.sqrt(dr.dot(dr)) ** 3)) * dr
			# add / subtract the force
			# the force is a 3d vector, with acceleration in x, y, z direction
			a[i] += force_i / m[i]
			a[j] -= force_i / m[j]

	return a



h = 4

stuff_for_plotting = []

# time steps are 4 -> divide days by h to get the correct
# amount of calculations
for day in range(0, int(days/h)):

	# calculate first drift to give r1/2
	r_step_half = r + h/2 * v

	# calculate acceleration
	acceleration = accel(k, m, r_step_half)

	# calculate the kick
	v_step_one = v + h * acceleration

	# calculate the second drift
	r_step_one = r_step_half + h/2 * v_step_one

	# reassign the values for r and v so that the next step is calculated
	r = r_step_one
	v = v_step_one

	# get a list of all values for plotting
	stuff_for_plotting.append(r_step_one.copy())

fig, axs = plt.subplots(2, 1)
fig.set_size_inches(11.69,8.27)

colors = itertools.cycle(["goldenrod", "darkorange", "cornflowerblue", "red", "khaki", "hotpink", "blue", "violet"])
for r in stuff_for_plotting:
	axs[0].scatter(r[:, 0], r[:, 1], color = ("firebrick", "gray", "tan", "cornflowerblue", "chocolate", "orange", "burlywood", "darkturquoise", "steelblue"), s = .3)
	axs[1].scatter(r[0:5, 0], r[0:5, 1], color = ("firebrick", "gray", "tan", "cornflowerblue", "chocolate"), s = .3)
fig.savefig("solar_system_plot_x_y.png", dpi=300)



# 1. Read Solar system data

# h = 4
# end time = days (see above: time2.jd - time1.jd)
# 2. for time in 0 to end time:
#    3. First drift to give r1/2
#    4. Calculate forces and accels at position r1/2
#    5. Kick: v1 = v0 + ... using accels from step 4.
#    6. Second drift to give r1
#    7. Store r1 for graphics
# 8. Plot all planets (x vs y) (x vs z -> side view)
#
