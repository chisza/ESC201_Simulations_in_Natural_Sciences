import numpy as np
import matplotlib.pyplot as plt

# with arrays


def LogEq(a, x):
	"""Define the logistic equation"""
	for k in range(1000):
		x = a * x * (1 - x)
	return x


# get the values
aListLen = 1000
xListLen = 1000

aList = np.linspace(0.0, 4, aListLen)
xList = np.linspace(0.0, 1, xListLen)

A, X = np.meshgrid(aList, xList)

# apply the function
feigenbaum = LogEq(A, X)

# calculate the number of values that are should be plotted
nplot = int(len(A)/2)

# generate the plot
for i in range(len(A)):
	plt.scatter(aList[i] * np.ones_like(feigenbaum[i][nplot:]),
				[item[i] for item in feigenbaum][nplot:], s = 0.1, c = 'k')

plt.show()
