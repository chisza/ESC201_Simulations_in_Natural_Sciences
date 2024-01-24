# ESC 201 Exercise 3
# Plot the Feigenbaum Diagram with the Logistic Equation

import numpy as np
import matplotlib.pyplot as plt

# with a loop
def LogEq(a, xList):
	"""Define the logistic equation"""

	# start the list
	xinfList = []

	# Iterate over the list -> reimplement the found value every time
	for x in xList:
		xn = x
		# Solve the logistic equation a thousand times
		for n in range(1000):
			xn = a * xn * (1-xn)
		# create a list with all the calculated values
		print(xn)
		xinfList.append(xn)
		print(xinfList)
	return xinfList

# define the length of a and x
aListLen = 100
xListLen = 100

# create an array with evently spaced numbers in the allowed boundary
aList = np.linspace(0, 4, aListLen)
xList = np.linspace(0, 1, xListLen)

# run the function and fill the array with values for the plot, then plot
# the values for each a
for a in aList:
	xinfList = LogEq(a, xList)

	for i in range((len(xinfList)//2), len(xinfList)):
		plt.scatter(a, xinfList[i], s = 0.1, c = 'k')
plt.show()

# 0 is displayed on the plot because it is also a solution

