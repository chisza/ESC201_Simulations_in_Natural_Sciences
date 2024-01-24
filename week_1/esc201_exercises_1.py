# ESC201 Exercise Week 1, Chiara Zanatta

# Define the input function
# Function is 0 = x^x - 100
def func(x):
	"""Define a Function"""
	return x ** x - 100


# Do the calculation with the bisection method
# Define the section of interest with values a and b
# Condition: a has to be left, b has to be right
# f(a) has to be negative, f(b) has to be positive
def bisection_method(a, b):
	"""Use the Bisection Method"""

	# check first, that the conditions are correct
	if func(a) > func(b):
		print("Not the right starting conditions for a or b.")
		return

# Define c
	c = 1 / 2 * (a + b)

# Have a condition for a stop
	while abs((a-b)/c) > 10e-6:
		# find the middle point of a and b
		c = 1 / 2 * (a + b)

		# check if f(c) > 0 -> yes: set b=c, no: set a=c
		if func(c) > 0:
			b = c
		if func(c) < 0:
			a = c

	return c


# enter values that useful, 1^1 - 100 = -99, 4^4 - 100 = 156
root = bisection_method(1, 4)
print("The value of root is: ", root)
