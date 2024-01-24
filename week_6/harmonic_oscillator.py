from math import sin
import matplotlib.pyplot as plt
import numpy as np

def odeHO(p,q):
	"""
	Harmonic oscillator
	Do not want to look at that in this exercise!
	"""
	dpdt = -q
	dqdt = p
	return (dpdt, dqdt)

def odePendulum(p,q):
	"""Pendulum
	Is the derivative of H = 1/2p^2 -epsilon*cos*q"""
	dpdt = - epsilon*sin(q)
	dqdt = p
	return (dpdt, dqdt)

def leapFrog(p0, q0, h, odeSystem):
	"""Leap Frog Method to solve ODE"""
	q12 = q0 + 0.5*h*p0         # first drift
	dp, dq = odeSystem(p0, q12)
	p1 = p0 + h*dp              # kick
	q1 = q12 + 0.5*h*p1         # second drift
	return (p1, q1)

def euler(p, q, h, odeSystem):
	"""Forward Euler Method to solve ODE"""
	dp, dq = odeSystem(p, q)
	p += h*dp
	q += h*dq
	return (p, q)

def rk2MidPoint(p, q, h, odeSystem):
	"""Runge-Kutta 2 Method to solve ODE"""
	dp = odeSystem(p + h/2., q + h/2. * float(odeSystem(p, q)[0]))
	dq = odeSystem(p + h / 2., q + h / 2. * float(odeSystem(p, q)[1]))
	p += h*dp[0]
	q += h*dq[1]
	return (p, q)

def rk4(p, q, h, odeSystem):
	"""Runge-Kutta 4 Method to solve ODE"""
	k1_p = h * odeSystem(p, q)[0]
	k1_q = h * odeSystem(p, q)[1]
	k2_p = h * odeSystem(p + h / 2, q + k1_p / 2)[0]
	k2_q = h * odeSystem(p + h / 2, q + k1_q / 2)[1]
	k3_p = h * odeSystem(p + h / 2, q + k2_p / 2)[0]
	k3_q = h * odeSystem(p + h / 2, q + k2_q / 2)[1]
	k4_p = h * odeSystem(p + h, q + k3_p)[0]
	k4_q = h * odeSystem(p + h, q + k3_q)[1]

	#p1 = p + (k1_p/6 + k2_p/3 + k3_p/3 + k4_p/6)
	#q1 = q + (k1_q/6 + k2_q/3 + k3_q/3 + k4_q/6)

	p += (k1_p/6 + k2_p/3 + k3_p/3 + k4_p/6)
	q += (k1_q/6 + k2_q/3 + k3_q/3 + k4_q/6)

	return (p, q)


def odeSolver(p0, q0, dfFunction, h, nSteps, odeSystem):
	"""A solver for differential equations"""

	# set starting conditions
	pn = p0 # starting momentum
	qn = q0 # starting angle
	plist = [p0] # list of the calculated momentum
	qlist = [q0] # list of the calculated angle

	# Solve the system for a specific number of times
	# get all the calculated values in the lists
	for step in range(nSteps):
		pn1, qn1 = odeSystem(pn, qn, h, dfFunction)
		plist.append(pn1)
		qlist.append(qn1)
		pn, qn = pn1, qn1
	return (np.array(plist), np.array(qlist))


# Starting conditions
h = 0.1
nSteps = 200
p0 = 1
#q0 = 1
epsilon = 1

for i in range(-25, 25):
	q0 = i * 0.1
	p_calculated, q_calculated = odeSolver(p0, q0, odeHO, h, nSteps, euler)
	plt.plot(q_calculated, p_calculated, label = f"{q0:.2f}")


plt.legend()
plt.show()




#initial conditions
# HO: p=0, q=1, q=2, q=3
#
# Pendulum: q is angle, range is [-pi,pi]
# For circulation you need an additional momentum p /= 0
#
# axs[0,0]: phase plot HO using leap frog
# axs[0,1]: phase plot HO using euler or rk2 or rk4
# axs[1,0]: phase plot pendulum using leap frog
# axs[1,1]: phase plot pendulum using euler or rk2 or rk4


