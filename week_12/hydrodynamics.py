import numpy as np
import matplotlib.pyplot as plt

# Something is wrong with it, and I don't know what it is

N = 1000
degrees_of_freedom = 3
gamma = degrees_of_freedom + 2 / degrees_of_freedom
dx = 0.1

def initial_condition():
	"""Initial conditions for shock tube problems"""
	U = np.empty((N, 3))
	rho = 1
	u = 0
	e = 1e-5
	E = 0.5 * rho * u**2 + rho * e
	U[:int(N/2), 0] = rho
	U[:int(N/2), 1] = rho * u
	U[:int(N/2), 2] = E
	rho = 4
	E = 0.5 * rho * u ** 2 + rho * e
	U[int(N / 2):, 0] = rho
	U[int(N / 2):, 1] = rho * u
	U[int(N / 2):, 2] = E
	return U

def inital_Sedov_Taylor():
	"""Initial conditions for shock blast problems"""
	U = np.empty((N, 3))
	rho = 1
	u = 0
	e = 1e-5
	E = 0.5 * rho * u ** 2 + rho * e
	U[:, 0] = rho
	U[:, 1] = rho * u
	U[:, 2] = E
	e = 1
	E = 0.5 * rho * u ** 2 + rho * e
	U[int(N/2), 2] = E
	return U

def F(U):
	rho = U[:, 0]
	rho_u = U[:, 1]
	E = U[:, 2]
	u = rho_u / rho
	P = (gamma - 1) * (E - 0.5 * rho_u * u)
	F0 = rho_u
	F1 = rho_u * u + P
	F2 = u * (E + P)
	f = np.array([F0, F1, F2])
	return f.T

def method_a(U, dt):
	"""Original LAX method"""
	f = F(U)
	f_left = np.roll(f, -1, axis=0)
	f_right = np.roll(f, 1, axis=0)
	Unew = 0.5 * (np.roll(U, -1, axis=0) + np.roll(U, 1, axis=0))
	Unew -= dt/(2*dx) * (f_right - f_left)
	return Unew

def D(U):
	rho = U[:, 0]
	rho_u = U[:, 1]
	E = U[:, 2]
	u = rho_u / rho
	P = (gamma - 1) * (E - 0.5 * rho_u * u)
	cs = np.sqrt(gamma * P / rho)
	d = np.abs(u) + cs
	return d

def F_left(U):
	U_left = np.roll( U, 1, axis=0)
	DL = D(U_left)
	DR = D(U)
	Dmax = np.maximum(DL, DR)
	dt = dx / np.max(Dmax)
	Dmax = np.array([Dmax, Dmax, Dmax]).T
	f = 0.5 * (F(U) + F(U_left)) - 0.5 * Dmax * (U - U_left)
	return (f, dt)

def method_b(U):
	"""LAX Friedrich Rieman Solver"""
	fl2, dt1 = F_left(U)
	fr2, dt2 = F_left(np.roll(U, -1, axis=0))
	dt = 0.9 * min(dt1, dt2)
	Unew = U - dt/dx * (fr2 - fl2)
	return Unew

def method_c(U):
	"""LAX Friedrich Riemann Solver, but with half-steps for approximation"""
	# prediction step
	fl2, dt1 = F_left(U)
	fr2, dt2 = F_left(np.roll(U, -1, axis=0))
	dt = 0.99 * min(dt1, dt2)
	Ustar = U - 0.5 * dt / dx * (fr2 - fl2)

	#update step
	fl2, dt1 = F_left(Ustar)
	fr2, dt2 = F_left(np.roll(Ustar, -1, axis=0))
	dt = 0.99 * min(dt1, dt2)
	U_new = U - dt / dx * (fr2 - fl2)
	return U_new

#U = initial_condition() # Shock Tube
U = inital_Sedov_Taylor() # Shock Blast

fig, ax = plt.subplots(3, 1, constrained_layout = True)
fig.suptitle("Shock Blast, method C")

ax[0].plot(U[:, 0])
ax[0].set_title("rho")

ax[1].set_title("rho*u")
ax[1].plot(U[:, 1])

ax[2].set_title("E")
ax[2].plot(U[:, 2])

for n in range (2000):
	U = method_a(U, 0.1)
	#U = method_b(U)
	#U = method_c(U)

	if n % 200 == 0:
		ax[0].plot(U[:, 0])
		ax[1].plot(U[:, 1])
		ax[2].plot(U[:, 2])

fig.savefig("Schock_blast_c.pdf")
fig.show()

plt.close()