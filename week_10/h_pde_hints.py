# 1) Using simple Python loop:
#
# def timeStep1(rho):
# rhoNew = np.zeros_like(rho)
# N = len(rho)
# for j in range(N):
#   jp1 = (j + 1) % N
#   jm1 = (j - 1) % N
#   rhoNew(j) = ... method ...
# return rhoNew
#
#
# 2) Do not use a for loop!
#
# def timeStep2(rho):
# np.roll(rho, 1):  Roll the array by one position to the right [1,2,3] -> [3,1,2]
# np.roll(rho, -1):  Roll the array by one position to the left  [1,2,3] ->j [2,3,1]
#
# rhoJplus1 = np.roll(rho, ????) Think twice!
# rhoJminus1 = ...
#
# rhoNew = 0.5*rhoJplus1 - 0.5*rhoJminus1
# return rhoNew
#
# 3) Optional: Use convolve (see SOR)
#
# def timeStep3(rho):
# a_j = c * a_{j-1} + b * a_{j+1}
#
# W = np.array([b 0 c])
# aNew = convolve(a, W, mode="periodic")
#
#
# 4) Initialziation:
#
# N = 400
# rho = np.zeros(N)
#
# for j in range(int(N/2)-10, int(N/2)+10):
#   rho[j] = 1
