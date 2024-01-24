#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 14:42:39 2018

https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.convolve.html

@author: Stefan
"""

import numpy as np
from scipy import ndimage, signal
import math
import matplotlib.pyplot as plt


def ex1():
    U = np.ones((5,5))
    W = np.array([[0,0,0],[0,2,0],[0,0,0]])
    Out = np.zeros_like(U)
    ndimage.convolve(U, W, output=Out, mode="constant", cval=0)
    print(U)
    print(W)
    print(Out)


def ex1a():
    U = np.ones((5,5))
    W = np.array([[0,0,0],[0,0,0],[0,0,2]])
    Out = np.zeros_like(U)
    ndimage.convolve(U, W, output=Out, mode="constant", cval=0)
    print(U)
    print(W)
    print(Out)


def ex2():
    U = np.ones((5,5))
    W = np.array([[0,0,0],[1,1,1],[0,0,0]])
    Out = np.zeros_like(U)
    ndimage.convolve(U, W, output=Out, mode="constant", cval=0)
    print(U)
    print(W)
    print(Out)


def ex3():
    U = np.ones((5,5))
    W = np.array([[0,1,0],[0,1,0],[0,1,0]])
    Out = np.zeros_like(U)
    ndimage.convolve(U, W, output=Out, mode="constant", cval=0)
    print(U)
    print(W)
    print(Out)
    

def ex4():
    U = np.ones((5,5))
    W = np.array([[0,1,0],[1,-4,1],[0,1,0]])
    Out = np.zeros_like(U)
    ndimage.convolve(U, W, output=Out, mode="constant", cval=0)
    print(U)
    print(W)
    print(Out)


def ex4corr():
    U = np.ones((5,5))
    W = np.array([[0,1,0],[1,-4,1],[0,1,0]])
    Out = signal.correlate2d(U, W, mode='same', boundary='fill', fillvalue=0)
    print(U)
    print(W)
    print(Out)
    
    #chess board
def ex5():
    C = np.zeros((5, 5), dtype=bool)
    print(C, "\n")
    C[::2, ::2] = True
    print(C, "\n")
    C[1::2, 1::2] = True
    print(C)


height, width = 100, 100

# omgea definition
w = 2 / (1 + math.pi / height)
print(w)

# definition of empty grid
U = np.zeros((height, width))

# definition stencil
stencil = np.array([[0,1,0],[1,-4,1],[0,1,0]])

# initialize boundary conditions -> R -> outside 0, inside w/4
R = np.zeros((height, width))
R.fill(w/4)
R[1, :] = 0
R[-1, :] = 0
R[:, 1] = 0
R[:, -1] = 0

# put center piece it in U / R entries at the same place are 0
U[int(height/2), int(width/3):-int(width/3)] = 1000
R[int(height/2), int(width/3):-int(width/3)] = 0

print(U)

#chess board
C = np.zeros((height, width), dtype=bool)
C[::2, ::2] = True
C[1::2, 1::2] = True
print(C)

#loop
correction_value = 2

while abs(correction_value) > 1:
    Out = np.zeros_like(U)
    ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0) #stencil angewwendet azu fu
    U[C] = U[C] + R[C] * Out[C]

    correction = R[C] * Out[C] # this is the correction that is calculated and added every time, it gets so small over time, that is shows that it converged
    correction_value = np.max(correction)

    ndimage.convolve(U, stencil, output=Out, mode="constant", cval=0) #stencil angewwendet azu fu
    U[~C] = U[~C] + R[~C] * Out[~C]


print(U)

plt.imshow(U, cmap = "plasma")
plt.contour(U, colors = "black")
cs = plt.contour(U, colors = "black")
plt.clabel(cs)
plt.show()
