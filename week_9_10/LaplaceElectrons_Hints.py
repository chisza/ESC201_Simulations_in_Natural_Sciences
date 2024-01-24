#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 8 09:42:39 2019

@author: Stefan
"""

import numpy as np


def placePlate(p1, p2, R, U, potential, ax):
    """
    Parameters
    ----------
    p1 : tuple (x,y)
        Starting point of line segment
    p2 : tuple (x,y)
        End point of line segment
    R : (N,N) matrix
        Defines the boundary conditions, zero entries at the boundary
    U : (N,N) matrix
        Potential matrix.
    potential : float
        Potential value of the plate.
    ax : axis object
        Used for plotting the plate

    Returns
    -------
    None.

    """
    x1, y1 = p1
    x2, y2 = p2
    if x2 < x1:
        t = x2
        x2 = x1
        x1 = t
        t = y2
        y2 = y1
        y1 = t
    a = (y2 - y1) / (x2 - x1)
    j1 = int(x1 / DELTA_X)
    j2 = int(x2 / DELTA_X)
    l1 = int(y1 / DELTA_Y)
    l2 = int(y2 / DELTA_Y)

    n = max(j2-j1+1, l2-l1+1)
    for i in range(n+1):
        x = x1 + i*(x2 - x1)/n
        y = y1 + a*(x - x1)
        j = int(x / DELTA_X)
        l = int(y / DELTA_Y)
        R[l, j] = 0
        U[l, j] = potential
    ax.plot([x1, x2], [y1, y2])


def generateElectrons(n):
    y = np.linspace(0.7*BOX_HEIGHT, 0.9*BOX_HEIGHT, n)
    x = np.zeros_like(y)
    # Take a random angle phi in the range -pi/2 to pi/2,
    # either uniform or normal distributed
    # vx = 1e6 * cos(phi)
    # vy = 1e6 * sin(phi)
    vx = 1e6*np.ones_like(x)
    vy = np.zeros_like(y)
    # can be vectors -> can work with an array of electrons
    return (x, y, vx, vy)

def coordToIndex(x, y):
    # delta_x and delta_y -> grid spacing 1 cm /100 grid cells
    j = np.array(x / DELTA_X, dtype='int')
    l = np.array(y / DELTA_Y, dtype='int')
    return (j, l)
    
# How to calculate t and u and calculate the acceleration
# -------------------------------------------------------
#
#    j, l = coordToIndex(x, y)
#
#    Make sure j,l are inside grid    
#    L, J = U.shape
#    j = np.maximum(j, np.zeros_like(j))
#    j = np.minimum(j, (J-2)*np.ones_like(j))
#    l = np.maximum(l, np.zeros_like(l))
#    l = np.minimum(l, (L-2)*np.ones_like(l))
#    
#    t = (x - j*DELTA_X) / DELTA_X
#    u = (y - l*DELTA_Y) / DELTA_Y
#    are the phi values
#    U1 = U[l,j]
#    U2 = ...
#    U3 = ...
#    U4 = ...
#    ax = -(e/me)/DELTA_X*((1-u)*(U2-U1) + u*(U4-U3))
#    ay = ...


# Electron charge, mass and time units
e = -1.6e-19 # Coulomb
me = 9.11e-31 # kg
h = 10e-12 # s, time step
tEnd = 6e-9 # s
    
    
