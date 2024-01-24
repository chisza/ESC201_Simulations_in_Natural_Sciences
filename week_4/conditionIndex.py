# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 09:18:32 2022

@author: Stefan
"""

import numpy as np

a = np.array([1, 2, 3])
c = np.zeros_like(a, dtype=bool)

print("a=", a)
print("c=", c)

print("a[c]=", a[c])
print("~c=", ~c)
print("a[~c]=", a[~c])

c = a > 1
print("a > 1:", c)
print("a[a > 1]:", a[c])

a2 = 2 * a[c]
print("2 * a[c]=", a2)

x = np.linspace(0, 2, 3)
y = np.linspace(0, -2, 3)

X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

print("Z=", Z)
c = abs(Z) > 1
print("abs(Z) > 1=", c)

import numpy as np
from imageio import imwrite

m = 480
n = 320

x = np.linspace(-2, 1, num=m).reshape((1, m))
y = np.linspace(-1, 1, num=n).reshape((n, 1))
C = np.tile(x, (n, 1)) + 1j * np.tile(y, (1, m))

Z = np.zeros((n, m), dtype=complex)
M = np.full((n, m), True, dtype=bool)
for i in range(100):
	Z[M] = Z[M] * Z[M] + C[M]
	M[np.abs(Z) > 2] = False

imwrite('mandelbrot.png', np.uint8(np.flipud(1 - M) * 255))


