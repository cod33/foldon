#!/usr/bin/env python
#Corinn whatever it's june 2017
#Ali must have helped you with this cause you don't know what's happening (3.12.18)

import numpy as np

arr = np.loadtxt('system.go.parameters',skiprows=1)
length = len(arr[:,0])
A=0
B=0
C=0
AB = 0
AC = 0
BC = 0
sigmas = []
#sigmas represent equilibrium distances probably?
for i in xrange(length):
    one = arr[i,0]
    two = arr[i,2]
    sigma = arr[i,4]
    if one == two and sigma not in sigmas:
        if one == 1:
            A += 1
        if one == 2:
            B += 1
        if one ==3:
            C += 1
    if one == 1 and two == 2 and sigma not in sigmas:
        AB += 1
    if one ==1 and two == 3  and sigma not in sigmas:
        AC += 1
    if one ==2 and two ==3 and sigma not in sigmas:
        BC += 1
    sigmas.append(sigma)

print('A', A, 'B', B, 'C', C)
print('AB',AB,'AC',AC,'BC',BC)

