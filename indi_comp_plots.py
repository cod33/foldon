#!/usr/env/bin python

import numpy as np
import matplotlib.pyplot as plt

todos = np.load('todos.npy')
todosm = np.min(todos[0:6],axis=0)
plt.figure(1)
plt.plot(todosm[0],todosm[1],label="Total")
plt.plot(todosm[0],todosm[2],label="A")
plt.plot(todosm[0],todosm[3],label="B")
plt.plot(todosm[0],todosm[4],label="C")
plt.legend()
plt.savefig('rmsd.png')

nc = np.load('nc_arr.npy')
xs = range(len(nc[:,0]))
plt.figure(2)
plt.plot(xs,nc[:,0],label="A")
plt.plot(xs,nc[:,1],label="B")
plt.plot(xs,nc[:,2],label="C")
plt.plot(xs,nc[:,3],label="AB")
plt.plot(xs,nc[:,4],label="AC")
plt.plot(xs,nc[:,5],label="BC")
plt.legend()
plt.savefig('nc.png')


