#!/usr/bin/env python
#should plot one round of native contacts. How do I not have this, still

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

arr = np.load('nc_arr.npy')

xs = range(len(arr[:,0]))
step = 250
plt.plot(xs,arr[:,3],alpha=0.5)
plt.plot(xs[step+1:],[arr[:,3][i:i+step].mean() for i in range(int(len(arr[:,3]))-int(step+1))])
print(i,i+step)
print(len(arr[:,0]-(step+1)))
plt.title('native contacts bound on')
plt.xlabel('time (uS)')
plt.ylabel('number of native contacts')
plt.savefig('nc_bound_on.png')
