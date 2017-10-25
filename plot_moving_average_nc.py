#!/usr/bin/env python
#should plot one round of native contacts. How do I not have this, still

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

files=['1_nc.npy','2_nc.npy','3_nc.npy','4_nc.npy','5_nc.npy']

f,((ax1,ax2,ax3),(ax4, ax5, ax6)) = plt.subplots(2,3,sharex='col',sharey='row')

nc = np.load(files[0])
xs = np.arange(0, nc.shape[0], 1, dtype=float)
xs *= 10**-4
axs = [ax1,ax2,ax3,ax4,ax5]
step = 250
arrs = []
for idx,fyle in enumerate(files):
    arr = np.load(fyle)
    arrs.append(arr)
for idx,ax in enumerate(axs):
    ax.plot(xs,arrs[idx][:,3],'r',alpha=0.5)
    ax.plot(xs,arrs[idx][:,4],'b',alpha=0.5)
    ax.plot(xs,arrs[idx][:,5],'k',alpha=0.5)
    ax.plot(xs[step+1:],[arrs[idx][:,3][i:i+step].mean() for i in range(int(len(arrs[idx][:,3]))-int(step+1))],'r',label='AB')
    ax.plot(xs[step+1:],[arrs[idx][:,4][i:i+step].mean() for i in range(int(len(arrs[idx][:,4]))-int(step+1))],'b',label='AC')
    ax.plot(xs[step+1:],[arrs[idx][:,5][i:i+step].mean() for i in range(int(len(arrs[idx][:,5]))-int(step+1))],'k',label='BC')

plt.legend()
plt.title('native contacts bound trimer')
ax4.set_xlabel('time (uS)')
ax4.set_ylabel('number of native contacts')
plt.savefig('nc_bound.png')
