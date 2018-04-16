import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_histogram(hist, edges, ax):
    ys = np.repeat(hist,2)
    xs = np.hstack((edges[0], np.repeat(edges[1:-1],2),edges[-1]))
    ax.plot(xs,ys)

SCF_bound = np.hstack((np.load('SCF_bound_1.npy'),np.load('SCF_bound_2.npy'),np.load('SCF_bound_3.npy'),np.load('SCF_bound_4.npy'),np.load('SCF_bound_5.npy')))
SCF_unbound = np.hstack((np.load('SCF_unbound_1.npy'),np.load('SCF_unbound_2.npy'),np.load('SCF_unbound_3.npy'),np.load('SCF_unbound_4.npy'),np.load('SCF_unbound_5.npy')))
WT_bound = np.hstack((np.load('WT_bound_1.npy'),np.load('WT_bound_2.npy'),np.load('WT_bound_3.npy'),np.load('WT_bound_4.npy'),np.load('WT_bound_5.npy')))
WT_unbound = np.hstack((np.load('WT_unbound_1.npy'),np.load('WT_unbound_2.npy'),np.load('WT_unbound_3.npy'),np.load('WT_unbound_4.npy'),np.load('WT_unbound_5.npy')))

todos=[SCF_bound,SCF_unbound,WT_bound,WT_unbound]
todosstr = ['SCF_bound','SCF_unbound','WT_bound','WT_unbound']

f,((ax1,ax2,ax3,ax4)) = plt.subplots(4,sharex='col',sharey='row')
axes = [ax1,ax2,ax3,ax4]
ax = plt.gca()

for idx,arr in enumerate(todos):
    edges = np.arange(0,1300,20)
    hist,_ = np.histogram(arr,bins=edges,density=True)
    plot_histogram(hist,edges,axes[idx])
    axes[idx].set_title(str(todosstr[idx]))

plt.savefig('hist_analysis')
