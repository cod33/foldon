#!/usr/bin/env python
#Corinn 17.5.25
#this should convert output data to numpy array, calculate RMSDs six times, take the minimum. 

import MDAnalysis as md
import MDAnalysis.analysis.rms as rms
import numpy as np
import matplotlib.pyplot as plt
import os

file_coords = ['MOVIE/movie.000000000.pdb','testout.xtc']
todos= []
ref = md.Universe(os.path.join(file_coords[0]))
coords = md.Universe(os.path.join(file_coords[0]),os.path.join(file_coords[1]))
rmsd = rms.RMSD(coords,ref)
rmsd_A = rms.RMSD(coords,ref,select='bynum 1-70')
rmsd_B = rms.RMSD(coords,ref,select='bynum 91-161')
rmsd_C = rms.RMSD(coords,ref,select='bynum 182-252')
for R in [rmsd,rmsd_A,rmsd_B,rmsd_C]:
    R.run()
full_arr = np.vstack((rmsd.rmsd[:,1],rmsd.rmsd[:,2],rmsd_A.rmsd[:,2],rmsd_B.rmsd[:,2],rmsd_C.rmsd[:,2]))
print(full_arr.shape)
#todos.append(full_arr)
#todoss = np.array(todos)
np.save('todos.npy',full_arr)
#ok so todos has shape (96,5,10001). it goes bound, unbound by increments of six. 
#bound_25 = np.min(todos[0:6], axis=0)
#unbound_25 = np.min(todos[6:12], axis=0)
#bound_27 = np.min(todos[12:18], axis=0)
#unbound_27 = np.min(todos[18:24], axis=0)
#bound_29 = np.min(todos[24:30], axis=0)
#unbound_29 = np.min(todos[30:36], axis=0)
#bound_31 = np.min(todos[36:42], axis=0)
#unbound_31 = np.min(todos[42:48], axis=0)
#bound_33 = np.min(todos[48:54], axis=0)
#unbound_33 = np.min(todos[54:60], axis=0)
#bound_35 = np.min(todos[60:66], axis=0)
#unbound_35 = np.min(todos[66:72], axis=0)
#bound_37 = np.min(todos[72:78], axis=0)
#unbound_37 = np.min(todos[78:84], axis=0)
#bound_39 = np.min(todos[84:90], axis=0)
#unbound_39 = np.min(todos[90:96], axis=0)
#
#
#ugh = [bound_25,unbound_25,bound_27,unbound_27,bound_29,unbound_29,bound_31,unbound_31,bound_33,unbound_33,bound_35,unbound_35,bound_37,unbound_37,bound_39,unbound_39]
#
##for item in ugh:
#    plt.plot(item[0],item[1],label="tots")
#    plt.plot(item[0],item[2],label="A")
#    plt.plot(item[0],item[3],label="B")
#    plt.plot(item[0],item[4],label="C")
#    plt.legend()
#    plt.savefig(item)
