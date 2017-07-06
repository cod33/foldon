#!/usr/bin/env python
#generates error plot for epsilon values. 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

todos = np.load('todos.npy')

#loads the individual simulations. 0 = time, 1 = total, 2 = A, 3 = B, 4 = C
#the function which generated todos concatenated all the RMSDs. this takes them all and takes the minimum along the 0 axis (number of RMSD calculations). You're left  with a (5,10001) array for each with the above values. 
bound_25 = np.min(todos[0:6], axis=0)
unbound_25 = np.min(todos[6:12], axis=0)
bound_27 = np.min(todos[12:18], axis=0)
unbound_27 = np.min(todos[18:24], axis=0)
bound_29 = np.min(todos[24:30], axis=0)
unbound_29 = np.min(todos[30:36], axis=0)
bound_31 = np.min(todos[36:42], axis=0)
unbound_31 = np.min(todos[42:48], axis=0)
bound_33 = np.min(todos[48:54], axis=0)
unbound_33 = np.min(todos[54:60], axis=0)
bound_35 = np.min(todos[60:66], axis=0)
unbound_35 = np.min(todos[66:72], axis=0)
bound_37 = np.min(todos[72:78], axis=0)
unbound_37 = np.min(todos[78:84], axis=0)
bound_39 = np.min(todos[84:90], axis=0)
unbound_39 = np.min(todos[90:96], axis=0)

#this takes the above values and puts all three chains A,B,C together in a tuple. You can plot individually, or, as the script has already done, plot the average for all three.

pugh= [(bound_25[2],bound_25[3],bound_25[4]),
       (bound_27[2],bound_27[3],bound_27[4]),
       (bound_29[2],bound_29[3],bound_29[4]),
       (bound_31[2],bound_31[3],bound_31[4]),
       (bound_33[2],bound_33[3],bound_33[4]),
       (bound_35[2],bound_35[3],bound_35[4]),
       (bound_37[2],bound_37[3],bound_37[4]),
       (bound_39[2],bound_39[3],bound_39[4]),
       (unbound_25[2],unbound_25[3],unbound_25[4]),
       (unbound_27[2],unbound_27[3],unbound_27[4]),
       (unbound_29[2],unbound_29[3],unbound_29[4]),
       (unbound_31[2],unbound_31[3],unbound_31[4]),
       (unbound_33[2],unbound_33[3],unbound_33[4]),
       (unbound_35[2],unbound_35[3],unbound_35[4]),
       (unbound_37[2],unbound_37[3],unbound_37[4]),
       (unbound_39[2],unbound_39[3],unbound_39[4])]

epsilons = [0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39,
            0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39]
means = np.empty((16,3)) 
stds = np.empty((16,3)) 

for idx,item in enumerate(pugh):
    A = item[0]
    B = item[1]
    C = item[2]
    means[idx]=np.mean(A),np.mean(B),np.mean(C)
    stds[idx]=(np.std(A)),(np.std(B)),(np.std(C))


epsilons_b=[]
for item in epsilons:
    epsilons_b.append(item - 0.0033)
plt.errorbar(epsilons_b,means[:,1],yerr=stds[:,1],fmt='o',label='Chain A')
plt.errorbar(epsilons,means[:,0],yerr=stds[:,0],fmt='o',label='Chain B')
epsilons_c = []
for item in epsilons:
    epsilons_c.append(item + 0.0033)
plt.errorbar(epsilons_c,means[:,2],yerr=stds[:,2],fmt='o',label='Chain C')
plt.legend()
plt.xlabel(r'$\epsilon$')
plt.ylabel(r'RMSD $\AA$')
plt.title('68% confidence of monomers, Langevin, standard dev')
plt.xticks(epsilons,epsilons)
plt.savefig('new_std_plot_indi_mono')
plt.show()

