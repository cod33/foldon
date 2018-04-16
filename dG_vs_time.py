import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

bound = ['19_0.74_10_bound.npy','19_0.74_1_bound.npy','19_0.74_2_bound.npy','19_0.74_3_bound.npy','19_0.74_4_bound.npy','19_0.74_5_bound.npy','19_0.74_6_bound.npy','19_0.74_7_bound.npy','19_0.74_8_bound.npy','19_0.74_9_bound.npy','22_0.74_10_bound.npy','22_0.74_1_bound.npy','22_0.74_2_bound.npy','22_0.74_3_bound.npy','22_0.74_4_bound.npy','22_0.74_5_bound.npy','22_0.74_6_bound.npy','22_0.74_7_bound.npy','22_0.74_8_bound.npy','22_0.74_9_bound.npy']
unbound = ['19_0.74_10_unbound.npy','19_0.74_1_unbound.npy','19_0.74_2_unbound.npy','19_0.74_3_unbound.npy','19_0.74_4_unbound.npy','19_0.74_5_unbound.npy','19_0.74_6_unbound.npy','19_0.74_7_unbound.npy','19_0.74_8_unbound.npy','19_0.74_9_unbound.npy','22_0.74_10_unbound.npy','22_0.74_1_unbound.npy','22_0.74_2_unbound.npy','22_0.74_3_unbound.npy','22_0.74_4_unbound.npy','22_0.74_5_unbound.npy','22_0.74_6_unbound.npy','22_0.74_7_unbound.npy','22_0.74_8_unbound.npy','22_0.74_9_unbound.npy']
all_74 = ['19_0.74_10_unbound.npy','19_0.74_1_unbound.npy','19_0.74_2_unbound.npy','19_0.74_3_unbound.npy','19_0.74_4_unbound.npy','19_0.74_5_unbound.npy','19_0.74_6_unbound.npy','19_0.74_7_unbound.npy','19_0.74_8_unbound.npy','19_0.74_9_unbound.npy','22_0.74_10_unbound.npy','22_0.74_1_unbound.npy','22_0.74_2_unbound.npy','22_0.74_3_unbound.npy','22_0.74_4_unbound.npy','22_0.74_5_unbound.npy','22_0.74_6_unbound.npy','22_0.74_7_unbound.npy','22_0.74_8_unbound.npy','22_0.74_9_unbound.npy','19_0.74_10_bound.npy','19_0.74_1_bound.npy','19_0.74_2_bound.npy','19_0.74_3_bound.npy','19_0.74_4_bound.npy','19_0.74_5_bound.npy','19_0.74_6_bound.npy','19_0.74_7_bound.npy','19_0.74_8_bound.npy','19_0.74_9_bound.npy','22_0.74_10_bound.npy','22_0.74_1_bound.npy','22_0.74_2_bound.npy','22_0.74_3_bound.npy','22_0.74_4_bound.npy','22_0.74_5_bound.npy','22_0.74_6_bound.npy','22_0.74_7_bound.npy','22_0.74_8_bound.npy','22_0.74_9_bound.npy']

R = 1.987e-3
T = 298
moles = 2/6.02E23
cub_A = (4/3)*np.pi*(125)**3
Liters= cub_A * 1E-27

conc = moles/Liters 
#each array is (40001,6) ish, where first dimension is number of timesteps. 
#all are 4 us

def is_bound(nativecontactsarray,timepoint):
    nc_cutoff = 50
    AB = int(nativecontactsarray[timepoint,3])
    if AB > nc_cutoff:
        return True
    else:
        return False

def is_unbound(nativecontactsarray,timepoint):
    AB = int(nativecontactsarray[timepoint,3])
    if AB == 0:
        return True
    else:
        return False

statearrays = {}

#f,((ax1,ax2,ax3,ax4,ax5),(ax6,ax7,ax8,ax9,ax10),(ax11,ax12,ax13,ax14,ax15),(ax16,ax17,ax18,ax19,ax20))=plt.subplots(4,5,sharex='col',sharey='row')
#
#axes = [ax1,ax2,ax3,ax4,ax5,
#        ax6,ax7,ax8,ax9,ax10,
#        ax11,ax12,ax13,ax14,ax15,
#        ax16,ax17,ax18,ax19,ax20]

kons = []
koffs = []
erron =[]
erroff = []

for i,path in enumerate(unbound):
    nativecontactsarray = np.load(path)
    statearray = np.ones(nativecontactsarray.shape[0])*-1
    timepoints = nativecontactsarray.shape[0]
    for timepoint in range(timepoints):
        if is_bound(nativecontactsarray,timepoint):
            statearray[timepoint] = 0
        elif is_unbound(nativecontactsarray,timepoint):
            statearray[timepoint] = 1
        else: # we cannot currently say whether the timepoint fits the definition for the boudn or unbound state
            if timepoint > 0:
                statearray[timepoint] = statearray[timepoint-1]
            else:
                statearray[timepoint] = -5
    statearrays[path] = statearray
    #make an array of size (20,40001)
statearrays_all = np.vstack(([statearrays[x] for x in statearrays.keys()]))
nsims = len(unbound)

transitions = statearrays_all[:,1:statearrays_all.shape[1]] - statearrays_all[:,0:statearrays_all.shape[1]-1]
transitions = numpy.concatenate((numpy.zeros((nsims,1)), transitions), axis=1)

#time is in 1E-10 seconds
dG_vs_time = []
dGerr_vs_time = []
timepoints = 40001
for endpoint in range(0, timepoints):
    try:
        k_on = numpy.count_nonzero(transitions[:,:endpoint] == 1)*1e10/\
                 float(numpy.count_nonzero(statearrays_all[:,:endpoint]==0)*conc)
        se_k_on = float((np.count_nonzero(transitions[:,:endpoint]==1)**0.5)*1e10/(np.count_nonzero(statearrays_all[:,:endpoint]==0)*conc))
    except ZeroDivisionError:
        k_on = np.nan
        se_k_on = np.nan
    try:
        k_off = numpy.count_nonzero(transitions[:,:endpoint] == -1)*1e10/\
                 float(numpy.count_nonzero(statearrays_all[:,:endpoint]==1))
        se_k_off = float(np.count_nonzero(transitions[:,:endpoint]==-1)**0.5)*1e10/(np.count_nonzero(statearrays_all[:,:endpoint]==1))
    except ZeroDivisionError:
        k_off=np.nan
        se_k_off = np.nan

    
    try:
        keq = (k_on*conc)/k_off
        dG = R*T*numpy.log(keq)
        eq_error = keq*(((se_k_on/k_on)**2 + (se_k_off/k_off)**2)**0.5)
        se_dG = R*T*(eq_error/(keq))
        dG_vs_time.append(dG)
        dGerr_vs_time.append(se_dG)
    except ZeroDivisionError:
        dG = np.nan
        dG_vs_time.append(dG)
        se_dG = np.nan
        dGerr_vs_time.append(se_dG)

dG_arr = np.array(dG_vs_time)
se_dG_arr = np.array(dGerr_vs_time)
full_arr = np.array((dG_arr,se_dG_arr))
np.save('dG_arr_unbound.npy',full_arr)
print(full_arr.shape)
#step = 250
#xs = range(len(dG_vs_time))
#plt.plot(xs,dG_vs_time,'ro')
#arr = dG_vs_time
#plt.plot(xs[step+1:],[arr[i:i+step].mean() for i in range(int(len(arr))-int(step+1))],'r')
#plt.errorbar(xs,dG_vs_time,yerr=dGerr_vs_time)
#plt.savefig('pls_god.png')
#    else:
#        k_on = float((utb_transition_count)*1e10/(total_time_unbound*conc))
#        se_k_on = float((utb_transition_count**0.5)*1e10/(total_time_unbound*conc))
#        kons.append(k_on)
#        erron.append(se_k_on)
#
#    if total_time_bound == 0:
#        koffs.append(np.nan)
#        erroff.append(np.nan)
#    else:
#        k_off = float((btu_transition_count)*1e10/total_time_bound)
#        se_k_off = float((btu_transition_count**0.5)*1e10/total_time_bound)
#        koffs.append(k_off)
#        erroff.append(se_k_off)
#
#plt.suptitle('state definitions; unbound initial')
#plt.savefig('states_from_unbound.png')
#
#
#
