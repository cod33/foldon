#!usr/bin/env python
# should plot fraction of native contacts for the dimer. maybe eventually I'll
# make it count the 'events' and sort them, but who knows. 

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#17.7.6_fixed_normal_....npy
#17.7.10/12_more_normal_....npy
# sort by epsilon value for normal trials
conc = (0.000541) #molar
all_eps = {}


all_eps['20']=['0.20_1_nc.npy','0.20_2_nc.npy','0.20_3_nc.npy','0.20_4_nc.npy','0.20_5_nc.npy','0.20_6_nc.npy','0.20_7_nc.npy','0.20_8_nc.npy','0.20_9_nc.npy','0.20_10_nc.npy','0.20_11_nc.npy','0.20_12_nc.npy','0.20_13_nc.npy','0.20_14_nc.npy','0.20_15_nc.npy']
all_eps['25']=['0.25_1_nc.npy','0.25_2_nc.npy','0.25_3_nc.npy','0.25_4_nc.npy','0.25_5_nc.npy','0.25_6_nc.npy','0.25_7_nc.npy','0.25_8_nc.npy','0.25_9_nc.npy','0.25_10_nc.npy','0.25_11_nc.npy','0.25_12_nc.npy','0.25_13_nc.npy','0.25_14_nc.npy','0.25_15_nc.npy']
all_eps['30']=['0.30_1_nc.npy','0.30_2_nc.npy','0.30_3_nc.npy','0.30_4_nc.npy','0.30_5_nc.npy','0.30_6_nc.npy','0.30_7_nc.npy','0.30_8_nc.npy','0.30_9_nc.npy','0.30_10_nc.npy','0.30_11_nc.npy','0.30_12_nc.npy','0.30_13_nc.npy','0.30_14_nc.npy','0.30_15_nc.npy']
all_eps['35']=['0.35_1_nc.npy','0.35_2_nc.npy','0.35_3_nc.npy','0.35_4_nc.npy','0.35_5_nc.npy','0.35_6_nc.npy','0.35_7_nc.npy','0.35_8_nc.npy', '0.35_10_nc.npy','0.35_11_nc.npy','0.35_12_nc.npy','0.35_13_nc.npy','0.35_14_nc.npy','0.35_15_nc.npy']
all_eps['40']=['0.40_1_nc.npy','0.40_2_nc.npy','0.40_3_nc.npy','0.40_4_nc.npy','0.40_5_nc.npy','0.40_6_nc.npy','0.40_7_nc.npy','0.40_8_nc.npy', '0.40_10_nc.npy','0.40_11_nc.npy','0.40_12_nc.npy','0.40_13_nc.npy','0.40_14_nc.npy','0.40_15_nc.npy']
all_eps['45']=['0.45_1_nc.npy','0.25_2_nc.npy','0.45_3_nc.npy','0.45_4_nc.npy','0.45_5_nc.npy','0.45_6_nc.npy','0.45_7_nc.npy','0.45_8_nc.npy', '0.45_10_nc.npy','0.45_11_nc.npy','0.45_12_nc.npy','0.45_13_nc.npy','0.45_14_nc.npy','0.45_15_nc.npy']
all_eps['50']=['0.50_1_nc.npy','0.50_2_nc.npy','0.50_3_nc.npy','0.50_4_nc.npy','0.50_5_nc.npy','0.50_6_nc.npy','0.50_7_nc.npy','0.50_8_nc.npy', '0.50_10_nc.npy','0.50_11_nc.npy','0.50_12_nc.npy','0.50_13_nc.npy','0.50_14_nc.npy','0.50_15_nc.npy']
all_eps['55']=['0.55_1_nc.npy','0.55_2_nc.npy','0.55_3_nc.npy','0.55_4_nc.npy','0.55_5_nc.npy','0.55_6_nc.npy','0.55_7_nc.npy','0.55_8_nc.npy', '0.55_10_nc.npy','0.55_11_nc.npy','0.55_12_nc.npy','0.55_13_nc.npy','0.55_14_nc.npy','0.55_15_nc.npy']
all_eps['60']=['0.60_1_nc.npy','0.60_2_nc.npy','0.60_3_nc.npy','0.60_4_nc.npy','0.60_5_nc.npy','0.60_6_nc.npy','0.60_7_nc.npy','0.60_8_nc.npy', '0.60_10_nc.npy','0.60_11_nc.npy','0.60_12_nc.npy','0.60_13_nc.npy','0.60_14_nc.npy','0.60_15_nc.npy']
all_eps['65']=['0.65_1_nc.npy','0.65_2_nc.npy','0.65_3_nc.npy','0.65_4_nc.npy','0.65_5_nc.npy','0.65_6_nc.npy','0.65_7_nc.npy','0.65_8_nc.npy', '0.65_10_nc.npy','0.65_11_nc.npy','0.65_12_nc.npy','0.65_13_nc.npy','0.65_14_nc.npy','0.65_15_nc.npy']
#removed all the 9's cause they didn't transfer
epsilons = []
kons = []
koffs = []
erron = []
erroff = []
def is_bound(nativecontacts):
    nc_cutoff = 50
    if nativecontacts > nc_cutoff:
        return True
    else:
        return False

def is_unbound(nativecontacts):
    if nativecontacts == 0:
        return True
    else:
        return False

for eps in all_eps.keys():
    btu_transition_count = 0
    utb_transition_count = 0
    # times in 10E-10 seconds
    total_time_bound = 0
    total_time_unbound = 0
    for isim, nativecontactspath in enumerate(all_eps[eps]):
        nativecontactsarray = np.load(nativecontactspath)
        timepoints = (len(nativecontactsarray[:,0]))
        # State reference info
        # state -1: unknown (we don't know if its bound or unbound)
        # state 0: bound (nc > 50)
        # state 1: unbound (nc = 0)
        statearray = np.ones(nativecontactsarray.shape[0])*-1
        for timepoint in range(timepoints):
            if is_bound(nativecontactsarray[timepoint,3]):
                statearray[timepoint] = 0
            elif is_unbound(nativecontactsarray[timepoint,3]):
                statearray[timepoint] = 1
            else: # we cannot currently say whether the timepoint fits the definition for the boudn or unbound state
                if timepoint > 0:
                    statearray[timepoint] = statearray[timepoint-1]
                else:
                    statearray[timepoint] = -1
            if timepoint > 0:
                if statearray[timepoint] == 1 and statearray[timepoint-1] == 0:
                    btu_transition_count += 1
                if statearray[timepoint] == 0 and statearray[timepoint-1] == 1:
                    utb_transition_count += 1

        # Now we have each timepoint labeled.  We now count how many timepoints are in each state!
        # The total number of timepoints in the bound state
        bound_tp_count = np.count_nonzero(statearray == 0)
        total_time_bound += bound_tp_count

        # The total number of timepoints in the unbound state
        unbound_tp_count = np.count_nonzero(statearray == 1)
        total_time_unbound += unbound_tp_count

    epsilons.append(int(eps))
    k_on = float((utb_transition_count)*1e10/(total_time_unbound*conc))
    se_k_on = ((utb_transition_count**0.5)*1e10/(total_time_unbound*conc))
#   #nt(eps)
#   #nt("k_on: {:.03e}".format(k_on))
#   #nt("k_off: {:.03e}".format(k_off))
#   # andard error of estimate of k_on

    k_off = float(btu_transition_count)*1e10/total_time_bound
    se_k_off = (btu_transition_count**0.5)*1e10/total_time_bound
#   nt("standard error of k_on: {:.03e}".format(se_k_on))
#   nt("standard error of k_off: {:.03e}".format(se_k_off))

    kons.append(k_on)
    koffs.append(k_off)
    erron.append(se_k_on)
    erroff.append(se_k_off)
errors = []
keqs = []
for i,thing in enumerate(kons):
    keq = thing/koffs[i]
    keqs.append(keq)
for i, thing in enumerate(erron):
    kon = kons[i]
    errkon = thing
    koff = koffs[i]
    errkoff = erroff[i]
    keq = keqs[i]
    if kon > 0 and koff > 0:
        error = keq*((errkon/kon)**2 + (errkoff/koff)**2)**0.5
        errors.append(error)
    else:
        error = 0
        errors.append(error)
exp_on = 1.9e6
exp_eq = 1.9e6/59
exp_off = 59
plt.figure(0)
plt.errorbar(epsilons,kons,yerr=erron,fmt='ro')
plt.yscale('log')
plt.axhspan(1.4e6,2.4e6,alpha=0.5,color='red')
plt.title('Estimates for kon, 3uS from bound')
plt.ylabel('k_on')
plt.xlabel('Epsilon')
plt.xticks(epsilons,epsilons)
plt.savefig('kons_3_b.png')
plt.figure(1)
plt.yscale('log')
keq_exp_error = exp_eq * np.sqrt((5e5/exp_on)**2 + (5/exp_off)**2)
experr_high = exp_eq - keq_exp_error
experr_low = exp_eq - keq_exp_error
plt.errorbar(epsilons,keqs,yerr=errors, fmt='ro')
plt.axhspan(experr_low,experr_high,alpha=0.5,color='red')
plt.title('Estimates for equilibrium constants,3 uS from bound, one std error')
plt.ylabel('Keq')
plt.xlabel('Epsilon')
plt.xticks(epsilons,epsilons)
plt.savefig('keqs_3_b.png')
plt.figure(2)
plt.yscale('log')
plt.errorbar(epsilons,koffs,yerr=erroff,fmt='ro')
plt.title('Estimates for koff, 3uS from bound')
plt.axhspan(55,63,alpha=0.5,color='red')
plt.ylabel('k_off')
plt.xlabel('Epsilon')
plt.xticks(epsilons,epsilons)
plt.savefig('koffs_3_b.png')

#plt.figure(1)
#plt.errorbar(epsilons,kons,yerr=erron,fmt='ro')
#plt.savefig('willthiswork.png')
#plt.figure(2)
#plt.errorbar(epsilons,koffs,yerr=erroff,fmt='ro')
#plt.savefig('willthiswork2.png')
