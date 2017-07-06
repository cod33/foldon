#!/usr/bin/sh python
#Corinn 17.5.9
import MDAnalysis as md
import MDAnalysis.analysis.distances as dist 
import numpy as np


def load_reference_map(filepath):
    '''
    reference map for native contacts in chain A
    <-chainA-><--chainB-->
    -------------------------------
    |x 4 6 7 8 0 0 0 
    |4 x 9 8 5
    |6 9 x 6 5
    |7 8 6 x 6
    |8 5 5 6 x
    |0
    |0
    |0
    |
    |
    |
    |
    '''
    native = np.loadtxt(filepath, skiprows = 1)
    chain_1 = native[:,0]
    atom_1 = native[:,1]
    chain_2 = native[:,2]
    atom_2 = native[:,3]
    sigma = native[:,4]
    reference_map = np.zeros((210,210))
    for i_nc in xrange(atom_1.shape[0]):
        a1 = atom_1[i_nc] 
        c1 = chain_1[i_nc] 
        a2 = atom_2[i_nc] 
        c2 = chain_2[i_nc] 
        s = sigma[i_nc]
        idx1 = a1 + 70*(c1-1)
        idx2 = a2 + 70*(c2-1)
        reference_map[min(idx1, idx2), max(idx1, idx2)] = s
    #print(np.count_nonzero(reference_map[:70,:70]))
    return reference_map

def count_nc_chain_a(compare_arr):
    '''
    Counts the number of native contacts within chain A.
    '''
    return compare_arr[:70,:70].sum()

def count_nc_chain_b(compare_arr):
    '''
    Counts the number of native contacts within chain B.
    '''
    return compare_arr[70:140,70:140].sum()

def count_nc_chain_c(compare_arr):
    '''
    Counts the number of native contacts within chain C.
    '''
    return compare_arr[140:210,140:210].sum()

def count_nc_chain_ab(compare_arr):
    '''
    Counts the number of native contacts within chain AB.
    '''
    return compare_arr[:70,70:140].sum()

def count_nc_chain_ac(compare_arr):
    '''
    Counts the number of native contacts within chain AC.
    '''
    return compare_arr[:70,140:210].sum()

def count_nc_chain_bc(compare_arr):
    '''
    Counts the number of native contacts within chain BC.
    '''
    return compare_arr[70:140,140:210].sum()

reference_map = load_reference_map('system.go.parameters')
u = md.Universe('MOVIE/movie.000000000.pdb', 'testout.xtc')
nc_arr = np.empty((len(u.trajectory),6), dtype=float)
for i, ts in enumerate(u.trajectory):
    dist_arr = dist.distance_array(u.atoms.positions,u.atoms.positions)
    compare_arr = dist_arr < reference_map*2
    nc_arr[i,0]  = count_nc_chain_a(compare_arr)
    nc_arr[i,1]  = count_nc_chain_b(compare_arr)
    nc_arr[i,2]  = count_nc_chain_c(compare_arr)
    nc_arr[i,3]  = count_nc_chain_ab(compare_arr)
    nc_arr[i,4]  = count_nc_chain_ac(compare_arr)
    nc_arr[i,5]  = count_nc_chain_bc(compare_arr)
np.save('nc_arr', nc_arr)

max_nc_a = 366/2
max_nc_b = 354/2
max_nc_c = 354/2
max_nc_ab = 342/2
max_nc_ac = 318/2
max_nc_bc = 324/2


#<chainID 1> <atomID 1> <chainID 2> <atomID 2> <sigma> <epsilon>
#native_A = ()
#native_B = ()
#native_C = ()
#native_AB = ()
#native_BC = ()
#native_AC = ()
#
#
#
#for i in range(0, len(chain_1)+1):
#    if chain_1[i] == chain_2[i] == 1:
#        native_A += 1
#    if chain_1[i] == chain_2[i] == 2:
#        native_B += 1
#    if chain_1[i] == chain_2[i] == 3:
#        native_c += 1
#    if chain_1[i] == 1 and chain_2[i] == 2:
#        native_AB += 1
#    if chain_1[i] ==1 and chain_2[i] == 3:
#        native_AC += 1
#    if chain_1[i] == 2 and chain_2[i] ==3:
#        native_BC += 1
#
#
#A = u.select_atoms('segid A')
#B = u.select_atoms('segid B')
#C = u.select_atoms('segid C')
#
#pos_A = A.positions
#pos_B = B.positions
#pos_C = C.positions
#
#for ts in u.trajectory:
#    for i in range(0, 69):
#        dif = (pos_A[i] - pos_A[i+1]) * 1.2
#        dif_A.append(dif)
#
#    
##YOU CAN CHOOSE ATOM_1 = NATIVE_A[:,0] ETC. USE FULL NATIVE FILE PROBABLY. 
