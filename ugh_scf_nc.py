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
    reference_map = np.zeros((246,246))
    for i_nc in xrange(atom_1.shape[0]):
        a1 = atom_1[i_nc]-1 
        c1 = chain_1[i_nc] 
        a2 = atom_2[i_nc]-1 
        c2 = chain_2[i_nc] 
        s = sigma[i_nc]
        idx1 = (a1 + 70*(c1-1))
        idx2 = (a2 + 70*(c2-1))
        reference_map[int(min((idx1), (idx2))), int(max((idx1), (idx2)))] = s
    return reference_map 

def count(compare_arr):
    '''
    Counts the number of native contacts within chain BC.
    '''
    return compare_arr[:,:].sum()

reference_map = load_reference_map('system.go.parameters')
u = md.Universe('MOVIE/movie.000000000.pdb', 'testout.xtc')
print(reference_map)
nc_arr = np.empty((len(u.trajectory)), dtype=float)
for i, ts in enumerate(u.trajectory):
    dist_arr = dist.distance_array(u.atoms.positions,u.atoms.positions)
    compare_arr = dist_arr < reference_map*1.2
    nc_arr[i]  = count(compare_arr)
print(nc_arr.shape)
np.save('nc_arr', nc_arr)



