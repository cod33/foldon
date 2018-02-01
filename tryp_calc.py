
import MDAnalysis as md
import MDAnalysis.analysis.distances as dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#u = md.Universe('%TYPE%_%STATE%_%ORIENT%.pdb','%TYPE%_%STATE%_%ORIENT%.xtc')
u = md.Universe('movie.000000000.pdb','testout.xtc')

AB = []
AC = []
BC = []
for ts in u.trajectory:
    dist_arr = dist.distance_array(u.atoms.positions,u.atoms.positions)
    A_B = dist_arr[48:52,118:122].min()
    A_C = dist_arr[48:52,188:192].min()
    B_C = dist_arr[118:122,188:192].min()
    AB.append(A_B)
    AC.append(A_C)
    BC.append(B_C)
AB_arr = np.asarray(AB)
AC_arr = np.asarray(AC)
BC_arr = np.asarray(BC)
all_arr = np.vstack((AB_arr,AC_arr,BC_arr))
np.save('tryp_arr',all_arr)
#np.save('%TYPE%_%STATE%_%ORIENT%',all_arr)

#SCF
u_SCF = md.Universe('%TYPE%_%STATE%_%ORIENT%.pdb','%TYPE%_%STATE%_%ORIENT%.xtc')

AB_SCF = []
AC_SCF = []
BC_SCF = []
for ts in u_SCF.trajectory:
    dist_arr = dist.distance_array(u.atoms.positions,u.atoms.positions)
    A_B = dist_arr[48:52,139:143].min()
    A_C = dist_arr[48:52,230:233].min()
    B_C = dist_arr[139:143,230:233].min()
    AB_SCF.append(A_B)
    AC_SCF.append(A_C)
    BC_SCF.append(B_C)
AB_arr = np.asarray(AB)
AC_arr = np.asarray(AC)
BC_arr = np.asarray(BC)
all_arr = np.vstack((AB_arr,AC_arr,BC_arr))
np.save('%TYPE%_%STATE%_%ORIENT%',all_arr)
