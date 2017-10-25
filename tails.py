import numpy as np
import MDAnalysis as md
import MDAnalysis.analysis.rms as rms
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

u = md.Universe('MOVIE/movie.000000000.pdb','testout.xtc')
ref = md.Universe('/home/cod33/foldon/templates/references/1rfo_beads.pdb')

#A = u.A.atoms
#B = u.B.atoms
#C = u.C.atoms
#A_nterm = A[0:10]
#A_cterm = A[-10:]
#B_nterm = B[0:10]
#B_cterm = B[-10:]
#C_nterm = C[0:10]
#C_cterm = C[-10:]

A_nterm = rms.RMSD(u,ref,select='bynum 0-5')
A_cterm = rms.RMSD(u,ref,select='bynum 65-70')
B_nterm = rms.RMSD(u,ref,select='bynum 70-75')
B_cterm = rms.RMSD(u,ref,select='bynum 135-140')
C_nterm = rms.RMSD(u,ref,select='bynum 140-145')
C_cterm = rms.RMSD(u,ref,select='bynum 205-210')

for obj in [A_nterm,A_cterm,B_nterm,B_cterm,C_nterm,C_cterm]:
    obj.run()
full_arr=[A_nterm.rmsd[:,0],A_nterm.rmsd[:,2],A_cterm.rmsd[:,2],B_nterm.rmsd[:,2],B_cterm.rmsd[:,2],C_nterm.rmsd[:,2],C_cterm.rmsd[:,2]]

np.save('terms.npy',full_arr)
