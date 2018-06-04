import numpy as np
import MDAnalysis as md
import os
import sys



todos = []
trajs = os.listdir('../../traj_segs/')
print(len(trajs))
for traj in trajs:
    ref = md.Universe('../../protein2.foldon.parm7')
    segs = os.listdir('../../traj_segs/'+str(traj))
    print(traj)
    traje = {seg: [] for seg in segs}
    for seg in segs:
        u = md.Universe('../../protein2.foldon.parm7','../../traj_segs/'+str(traj)+'/'+str(seg)+'/solute.nc')
        for ts in u.trajectory:
            traje[seg].append((u.atoms.radius_of_gyration()))
    todos.append(traje)
print(traje)

#todos becomes a list of dictionaries. the key of each dictionary is the segment name, and the values are the radius of gyration. hopefully.


np.save('rgyrs.npy',todos)




