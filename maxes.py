import MDAnalysis as md
import MDAnalysis.analysis.distances as dist
import numpy as np
import os

SCFs = os.listdir('../')
SCFs = filter(lambda x: x.startswith('SCF'),SCFs)
SCFs = filter(lambda x: x.endswith('.pdb'),SCFs)

todos = []
order = []
for pdb in SCFs:
    u = md.Universe('../'+pdb)
    distances = dist.self_distance_array(u.atoms.positions)
    todos.append(distances.max())
    order.append(pdb)

index = np.where(todos == max(todos))[0]
print('structure with max distance:',order[index])

np.save('max_distances',todos)
print('max distance',max(todos))
