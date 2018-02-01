import numpy as np

#dimer
moles = 2/6.02E23
cub_A = (4/3)*np.pi*(125)**3
Liters= cub_A * 1E-27
R = 1.987e-3
T = 298
conc = moles/Liters 
print('dimer',conc)

#trimer
moles = 3/6.02E23
cub_A = (4/3)*np.pi*(125)**3
Liters= cub_A * 1E-27
R = 1.987e-3
T = 298
conc = moles/Liters 
print('trimer',conc)
