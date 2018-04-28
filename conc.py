import numpy as np

#dimer
moles = 1/6.02E23
cub_A = (4/3)*np.pi*(125)**3
print(cub_A)
Liters= cub_A * 1E-27
R = 1.987e-3
T = 298
conc = moles/Liters 
print('1 molecules',conc)

#trimer
moles = 3/6.02E23
cub_A = (4/3)*np.pi*(125)**3
Liters= cub_A * 1E-27
R = 1.987e-3
T = 298
conc_m = moles/Liters 
print('3 molecules',conc_m)

C_eff = 0.0078

C_eff_m = C_eff/3
print('trimer as 1/3 monomer',C_eff_m)
#as if we multiplied by 1L, divided by simulation volume
print('volumechange',C_eff_m/Liters)
C_m = C_eff*6.02E23*Liters
print(C_m)

print(0.0026*6.02E23)
