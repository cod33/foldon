import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

moles = 2/6.02E23
cub_A = (4/3)*np.pi*(125)**3
Liters= cub_A * 1E-27
R = 1.987e-3
T = 298
conc = moles/Liters 
exp_off = 59
#per second
exp_on = 1.9e6*conc
exp_eq = exp_on/exp_off
keq_exp_error = exp_eq * np.sqrt((5e5*conc/exp_on)**2 + (5/exp_off)**2)
#this is unitless 
experr_high = (exp_eq + keq_exp_error)
experr_low = (exp_eq - keq_exp_error)
g_error = (-1)*R*T*((keq_exp_error/exp_eq))
g_high = (-1)*R*T*np.log(exp_eq)+ g_error
g_low = (-1)*R*T*np.log(exp_eq) - g_error
g_exp_estimate = (-1)*R*T*np.log(exp_eq)

k_arr = np.load('k_arr_vstime.npy')


step = 1000
x_range = range(len(k_arr[0]))
print(x_range)
higherr_eq = [k_arr[4,x]+k_arr[5,x] for x in x_range]
lowerr_eq = [k_arr[4,x]-k_arr[5,x] for x in x_range]
xs = [(x*1e-4) for x in x_range]
plt.scatter(xs,k_arr[0],s=2,c='b')
#plt.errorbar(xs,bound[0],yerr =bound[1],ecolor='blue',elinewidth=2,ma=0.25,marker='o',fmt='none')
plt.fill_between(xs,higherr_eq,lowerr_eq,facecolor='blue',alpha=0.25)
plt.axhspan(g_low,g_high,alpha=0.25,color='black')
plt.axhline(y=g_exp_estimate,color='black')
#plt.plot(xs[step+1:],[(bound[0,i:i+step].mean()) for i in range(int(len(bound[0]))-int(step+1))],'k')
#step = 1000
##unbound doesn't have events til ~5000fs, is ok
#plt.scatter(xs,unbound[0],s=2,c='r')
#higherr_unbound = [unbound[0,x]+unbound[1,x] for x in x_range]
#lowerr_unbound = [unbound[0,x]-unbound[1,x] for x in x_range]
#xs = [(x*1e-4) for x in x_range]
##plt.errorbar(xs,unbound[0],yerr =unbound[1],ecolor='red',elinewidth=2,ma=0.25,marker='o',fmt='none')
#plt.fill_between(xs,higherr_unbound,lowerr_unbound,facecolor='red',alpha=0.25)
#plt.plot(xs[step+1:],[(unbound[0,i:i+step].mean()) for i in range(int(len(unbound[0]))-int(step+1))],'k')
plt.xlabel(r'time ($\mu s)$')
plt.ylabel(r'$\Delta G$')
plt.savefig('pls_god.png')
