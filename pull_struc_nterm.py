#exes,whys,etc from plot_hist.py

import numpy as np

exes = np.load('../exes_nterm_wt.npy')
whys = np.load('../whys_nterm_wt.npy')
idxes = np.load('../idxes_nterm_wt.npy')
tryps = np.load('../tryps_nterm_wt.npy')
evens = list(range(0,168,2))

max_y = np.where(whys == max(whys))
print('max y =', max_y)

ex = exes[max_y]
#prints 2,2 array
print('x at max y =', ex)
exmin = ex[0,1]
exmax = ex[1,1]

goods = []

for idx in evens:
    goods.append((np.where(np.logical_and(tryps[idx]<=exmax,tryps[idx]>=exmin)),idx))

goods = np.array(goods)
#goods = array (84,2) where goods[x,0][0]=list of ts
sims = range(83)
#print(goods[1])
#print(goods[1,0][0])
#print(goods[1,1])
pullable = []
for sim in sims:
    tss = goods[sim,0][0]
    idx = goods[sim,1]
    if len(tss)>1:
        for ts in tss:
            if ts < 30000:
                if str(ts)[-1]==str(0):
                    pullable.append((ts,idxes[int(idx/2)][1][3:6]))
                else:
                    pass
            else:
                pass
        
    else:
        pass

pullable = np.array(pullable)

np.save('pullable_1st_wt_nterm',pullable)
#takes xs, cuts before global max
xs_2 = np.where(exes[:,1]<3.5)[0]

ys_2 = []
for x in xs_2:
    ys_2.append(whys[x])

ys_2_max = np.where(ys_2 == max(ys_2))
print('2nd max y =', max(ys_2))

x2 = exes[ys_2_max]
print('x at 2nd max', x2)
x2min = x2[0,1]
x2max = x2[1,1]

goods_2 = []
for idx in evens:
    goods_2.append((np.where(np.logical_and(tryps[idx]<=x2max,tryps[idx]>=x2min)),idx))

goods2 = np.array(goods_2)
pullable2 = []

for sim in sims:
    tss = goods2[sim,0][0]
    idx = goods2[sim,1]
    if len(tss)>1:
        for ts in tss:
            if ts < 30000:
                if str(ts)[-1]==str(0):
                    pullable2.append((ts,idxes[int(idx/2)][1][3:6]))
                else:
                    pass
            else:
                pass
        
    else:
        pass

pullable2 = np.array(pullable2)

np.save('pullable_2nd_wt_nterm',pullable2)

xs_3 = np.where(exes[:,1]<3.2)[0]

ys_3 = []
for x in xs_3:
    ys_3.append(whys[x])

ys_3_max = np.where(ys_3 == max(ys_3))
print('3rd max y =', max(ys_3))

x3 = exes[ys_3_max]
print('x at 3rd max', x3)
x3min = x3[0,1]
x3max = x3[1,1]

goods_3 = []


for idx in evens:
    goods_3.append((np.where(np.logical_and(tryps[idx]<=x3max,tryps[idx]>=x3min)),idx))

goods3 = np.array(goods_3)
pullable3 = []
for sim in sims:
    tss = goods3[sim,0][0]
    idx = goods3[sim,1]
    if len(tss)>1:
        for ts in tss:
            if ts < 30000:
                if str(ts)[-1]==str(0):
                    pullable3.append((ts,idxes[int(idx/2)][1][3:6]))
                else:
                    pass
            else:
                pass
        
    else:
        pass

pullable3 = np.array(pullable3)

np.save('pullable_3rd_wt_nterm',pullable3)
