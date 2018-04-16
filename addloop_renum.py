#!/usr/bin/env python
# takes output of renum_orient and adds loop, formats, renumbers (again)
#verified good shit 2.20.18
import numpy as np
import sys


def _format_line(atom, atomid, atomtype, residue, resnum, x,y,z):
    space = " "
    if len(atomtype) >= 4:
        l = "{:<4}".format(atom) + \
            "{:>7d}".format(atomid) + \
            "{:^6}".format(atomtype) + \
            "{:<4}".format(residue) + space + \
            "{:>4d}".format(resnum) + \
            "{:>12.3f}".format(x) + \
            "{:>8.3f}".format(y) + \
            "{:>8.3f}".format(z) + \
            "\n"
    else:
        l = "{:<4}".format(atom) + \
            "{:>7d}".format(atomid) + space*2 + \
            "{:<4}".format(atomtype) + \
            "{:<4}".format(residue) + space + \
            "{:>4d}".format(resnum) + \
            "{:>12.3f}".format(x) + \
            "{:>8.3f}".format(y) + \
            "{:>8.3f}".format(z) + \
            "\n"
    return l

arr = np.loadtxt('1rfo_oriented_clean.pdb',dtype='string')
loop = np.loadtxt('loop.pdb',dtype='string')
lines = []
for i in range(len(arr[:432,0])):
    atom     = arr[i,0]
    atomid   = int(arr[i,1])
    atomtype = arr[i,2]
    residue  = arr[i,3]
    resnum   = int(arr[i,4])
    x        = float(arr[i,5])
    y        = float(arr[i,6])
    z        = float(arr[i,7])
    lines.append(_format_line(atom,atomid,atomtype,residue,resnum,x,y,z))
#good
resnum = 28
for i in range(len(loop)):
    atom     = loop[i,0]
    atomid   = 433+i
    atomtype = loop[i,2]
    residue  = loop[i,3]
    if loop[i-1,3] == loop[i,3]:
        pass
    else:
        resnum += 1
    x        = float(loop[i,5])
    y        = float(loop[i,6])
    z        = float(loop[i,7])
    lines.append(_format_line(atom,atomid,atomtype,residue,resnum,x,y,z))

#good
resnum = 42
for i in range(len(arr[432:864,0])):
    atom     = arr[i+432,0]
    atomid   = int(arr[i+432,1])+137
    atomtype = arr[i+432,2]
    residue  = arr[i+432,3]
    if arr[(i+432)-1,4] == arr[(i+432),4]:
        pass
    else:
        resnum += 1
    x        = float(arr[i+432,5])
    y        = float(arr[i+432,6])
    z        = float(arr[i+432,7])
    lines.append(_format_line(atom,atomid,atomtype,residue,resnum,x,y,z))
#good
resnum = 70
for i in range(len(loop)):
    atom     = loop[i,0]
    atomid   = 1002+i
    atomtype = loop[i,2]
    residue  = loop[i,3]
    if loop[i-1,3] == loop[i,3]:
        pass
    else:
        resnum += 1
    x        = float(loop[i,5])
    y        = float(loop[i,6])
    z        = float(loop[i,7])
    lines.append(_format_line(atom,atomid,atomtype,residue,resnum,x,y,z))
#good
resnum = 84
for i in range(len(arr[864:,0])):
    atom     = arr[i+864,0]
    atomid   = int(arr[i+864,1])+((137*2))
    atomtype = arr[i+864,2]
    residue  = arr[i+864,3]
    if arr[(i+864)-1,4] == arr[(i+864),4]:
        pass
    else:
        resnum += 1
    x        = float(arr[i+864,5])
    y        = float(arr[i+864,6])
    z        = float(arr[i+864,7])
    lines.append(_format_line(atom,atomid,atomtype,residue,resnum,x,y,z))


with open('1rfo_looped.pdb', 'w+') as outfile:
    for line in lines:
        outfile.write(line)
