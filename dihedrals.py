#usr/bin/env python
#this is written by Ali Sinan Saglam ??, edited by Corinn 12/13


import numpy as np 
import argparse, sys
class Ass(object):
    def __init__(self):
        self.par_file = 'internal.parameters'
    def mod_params(self,par_file):
        new_lines = []
        open_file = open(self.par_file)
        lines     = open_file.readlines()
        for line in lines:
            llen = len(line.split())
            if llen == int(5):
                #print('bond')
                new_lines.append(line)
                pass
                #bonds  
            elif llen == int(7):
                #print('angle')
                new_lines.append(line)
                pass
                #angles
            elif llen == int(15):
                #print('dihedral')
                #dihedral
                c,i1,i2,i3,i4,v1_1,v3_1,d1,d2,v1_2,d3,d4,v3_2,d5,d6 = line.split()
                #if either is between 28 and 42 or 69 and 83
                loop1 = np.arange(70,93)
                loop2 = np.arange(158,183)
                loop1 = set(loop1)
                loop2 = set(loop2)
                i11,i21,i31,i41 = map(int,(i1,i2,i3,i4))
                indices = {i11,i21,i31,i41}
                if bool(indices.intersection(loop1.union(loop2))):
                    v1_1  = 0
                    v3_1  = 0
                    v1_2 = 0
                    v3_2 = 0
                new_line = [c,i1,i2,i3,i4,v1_1,v3_1,d1,d2,v1_2,d3,d4,v3_2,d5,d6]
                new_lines.append(new_line)
            else: 
                print("Len too long in mod params")
                sys.exit()
        return new_lines
    def format_line(self,line):
        c,i1,i2,i3,i4,v1_1,v3_1,d1,d2,v1_2,d3,d4,v3_2,d5,d6 = line
        l = "{:d}".format(int(c)).rjust(10) +\
            "{:d}".format(int(i1)).rjust(10) +\
            "{:d}".format(int(i2)).rjust(10) +\
            "{:d}".format(int(i3)).rjust(10) +\
            "{:d}".format(int(i4)).rjust(10) +\
            "{:07.5f}".format(float(v1_1)).rjust(10) +\
            "{:07.5f}".format(float(v3_1)).rjust(10) +\
            "{:.5f}".format(float(d1)).rjust(10) +\
            "{:.5f}".format(float(d2)).rjust(10) +\
            "{:07.5f}".format(float(v1_2)).rjust(10) +\
            "{:.5f}".format(float(d3)).rjust(10) +\
            "{:.5f}".format(float(d4)).rjust(10) +\
            "{:07.5f}".format(float(v3_2)).rjust(10) +\
            "{:.5f}".format(float(d5)).rjust(10) +\
            "{:.5f}".format(float(d6)).rjust(10) +\
            "\n"
        
        return l

    def print_params(self,param):
        for i,line in enumerate(param):
            if len(line) == 15:
                param[i] = self.format_line(line)
            else:
                pass
        
        print(''.join(param))

    
    #bonds        96 (expected number:       97)
    #1         1         2 100.00000   3.81263
    #molid iind jind eps dist
    #angles       94 (expected number:       94)
    #molid iid jid kid epsilon dist angl
    #1         1         2         3  20.00000   1.77772 101.85600
    #dihedrals    92 (expected number:       92)
    #molid i j k l v1 v3 dist dist v1 dist dist v3 dist dist
    #1        14        15        16        17   1.00000   0.50000   1.22823   3.68470   1.00000   0.33590   0.94190   0.50000  -0.85611  -0.51680
    
    # Calculate and print out
    def go(self):
      # self.print_params(
       params = self.mod_params(self.par_file) 
       self.print_params(params)
if __name__ == "__main__":
    A = Ass()
    A.go()
    # Pull MDM2 params
    # Load param file
#    params = param_loader(args.param_file)
#    params = mod_params(args.param_file)
