import numpy as np 
import argparse, sys
from scipy.spatial.distance import cdist

# Parsing time
parser = argparse.ArgumentParser()
parser.add_argument("-ip", "--ipar", dest="param_file", default="internal.parameters",
                  help="internal parameter file")
parser.add_argument("-e", dest="eps", default=1.0, type=np.float64,
                  help="epsilon, everything is scaled off this one")
parser.add_argument("-a", dest="alpha", default=1.0, type=np.float64,
                  help="alpha, scales intramolecular epsilons with this")
args = parser.parse_args()

def param_loader(par_file):
    open_file = open(par_file)
    lines     = open_file.readlines()
    newlines = []
    for line in lines:
        split = line.split()
        #hard-coded; this just trying to see if it's a header
        if len(split[0]) == 1:
            if len(split) == 5:
                #must be bonds
                nline = map(int, split[0:3]) + map(np.float, split[3:])
#                nline = [int(split[0]), int(split[1]), int(split[2]),
#                         np.float(split[3]), np.float(split[4])]
            elif len(split) == 7:
                #must be angles
                nline = map(int, split[0:3]) + map(np.float, split[4:])
#                nline = [int(split[0]), int(split[1]), int(split[2]), int(split[3]),
#                         np.float(split[4]), np.float(split[5]), np.float(split[6])]
            elif len(split) == 15:
                #must be dihedrals
                nline = map(int, split[0:4]) + map(np.float, split[5:])
#                nline = [int(split[0]), int(split[1]), int(split[2]),
#                         int(split[3]), int(split[4]), 
#                         np.float(split[5]), np.float(split[6]),
#                         np.float(split[7]), np.float(split[8]), 
#                         np.float(split[9]), np.float(split[10]), 
#                         np.float(split[11]), np.float(split[12]), 
#                         np.float(split[13]), np.float(split[14])]
            else: 
                print("Too many fields")
                sys.exit()
        else:
            nline = line
            #this is if it's a header
        newlines.append(nline)
    return newlines

def mod_params(param, mps, pps):
    mk1,mk2,mv1,mv3 = mps
    pk1,pk2,pv1,pv3 = pps

    for line in param:
        try:
            ffield = int(line[0])
            llen = len(line)
            #lol doing the same check again cause he's a freaking dumbass
            if llen == 5:
                #again with the bond stuff
                fi = line[1]
                #fi is index for first atom
                si = line[2]
                #si is index for second atom
                if fi <= 85 and si <= 85:
                    #if it's within the first protein, this is intramolecular
                    line[3] = mk1
                elif fi > 85 and si > 85:
                    #^^intramolecular for the second protein`
                    line[3] = pk1
                elif (fi <= 85 and si > 85) or (fi > 85 and si <= 85):
                    #intermolecular -- shouldn't happen because two separate proteins  
                    print("Something is wrong with indices")
                    print("There shouldn't be bonds between peptides")
                    sys.exit()
            elif llen == 7:
                #okay this is angles
                fi = line[1]
                #first atom
                si = line[2]
                #second atom
                if fi <= 85 and si <= 85:
                    #angle for the first one
                    line[4] = mk2
                elif fi > 85 and si > 85:
                    #angle for the second one
                    line[4] = pk2
                elif (fi <= 85 and si > 85) or (fi > 85 and si <= 85):
                    #same thingy
                    print("Something is wrong with indices")
                    print("There shouldn't be angles between peptides")
                    sys.exit()
            elif llen == 15:
                #dihedral
                fi = line[1]
                si = line[2]
                if fi <= 85 and si <= 85:
                    #these are your dihedral barriers
                    line[5] = mv1
                    line[9] = mv1
                    line[6] = mv3
                    line[12] = mv3
                elif fi > 85 and si > 85:
                    #same thing for the second protein
                    line[5] = pv1
                    line[9] = pv1
                    line[6] = pv3
                    line[12] = pv3
                elif (fi <= 85 and si > 85) or (fi > 85 and si <= 85):
                    print("Something is wrong with indices")
                    print("There shouldn't be torsions between peptides")
                    sys.exit()
            else: 
                print(llen)
                print(line)
                print("Len too long in mod params")
                sys.exit()
        except:
            pass
    return param

def print_params(param):
    for line in param:
        try:
            ffield = int(line[0])
            llen = len(line)
            if llen == 5:
                print("{0[0]:>10}{0[1]:>10}{0[2]:>10}{0[3]:>10}{0[4]:>10}".format(line))
            elif llen == 7:
                print("{0[0]:>10}{0[1]:>10}{0[2]:>10}{0[3]:>10}{0[4]:>10}{0[5]:>10}{0[6]:>10}".format(line))
            elif llen == 15:
                print("{0[0]:>10}{0[1]:>10}{0[2]:>10}{0[3]:>10}{0[4]:>10}{0[5]:>10}{0[6]:>10}{0[7]:>10}{0[8]:>10}{0[9]:>10}{0[10]:>10}{0[11]:>10}{0[12]:>10}{0[13]:>10}{0[14]:>10}".format(line))
            else: 
                print("line has too many fields")
                sys.exit()
        except:
            print("".join(line[:-2]) + ")")

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
if __name__ == "__main__":
    # Pull MDM2 params
    eps = args.eps
    k1 = 100 * eps
    k2 = 20 * eps
    v1 = eps
    v3 = 0.5 * eps
    # Calc p53 params
    a = args.alpha
    pk1 = k1
    pk2 = a*k2
    pv1 = a*v1
    pv3 = a*v3
    # Load param file
    params = param_loader(args.param_file)
    params = mod_params(params, (k1,k2,v1,v3), (pk1,pk2,pv1,pv3))
    print_params(params)
