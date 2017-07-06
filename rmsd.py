#!/usr/bin/env python2
import MDAnalysis as MDA
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis.analysis.rms as RMS
import argparse

class RMSD(object):
    def __init__(self):
       self._parseargs()

    def _parseargs(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--coordinate-file", dest='coordinate_path', 
            default = "testout.xtc",
            #default = "temp.xtc",
            help = "potato, please input coordinate file to be calculated from. ")
        parser.add_argument("--reference-file", dest='reference_path',
            default = "/home/cod33/foldon/templates/reference/trimer_reference.pdb",
            help = "potato, please input reference file")
        self.args = parser.parse_args()

    def run(self):
        u_coord = MDA.Universe(self.args.reference_path, self.args.coordinate_path)
        u_ref = MDA.Universe(self.args.reference_path)
        R_tot = RMS.RMSD(u_coord, u_ref)
        R_A = RMS.RMSD(u_coord,u_ref,select = "bynum 1-70")
        R_B = RMS.RMSD(u_coord,u_ref,select = "bynum 71-140")
        R_C = RMS.RMSD(u_coord,u_ref,select = "bynum 141-210")
        for R in [R_tot, R_A, R_B, R_C]:
            R.run()
        tot_arr = R_tot.rmsd[:,2] 
        a_arr = R_A.rmsd[:,2] 
        b_arr = R_B.rmsd[:,2] 
        c_arr = R_C.rmsd[:,2] 
        time = R_tot.rmsd[:,1]
        full_arr = numpy.vstack((time, tot_arr, a_arr, b_arr, c_arr))
        numpy.save('output', full_arr)

if __name__ == "__main__":
    rmsd = RMSD()    
    rmsd.run()



arr = numpy.load('output.npy')
for inx in range(1,5):
    plt.plot(arr[0],arr[inx], label="{:d}".format(inx))
Total, = plt.plot(arr[0],arr[1], 'r', label="Total")
A, = plt.plot(arr[0],arr[2], 'y', label="Chain A")
B, = plt.plot(arr[0],arr[3], 'g', label="Chain B")
C, = plt.plot(arr[0],arr[4], 'b', label="Chain C")
plt.legend(handles = [Total, A, B, C])
plt.ylabel('RMSD')
plt.xlabel('time')
#plt.show()
plt.savefig('rmsd_plot.png')
