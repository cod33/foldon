#!/usr/bin/env python
#
# random_orientation.py
#
# Randomly orient 3 foldon trimers such that the centers of mass for each
# subunit forms an equilateral triangle with a specified edge length. 
#
# Written 17.4.27 by Alex DeGrave
#
# Adapted sometime in late 2017 by Corinn Durham to exclude 60 degree cone around axis of monomer AB
#
from __future__ import print_function
import argparse
import MDAnalysis
import MDAnalysis.analysis.distances
import numpy as np
import numpy
import sys

class RandomOrientation(object):
    def __init__(self):
        self.parse_args()
        self.edgelength = self.args.edgelength 
        self.load_chains()
        self.threshold=0.01

    def parse_args(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('--edgelength', dest='edgelength', default=20,
                                 type=float,
                                 help="Separate foldon subunits by the " 
                                      "specified length, in Angstroms" 
                                 )

        self.parser.add_argument('--inputpath', dest='inputpath', default='1rfo.pdb',
                                 help="Path to the PDB file for the foldon "
                                      "trimer (should be something like "
                                      "1rfo.pdb"
                                 )

        self.parser.add_argument('--outputpath', dest='outputpath', 
                                 default='1rfo_oriented.pdb',
                                 help="Path to which a new PDB file will be "
                                      "output, containing the 3 monomers at "
                                      "the requested separation and rotated "
                                      "randomly."
                                 )
     
        self.parser.add_argument('--surface', dest='surface', 
                                 default='surface',
                                 action='store_true',
                                 help="If this argument is supplied, calculate "
                                      "positions of the monomers such that the "
                                      "nearest atoms in each monomer are "
                                      "'edgelength' Angstroms apart.  If this "
                                      "argument is not supplied, base the "
                                      "calculation on the center of mass."
                                 )

        

        self.args = self.parser.parse_args()

    def random_rotation_matrix(self):
        '''
        Generate a random (uniform) 3 by 3 rotation matrix.  
        '''
        # Pick a random point on the unit sphere (S3)
        # See http://mathworld.wolfram.com/SpherePointPicking.html
        theta = numpy.random.random()*2*numpy.pi
        phi = numpy.arccos(2*numpy.random.random()-1)
        x = numpy.cos(theta)*numpy.sin(phi)
        y = numpy.sin(theta)*numpy.sin(phi)
        z = numpy.cos(phi)

        # The uniform random rotation matrix is then given by the rotation
        # matrix U such that Uv1=v2, as shown below.
        v1 = numpy.array((1,0,0))
        v2 = numpy.array((x,y,z))

        # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        v_cross = numpy.cross(v1, v2)        
        m = numpy.array(((0            ,-1*v_cross[2],    v_cross[1]),
                         (   v_cross[2], 0           , -1*v_cross[0]),
                         (-1*v_cross[1],   v_cross[0], 0            )))

        U = numpy.identity(3) + m + m.dot(m)*(1/(1+v1.dot(v2)))

        return U

    def load_chains(self):
        self.u = MDAnalysis.Universe(self.args.inputpath)
        self.chainA = self.u.select_atoms('bynum 1-432')
        self.chainB = self.u.select_atoms('bynum 433-864')
        self.chainC = self.u.select_atoms('bynum 865-1296')

    def translate_equilateral(self, edgelength):

        pos1 = numpy.array((0, 
                            -1*self.edgelength,
                            0)
                            #((4*np.pi/3. - np.pi/3.)* np.random.random() + np.pi/3.)*self.edgelength)
                          )

        pos2 = numpy.array((0,0,0)) 
 

#        x = ((1+1)*np.random.random() -(1*1))*self.edgelength
        x = np.random.random()*self.edgelength
        y = ((1+1)*np.random.random() -(1*1))*self.edgelength
        z = ((1+1)*np.random.random() -(1*1))*self.edgelength
        
        print(x,y,z)

        vect1 = self.edgelength/np.sqrt(float(x)**2+float(y)**2+float(z)**2)
        print('vect1',vect1)
        vect2 = [x,y,z]
        print('vect2',vect2)
        coords = np.array(np.multiply(vect1,vect2))
        print('coords',coords)

        self._check_it(coords)
        
        pos3 = numpy.array((coords))

        self.chainA.translate(pos1)
        self.chainB.translate(pos2)
        self.chainC.translate(pos3)

    def _check_it(self, coords):
        '''
        excluding the cone around A
        '''
        
        if (-1*self.edgelength*np.sqrt(3)/2. <= coords[0] <= self.edgelength*np.sqrt(3)/2.) and (-1.0*self.edgelength <= coords[1] <= -0.5*self.edgelength):
            sys.exit('nope')
        else:
            return(coords)
#        pos1 = numpy.array((self.edgelength/numpy.sqrt(3),
#                            0,
#                            0)
#                           )
#
#        pos2 = numpy.array((-1./2*self.edgelength/numpy.sqrt(3),
#                            self.edgelength/2,
#                            0)
#                           )
# 
#        pos3 = numpy.array((-1./2*self.edgelength/numpy.sqrt(3),
#                            -1*self.edgelength/2,
#                            0)
#                           )
#
#        self.chainA.translate(pos1)
#        self.chainB.translate(pos2)
#        self.chainC.translate(pos3)
#
    def separate_surfaces(self):
        distarray = MDAnalysis.analysis.distances.distance_array
        # Calculate the protein diameter
        diams = []
        for chain in [self.chainA, self.chainB, self.chainC]:
            dists = distarray(chain.positions,
                              chain.positions)
            diams.append(dists.max())
        # Take a guess at the correct separation
        
        self.translate_equilateral(max(diams)+self.edgelength)

        error = self.threshold*10 
        # Iteratively move the subunits to achieve the desired separation
        # between the surfaces.
        div = 2
        while error > self.threshold:

            # Calculate the vectors between the centers of geometry, so
            # that we know which way to move.
            p1 = self.chainA.center_of_geometry()
            p2 = self.chainB.center_of_geometry()
            v12 = p1-p2
            v12 /= numpy.linalg.norm(v12)
       
            # Start by moving chain A          
            distAB = distarray(self.chainA.positions,
                               self.chainB.positions).min()
            self.chainA.translate(v12*(self.edgelength-distAB)/div)


            # Calculate the vectors between the centers of geometry, so
            # that we know which way to move.
            p2 = self.chainB.center_of_geometry()
            p3 = self.chainC.center_of_geometry()
            v23 = p3-p2
            v23 /= numpy.linalg.norm(v23)
       
            # Move chain B 
            distCB = distarray(self.chainC.positions,
                               self.chainB.positions).min()
            self.chainC.translate(v23*(self.edgelength-distCB)/div)


            dist_AC = distarray(self.chainA.positions,self.chainC.positions).min()
            if dist_AC < self.edgelength:# Calculate the vectors between the centers of geometry, so
            # that we know which way to move.
                p1 = self.chainA.center_of_geometry()
                p3 = self.chainC.center_of_geometry()
                v31 = p3-p1
                v31 /= numpy.linalg.norm(v31)
           
                # Move chain C 
                distAC = distarray(self.chainA.positions,
                                   self.chainC.positions).min()
                self.chainC.translate(v31*(self.edgelength-distAC)/div)

            error  = abs(distarray(self.chainA.positions,
                                   self.chainB.positions).min() \
                         - self.edgelength)
            error += abs(distarray(self.chainB.positions,
                                   self.chainC.positions).min() \
                         - self.edgelength)
            if dist_AC < self.edgelength:
                error += abs(distarray(self.chainA.positions,
                                   self.chainC.positions).min() \
                         - self.edgelength)

            sys.stdout.write("\rSummed deviation from ideal separation: "
                             "{:.04f}".format(error))
            sys.stdout.flush()

    def go(self) :

        #### Testing ####
        for chain in [self.chainA, self.chainB, self.chainC]:
            chain.atoms.translate(-1*chain.center_of_geometry())
        oldApos = self.chainA.positions
        oldBpos = self.chainB.positions
        oldCpos = self.chainC.positions

        for chain in [self.chainA, self.chainB, self.chainC]:
            rand_rot_mat = self.random_rotation_matrix()
            chain.rotate(rand_rot_mat)

        newApos = self.chainA.positions
        newBpos = self.chainB.positions
        newCpos = self.chainC.positions

        #for chain in [self.chainA, self.chainB, self.chainC]:
        #    chain.atoms.translate(-1*chain.center_of_geometry())
        #    rand_rot_mat = self.random_rotation_matrix()
        #    chain.rotate(rand_rot_mat)

        if not self.args.surface:
            self.translate_equilateral(self.edgelength)
        else:
            self.separate_surfaces()

       
        self.u.atoms.write(self.args.outputpath)

if __name__ == "__main__":
    r = RandomOrientation() 
    r.go()
