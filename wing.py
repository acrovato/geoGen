#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package FPE grid creator
#
# Create a rectangular unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow Full Potential solver
# Adrien Crovato

import numpy as np

## Handle wing data
#
# Adrien Crovato
class Wing:
    def __init__(self, filenames, span, taper, sweep, dihedral, twist, rootChord):
        # Number of airfoils
        self.n = len(filenames)

        # Convert degrees to radians
        sweep = [x * np.pi/180 for x in sweep]
        dihedral = [x * np.pi/180 for x in dihedral]
        twist = [x * np.pi/180 for x in twist]

        # Compute wing shape parameters
        self.shape(span, taper, rootChord)

        # Create airfoil points and indices
        self.airfoil(filenames, span, twist, sweep, dihedral)
        #tip()

        # find special numbers

        # Line numbering
        #self.linN = [np.zeros([6,1])]*(n-1)
        # Surface numering
        #self.surN = [np.zeros([6,1])]*(n-1)

    def shape(self, span, taper, rootChord):
        """Compute basic shape parameters of the wing
        """
        areas = []
        self.chord = []
        self.spanPos = []
        self.chord.append(rootChord)
        self.spanPos.append(0.)
        for i in range(1, self.n):
            self.chord.append(self.chord[i-1]*taper[i-1])
            areas.append((self.chord[i-1]+self.chord[i])*span[i-1]/2)
            self.spanPos.append(self.spanPos[i-1]+span[i-1])
        self.S = sum(areas)
        self.b = sum(span)

    def airfoil(self, filenames, span, twist, sweep, dihedral):
        """Read, transform and store airfoil points
        """
        self.pts = []
        self.ptsN = []
        # read and store coordinates
        for i in range(0, self.n):
            aPts = self.read(filenames[i])
            size = aPts.shape[0]
            aPts = np.hstack((aPts, self.spanPos[i]*np.ones([size,1])))
            aPts[:,[1,2]] = np.fliplr(aPts[:,[1,2]])
            aIdx = np.transpose(np.arange(i*250+1, i*250+1+size))
            self.pts.append(aPts)
            self.ptsN.append(aIdx)
        # transform coordinates
        for i in range(0, self.n):
            # apply taper (scaling)
            self.pts[i][:, [0,2]] = self.pts[i][:, [0,2]]*self.chord[i]
            # aplly twist (rotation)
            self.pts[i][:, [0,2]] = np.dot(self.pts[i][:, [0,2]],np.array([[np.cos(twist[i]), -np.sin(twist[i])],[np.sin(twist[i]), np.cos(twist[i])]]))
        for i in range(1, self.n):
            # apply sweep (translation)
            self.pts[i][:, 0] = self.pts[i][:, 0] + np.min(self.pts[i-1][:, 0]) + np.tan(sweep[i-1])*span[i-1]
            # apply dihedral (translatation)
            self.pts[i][:, 2] = self.pts[i][:, 2] + sum(np.tan(dihedral[0:i-1])*span[0:i-1])

    def tip(self):
        """Create tip points
        """
        self.tip = []
        self.tipN = []

    def read(self,fname):
        """Read data from file and stroe in matrix
        """
        _file = file(fname)
        label = _file.next().split(',')
        _file.close()
        data = np.loadtxt(fname, skiprows=1)
        return data

    def writeInfo(self,fname):
        """Write wing geometrical parameters
        """
        file = open(fname, 'a')
        file.write('// --- Wing geometry ---\n')
        file.write('// Number of spanwise stations: {0:d}\n'.format(self.n))
        file.write('// Spanwise stations normalized coordinate: ')
        for p in self.spanPos:
            file.write('{0:f} '.format(p/self.b))
        file.write('\n')
        file.write('// Chord lengths: ')
        for c in self.chord:
            file.write('{0:f} '.format(c))
        file.write('\n')
        file.write('// Half-wing area: {0:f}\n'.format(self.S))
        file.write('// Half-wing span: {0:f}\n'.format(self.b))
        file.write('\n')
        file.close()

    def writeOpts(self, fname):
        """Write wing gmsh options
        """
        file = open(fname, 'a')
        file.write('// --- Wing options ---\n')
        for i in range(0, self.n):
            file.write('DefineConstant[ msLe{0:1d} = {{ {1:f}, Name "leading edge mesh size on {2:1d}th spanwise station" }} ];\n'.format(i, self.chord[i]/100, i))
            file.write('DefineConstant[ msTe{0:1d} = {{ {1:f}, Name "trailing edge mesh size on {2:1d}th spanwise station" }} ];\n'.format(i, self.chord[i]/100, i))
            file.write('DefineConstant[ gr{0:1d} = {{ {1:f}, Name "growth ratio for {2:1d}th spanwise station" }} ];\n'.format(i, 1.5, i))
        file.write('\n')
        file.close()

    def writePoints(self,fname):
        """Write wing points
        """
        file = open(fname, 'a')
        file.write('// --- Wing points ---\n')
        for i in range(0, self.n):
            file.write('// -- Airfoil {0:d}\n'.format(i))
            for j in range(0, self.pts[i].shape[0]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            file.write('\n')
        file.write('\n')
        file.close()

    def writeLines(self):
            """Write wing lines
            """

    def writeSurface(self):
            """Write wing line loops and surfaces
            """

    def writePhysical(self):
            """Write wing physical groups
            """