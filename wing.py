#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' 
Copyright 2020 University of Liege

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
'''

## @package GeoGen (CFD basic grid creator)
#
# Create an unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow or SU2 CFD solvers
# Adrien Crovato

import numpy as np

## Handle wing data
#
# Adrien Crovato
class Wing:
    def __init__(self, filenames, span, taper, sweep, dihedral, twist, rootChord, offset):
        # Number of airfoils
        self.n = len(filenames)
        if self.n > 10:
            raise Exception('Wing: read {0:d} airfoils but a maximum of 10 is suppored!\n'.format(self.n))

        # Convert degrees to radians
        sweep = [x * np.pi/180 for x in sweep]
        dihedral = [x * np.pi/180 for x in dihedral]
        twist = [x * np.pi/180 for x in twist]

        # Compute wing shape parameters
        self.compShape(span, taper, rootChord)

        # Create airfoil points and indices
        self.initData(filenames, span, twist, sweep, dihedral, offset)

    def compShape(self, span, taper, rootChord):
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
        self.AR = 2 * self.b*self.b/self.S

    def initData(self, filenames, span, twist, sweep, dihedral, offset):
        """Read, transform and store airfoil points, and define numbering
        """
        self.pts = []
        self.ptsN = []
        # read and store coordinates (10 airfoils of max. 499 points each: 1-5000)
        for i in range(0, self.n):
            aPts = self.read(filenames[i])
            size = aPts.shape[0]
            aPts = np.hstack((aPts, self.spanPos[i]*np.ones([size,1])))
            aPts[:,[1,2]] = np.fliplr(aPts[:,[1,2]])
            aIdx = np.arange(i*500+1, i*500+1+size)
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
            self.pts[i][:, 2] = self.pts[i][:, 2] + sum(np.tan(dihedral[0:i])*span[0:i])
        for i in range(0, self.n):
            # apply offset
            self.pts[i][:, 0] += offset[0] # x
            self.pts[i][:, 2] += offset[1] # z
        # get separation points numbering
        self.sptsNl = []
        self.sptsNg = [] # todo: remove since global index can be recovered from local index: ptsN[local]
        for i in range(0, self.n):
            numb = self.specPts(i)
            self.sptsNl.append(numb)
            self.sptsNg.append(numb+i*500+1)

        # define line numbering (6 lines per airfoil: 1-60)
        self.linaN = []
        for i in range(0, self.n):
            self.linaN.append(np.arange(i*6+1, (i+1)*6+1))
        # define line numbering (6 lines per wing station: 61-114)
        self.linpN = []
        for i in range(0, self.n-1):
            self.linpN.append(np.arange(i*6+61, (i+1)*6+61))
        # define surface numbering (6 per wing station: 1-55)
        self.surN = []
        for i in range(0, self.n-1):
            self.surN.append(np.arange(i*6+1, (i+1)*6+1))

    def specPts(self, idx):
        """Find (local) index and location of separation points
        """
        # fraction of the chord defining separation points (could be given as user-def params)
        sepFwd = 0.3
        sepAft = 0.9
        # trailing and leading edge
        te = 0
        le = np.argmin(self.pts[idx][:,0])
        # find upper separations
        teU = np.argmin(np.abs(self.pts[idx][te:le,0]-self.pts[idx][le,0] - sepAft*self.chord[idx]))
        leU = np.argmin(np.abs(self.pts[idx][te:le,0]-self.pts[idx][le,0] - sepFwd*self.chord[idx]))
        # find lower separations
        leL = le + np.argmin(np.abs(self.pts[idx][le:-1,0]-self.pts[idx][le,0] - sepFwd*self.chord[idx]))
        teL = le + np.argmin(np.abs(self.pts[idx][le:-1,0]-self.pts[idx][le,0] - sepAft*self.chord[idx]))

        return np.array([te , teU, leU, le, leL, teL])

    def read(self,fname):
        """Read data from file and stroe in matrix
        """
        _file = open(fname, 'r')
        label = next(_file).split(',')
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
        file.write('// Full-wing aspect ratio: {0:f}\n'.format(self.AR))
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
            # TE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msTe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][0]], self.pts[i][self.sptsNl[i][0],0], self.pts[i][self.sptsNl[i][0],1], self.pts[i][self.sptsNl[i][0],2], i))
            for j in range(self.sptsNl[i][0]+1, self.sptsNl[i][1]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            # upper TE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msTe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][1]], self.pts[i][self.sptsNl[i][1],0], self.pts[i][self.sptsNl[i][1],1], self.pts[i][self.sptsNl[i][1],2], i))
            for j in range(self.sptsNl[i][1]+1, self.sptsNl[i][2]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            # upper LE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msLe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][2]], self.pts[i][self.sptsNl[i][2],0], self.pts[i][self.sptsNl[i][2],1], self.pts[i][self.sptsNl[i][2],2], i))
            for j in range(self.sptsNl[i][2]+1, self.sptsNl[i][3]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            # LE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msLe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][3]], self.pts[i][self.sptsNl[i][3],0], self.pts[i][self.sptsNl[i][3],1], self.pts[i][self.sptsNl[i][3],2], i))
            for j in range(self.sptsNl[i][3]+1, self.sptsNl[i][4]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            # lower LE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msLe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][4]], self.pts[i][self.sptsNl[i][4],0], self.pts[i][self.sptsNl[i][4],1], self.pts[i][self.sptsNl[i][4],2], i))
            for j in range(self.sptsNl[i][4]+1, self.sptsNl[i][5]):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
            # lower TE
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msTe{4:d}}};\n'.format(self.ptsN[i][self.sptsNl[i][5]], self.pts[i][self.sptsNl[i][5],0], self.pts[i][self.sptsNl[i][5],1], self.pts[i][self.sptsNl[i][5],2], i))
            for j in range(self.sptsNl[i][5]+1, self.ptsN[i].shape[0]-1):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
        file.write('\n')
        file.close()

    def writeLines(self, fname):
        """Write wing lines
        """
        file = open(fname, 'a')
        file.write('// --- Wing lines ---\n')
        # airfoil lines
        for i in range(0, self.n):
            file.write('// -- Airfoil {0:d}\n'.format(i))
            for j in range(0, self.linaN[i].shape[0]-1):
                file.write('Spline({0:d}) = {{'.format(self.linaN[i][j]))
                for k in range(self.sptsNg[i][j], self.sptsNg[i][j+1]):
                    file.write('{0:d}, '.format(k))
                file.write('{0:d}'.format(self.sptsNg[i][j+1]))
                file.write('};\n')
            file.write('Spline({0:d}) = {{'.format(self.linaN[i][-1]))
            for k in range(self.sptsNg[i][self.linaN[i].shape[0]-1], self.sptsNg[i][0]+self.ptsN[i].shape[0]-1):
                file.write('{0:d}, '.format(k))
            file.write('{0:d}'.format(self.sptsNg[i][0]))
            file.write('};\n')
        # planform lines
        for i in range(0, self.n-1):
            file.write('// -- Planform {0:d}\n'.format(i))
            for j in range(0, self.linpN[i].shape[0]):
                file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linpN[i][j], self.sptsNg[i][j], self.sptsNg[i+1][j]))
        file.write('\n')
        file.close()

    def writeSurfaces(self, fname):
        """Write wing line loops and surfaces
        """
        file = open(fname, 'a')
        file.write('// --- Wing line loops and surfaces ---\n')
        for i in range(0, self.n-1):
            file.write('// -- Planform {0:d}\n'.format(i))
            for j in range(0, self.surN[i].shape[0]):
                file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[i][j], self.linaN[i][j], self.linpN[i][np.mod(j+1,self.linpN[i].shape[0])], -self.linaN[i+1][j], -self.linpN[i][j]))
            for j in range(0, self.surN[i].shape[0]):
                file.write('Surface({0:d}) = {{-{0:d}}};\n'.format(self.surN[i][j]))
        file.write('\n')
        file.close()

    def writePhysical(self, fname):
        """Write wing physical groups
        """
        import os
        file = open(fname, 'a')
        file.write('// --- Wing physical groups ---\n')
        file.write('Physical Surface("wing") = {')
        for i in range(0, self.n-1):
            for j in range(0, 3):
                file.write('{0:d},'.format(self.surN[i][j]))
        file.seek(0, os.SEEK_END) # first seek end of file; f.seek(0, 2) is legal
        file.seek(file.tell() - 1, os.SEEK_SET) # then go backward
        file.truncate()
        file.write('};\n')
        file.write('Physical Surface("wing_") = {')
        for i in range(0, self.n-1):
            for j in range(3, 6):
                file.write('{0:d},'.format(self.surN[i][j]))
        file.seek(0, os.SEEK_END)
        file.seek(file.tell() - 1, os.SEEK_SET)
        file.truncate()
        file.write('};\n')
        file.write('\n')
        file.close()
