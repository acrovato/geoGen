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

## Generic wingtip class
#
# Adrien Crovato
class Tip:
    def __init__(self, _wing):
        self.wing = _wing

        self.initData()

    def initData(self):
        """Initialize data, define numbering
        """
        # build mean line
        self.pts = [(np.zeros([(self.wing.pts[-1].shape[0]-3)//2,3]))]
        for i in range(0, self.pts[0].shape[0]):
            self.pts[0][i,:] = np.array([0.5*(self.wing.pts[-1][1+i,0]+self.wing.pts[-1][-2-i,0]), self.wing.pts[-1][0,1], 0.5*(self.wing.pts[-1][1+i,2]+self.wing.pts[-1][-2-i,2])])

        # define point numbering (max. (499-3)/2 more points: 5101-5348)
        self.ptsN = [np.arange(5101, 5101+self.pts[0].shape[0])]
        # define line numbering (7 lines: 121-127)
        self.linN = [np.arange(121, 128)]
        # define surface numbering (6 surfaces: 71-76)
        self.surN = [np.arange(71, 77)]

        # fraction of the chord defining separation points (could be given as user-def params)
        sepFwd = 0.3
        sepAft = 0.9
        # find and store separation poins
        orgn = np.min(self.wing.pts[-1][:,0])
        aft = np.argmin(np.abs(self.pts[0][:,0]-orgn - sepAft*self.wing.chord[-1]))
        fwd = np.argmin(np.abs(self.pts[0][:,0]-orgn - sepFwd*self.wing.chord[-1]))
        self.sptsN = [np.array([aft, fwd])]

## Handle cutoff wingtip data
#
# Adrien Crovato
class CTip(Tip):
    def __init__(self, _wing):

        Tip.__init__(self, _wing)

    def writeInfo(self,fname):
        """Write wing geometrical parameters
        """
        file = open(fname, 'a')
        file.write('// --- Wingtip geometry ---\n')
        file.write('// Cutoff wingtip\n')
        file.write('\n')
        file.close()

    def writePoints(self,fname):
        """Write wing points
        """
        file = open(fname, 'a')
        file.write('// --- Wingtip points ---\n')
        for i in range(0, self.sptsN[0][0]):
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[0][i], self.pts[0][i,0], self.pts[0][i,1], self.pts[0][i,2]))
        file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msTe{4:d}}};\n'.format(self.sptsN[0][0]+self.ptsN[0][0], self.pts[0][self.sptsN[0][0],0], self.pts[0][self.sptsN[0][0],1],self.pts[0][self.sptsN[0][0],2], self.wing.n-1))
        for i in range(self.sptsN[0][0]+1, self.sptsN[0][1]):
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[0][i], self.pts[0][i,0], self.pts[0][i,1], self.pts[0][i,2]))
        file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},gr{4:d}*msLe{4:d}}};\n'.format(self.sptsN[0][1]+self.ptsN[0][0], self.pts[0][self.sptsN[0][1],0], self.pts[0][self.sptsN[0][1],1],self.pts[0][self.sptsN[0][1],2], self.wing.n-1))
        for i in range(self.sptsN[0][1]+1, self.ptsN[0].shape[0]):
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[0][i], self.pts[0][i,0], self.pts[0][i,1], self.pts[0][i,2]))
        file.write('\n')
        file.close()

    def writeLines(self, fname):
        """Write wing lines
        """
        file = open(fname, 'a')
        file.write('// --- Wingtip lines ---\n')
        # midlines
        file.write('Spline({0:d}) = {{{1:d},'.format(self.linN[0][0], self.wing.ptsN[-1][self.wing.sptsNl[-1][0]]))
        for i in range(0, self.sptsN[0][0]):
            file.write('{0:d}, '.format(self.ptsN[0][i]))
        file.write('{0:d}}};\n'.format(self.ptsN[0][self.sptsN[0][0]]))
        file.write('Spline({0:d}) = {{{1:d},'.format(self.linN[0][1], self.ptsN[0][self.sptsN[0][0]]))
        for i in range(self.sptsN[0][0]+1, self.sptsN[0][1]):
            file.write('{0:d}, '.format(self.ptsN[0][i]))
        file.write('{0:d}}};\n'.format(self.ptsN[0][self.sptsN[0][1]]))
        file.write('Spline({0:d}) = {{{1:d},'.format(self.linN[0][2], self.ptsN[0][self.sptsN[0][1]]))
        for i in range(self.sptsN[0][1]+1, self.ptsN[0].shape[0]):
            file.write('{0:d}, '.format(self.ptsN[0][i]))
        file.write('{0:d}}};\n'.format(self.wing.ptsN[-1][self.wing.sptsNl[-1][3]]))
        # to-midlines
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[0][3], self.wing.sptsNg[-1][1], self.ptsN[0][self.sptsN[0][0]]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[0][4], self.wing.sptsNg[-1][2], self.ptsN[0][self.sptsN[0][1]]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[0][5], self.wing.sptsNg[-1][4], self.ptsN[0][self.sptsN[0][1]]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[0][6], self.wing.sptsNg[-1][5], self.ptsN[0][self.sptsN[0][0]]))
        file.write('\n')
        file.close()

    def writeSurfaces(self, fname):
        """Write wing line loops and surfaces
        """
        file = open(fname, 'a')
        file.write('// --- Wingtip line loops and surfaces ---\n')
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d}}};\n'.format(self.surN[0][0], self.wing.linaN[-1][0], self.linN[0][3], -self.linN[0][0]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][1], self.wing.linaN[-1][1], self.linN[0][4], -self.linN[0][1], -self.linN[0][3]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d}}};\n'.format(self.surN[0][2], self.wing.linaN[-1][2], -self.linN[0][2], -self.linN[0][4]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d}}};\n'.format(self.surN[0][3], self.wing.linaN[-1][3], self.linN[0][5], self.linN[0][2]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][4], self.wing.linaN[-1][4], self.linN[0][6], self.linN[0][1], -self.linN[0][5]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d}}};\n'.format(self.surN[0][5], self.wing.linaN[-1][5], self.linN[0][0], -self.linN[0][6]))
        for i in range(0, self.surN[0].shape[0]):
                file.write('Surface({0:d}) = {{-{0:d}}};\n'.format(self.surN[0][i]))
        file.write('\n')
        file.close()

    def writePhysical(self, fname):
        """Write wing physical groups
        """
        import os
        file = open(fname, 'a')
        file.write('// --- Wingtip physical groups ---\n')
        file.write('Physical Surface("wing") += {')
        for j in range(0, 3):
            file.write('{0:d},'.format(self.surN[0][j]))
        file.seek(0, os.SEEK_END)
        file.seek(file.tell() - 1, os.SEEK_SET)
        file.truncate()
        file.write('};\n')
        file.write('Physical Surface("wing_") += {')
        for j in range(3, 6):
            file.write('{0:d},'.format(self.surN[0][j]))
        file.seek(0, os.SEEK_END)
        file.seek(file.tell() - 1, os.SEEK_SET)
        file.truncate()
        file.write('};\n')
        file.write('\n')
        file.close()

## Handle rounded wingtip data
#
# Adrien Crovato
class RTip(Tip):
    def __init__(self, _wing):
        
        Tip.__init__(self, _wing)

        raise Exception('RTip: rounded wingtip not implemented yet!\n')

        # dummy centers
        # offset the poitn with a linear law
        # build bezier
