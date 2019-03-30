#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' 
Copyright 2019 University of Liege

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

## Generic wake class (does nothing)
#
# Adrien Crovato
class GWake:
    def __init__(self):
        """Desc.
        """

    def writePoints(self, fname):
        """Desc.
        """

    def writeLines(self, fname):
        """Desc.
        """

    def writeSurfaces(self, fname):
        """Desc.
        """

    def writePhysical(self, fname):
        """Desc.
        """

## Manage wake data
#
# Adrien Crovato
class Wake(GWake):
    def __init__(self, xO, xF, yF, nSlope, _wing, _tip):
        self.wing = _wing
        self.tip = _tip

        self.initData(xO, xF, yF, nSlope)

    def initData(self, xO, xF, yF, nSlope):
        """Initialize data, define numbering
        """
        n = self.wing.n
        pts = np.zeros([n*2+6,3])
        for i in range(0, n):
            slopeU = (self.wing.pts[i][nSlope,2] - self.wing.pts[i][0,2]) / (self.wing.pts[i][nSlope,0] - self.wing.pts[i][0,0])
            slopeL = (self.wing.pts[i][-nSlope-1,2] - self.wing.pts[i][0,2]) / (self.wing.pts[i][-nSlope-1,0] - self.wing.pts[i][0,0])
            zW = self.wing.pts[i][0,2] + 0.5*(slopeU+slopeL) * (xF - self.wing.pts[i][0,0])
            pts[i,:] = np.array([xF, self.wing.spanPos[i], zW])
        pts[n,:] = np.array([xF, yF, self.wing.pts[-1][0,2]])
        pts[n+1,:] = np.array([self.wing.pts[-1][0,0], yF, self.wing.pts[-1][0,2]])
        pts[n+2,:] = np.array([self.tip.pts[0][self.tip.sptsN[0][0],0], yF, self.tip.pts[0][self.tip.sptsN[0][0],2]])
        pts[n+3,:] = np.array([self.tip.pts[0][self.tip.sptsN[0][1],0], yF, self.tip.pts[0][self.tip.sptsN[0][1],2]])
        pts[n+4,:] = np.array([self.wing.pts[-1][self.wing.sptsNl[-1][3],0], yF, self.wing.pts[-1][self.wing.sptsNl[-1][3],2]])
        pts[n+5,:] = np.array([xO, yF, self.wing.pts[-1][self.wing.sptsNl[-1][3],2]])
        for i in range(0, n):
            pts[n+6+i,:] = np.array([xO, self.wing.spanPos[n-1-i], self.wing.pts[n-i-1][self.wing.sptsNl[n-i-1][3],2]])
        self.pts = [pts]
        self.ptsN = [np.arange(5351, 5351+n*2+6)]

        # line numbering (max. 2*9+7 lines: 131-155) AND (max 2*10+4 lines: 161-184)
        self.linN = [np.arange(131, 131+2*(n-1)+7), np.arange(161, 161+2*(n)+4)]

        # surface numbering (max. 2*9+5: 81-103)
        self.surN = [np.arange(81, 81+2*(n-1)+5)]

    def writePoints(self, fname):
        """Write wake points
        """
        file = open(fname, 'a')
        file.write('// --- Wake points ---\n')
        for i in range(0, self.pts[0].shape[0]):
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msF}};\n'.format(self.ptsN[0][i], self.pts[0][i,0], self.pts[0][i,1], self.pts[0][i,2]))
        file.write('\n')
        file.close()

    def writeLines(self, fname):
        """Write wake lines
        """
        file = open(fname, 'a')
        file.write('// --- Wake lines ---\n')
        # domain lines
        for i in range(0, self.linN[0].shape[0]):
            file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[0][i], self.ptsN[0][i], self.ptsN[0][i+1]))
        # wing-to-domain lines
        for i in range(0, self.wing.n):
            file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][i], self.wing.sptsNg[i][0], self.ptsN[0][i]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][self.wing.n], self.wing.sptsNg[-1][0], self.ptsN[0][self.wing.n+1]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][self.wing.n+1], self.tip.ptsN[0][self.tip.sptsN[0][0]], self.ptsN[0][self.wing.n+2]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][self.wing.n+2], self.tip.ptsN[0][self.tip.sptsN[0][1]], self.ptsN[0][self.wing.n+3]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][self.wing.n+3], self.wing.sptsNg[-1][3], self.ptsN[0][self.wing.n+4]))
        for i in range(0, self.wing.n):
            file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linN[1][-self.wing.n+i], self.wing.sptsNg[-i-1][3], self.ptsN[0][-self.wing.n+i]))
        file.write('\n')
        file.close()

    def writeSurfaces(self, fname):
        """Write wake surfaces
        """
        file = open(fname, 'a')
        file.write('// --- Wake surfaces ---\n')
        # line loops
        file.write('// -- Wake\n')
        for i in range(0, self.wing.n-1):
            file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][i], self.linN[1][i], self.linN[0][i], -self.linN[1][i+1], -self.wing.linpN[i][0]))
        file.write('// -- Side\n')
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][self.wing.n-1], self.linN[1][self.wing.n-1], self.linN[0][self.wing.n-1], self.linN[0][self.wing.n], -self.linN[1][self.wing.n]))
        for i in range(0, 3):
            file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][self.wing.n+i], self.linN[1][self.wing.n+i], self.linN[0][self.wing.n+i+1], -self.linN[1][self.wing.n+i+1], -self.tip.linN[0][i]))    
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][self.wing.n+3], self.linN[1][self.wing.n+3], self.linN[0][self.wing.n+4], self.linN[0][self.wing.n+5], -self.linN[1][self.wing.n+4]))
        file.write('// -- Front\n')
        for i in range(0, self.wing.n-1):
            file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][-self.wing.n+1+i], self.linN[1][-self.wing.n+i], self.linN[0][-self.wing.n+1+i], -self.linN[1][-self.wing.n+1+i], self.wing.linpN[-1-i][3]))
        # surfaces
        for i in range(0, self.surN[0].shape[0]):
            file.write('Surface({0:d}) = {{{0:d}}};\n'.format(self.surN[0][i]))
        file.write('\n')
        file.close()

    def writePhysical(self, fname):
        """Write wake physical groups
        """
        import os
        file = open(fname, 'a')
        file.write('// --- Wake physical groups ---\n')
        file.write('Physical Line("wakeTip") = {{{0:d}}};\n'.format(self.linN[1][self.wing.n-1]))
        file.write('Physical Line("teTip") = {{{0:d}'.format(self.linN[1][self.wing.n-1]))
        for i in range(0, self.wing.n-1):
            file.write('{0:d},'.format(self.wing.linpN[i][0]))
        file.seek(-1, os.SEEK_END)
        file.truncate()
        file.write('};\n')
        file.write('Physical Surface("wake") = {')
        for j in range(0, self.wing.n-1):
            file.write('{0:d},'.format(self.surN[0][j]))
        file.seek(-1, os.SEEK_END)
        file.truncate()
        file.write('};\n')
        file.write('\n')
        file.close()