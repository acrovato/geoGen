#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package GeoGen (CFD basic grid creator)
#
# Create an unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow or SU2 CFD solvers
# Adrien Crovato

import numpy as np

## Handle box data
#
# Adrien Crovato
class Box:
    def __init__(self, xO, xF, yF, zO, zF, _wing, _tip, _wake):
        self.wing = _wing
        self.tip = _tip
        self.wake = _wake

        self.initData(xO, xF, yF, zO, zF)

    def initData(self, xO, xF, yF, zO, zF):
        """Initialize data, define numbering
        """
        self.pts = [np.array([[xF, 0., zF],
                              [xO, 0., zF],
                              [xO, 0., zO],
                              [xF, 0., zO]]),
                    np.array([[xF, yF, zF],
                              [xO, yF, zF],
                              [xO, yF, zO],
                              [xF, yF, zO]])]
        self.ptsN = [np.arange(5000, 5004), np.arange(5004, 5008)]

        # line numbering (2*6 x lines: 191-203) AND (4 y lines: 205-209)
        self.linxN = [np.arange(191, 197), np.arange(197, 204)]
        self.linyN = [np.arange(205, 210)]

        # surface numbering (10 surfaces: 111-120)
        self.surN = [np.arange(111, 121)]

    def writeInfo(self, fname):
        """Write box geometrical parameters
        """
        file = open(fname, 'a')
        file.write('// --- Domain geometry ---\n')
        file.write('// Box length: {0:f}\n'.format(self.pts[0][0,0]-self.pts[0][1,0]))
        file.write('// Box width: {0:f}\n'.format(self.pts[1][0,1]))
        file.write('// Box height: {0:f}\n'.format(self.pts[0][0,2]-self.pts[0][3,2]))
        file.write('\n')
        file.close()

    def writePoints(self, fname):
        """Write box points
        """
        file = open(fname, 'a')
        file.write('// --- Box points ---\n')
        for i in range(0, 2):
            for j in range(0,4):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msF}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
        file.write('\n')

    def writeLines(self, fname):
        """Write box lines
        """
        file = open(fname, 'a')
        file.write('// --- Box lines ---\n')
        file.write('// -- Symmetry\n')
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][0], self.wake.ptsN[0][0], self.ptsN[0][0]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][1], self.ptsN[0][0], self.ptsN[0][1]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][2], self.ptsN[0][1], self.wake.ptsN[0][-1]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][3], self.wake.ptsN[0][-1], self.ptsN[0][2]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][4], self.ptsN[0][2], self.ptsN[0][3]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[0][5], self.ptsN[0][3], self.wake.ptsN[0][0]))
        file.write('// -- Back\n')
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][0], self.wake.ptsN[0][self.wing.n], self.ptsN[1][0]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][1], self.ptsN[1][0], self.ptsN[1][1]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][2], self.ptsN[1][1], self.wake.ptsN[0][self.wing.n+5]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][3], self.wake.ptsN[0][self.wing.n+5], self.ptsN[1][2]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][4], self.ptsN[1][2], self.ptsN[1][3]))
        file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linxN[1][5], self.ptsN[1][3], self.wake.ptsN[0][self.wing.n]))
        file.write('// -- Transverse\n')
        for i in range(0, 4):
            file.write('Line({0:d}) = {{{1:d},{2:d}}};\n'.format(self.linyN[0][i], self.ptsN[0][i], self.ptsN[1][i]))
        file.write('\n')
        file.close()

    def writeSurfaces(self, fname):
        """Write box surfaces
        """
        file = open(fname, 'a')
        file.write('// --- Box surfaces ---\n')
        # line loops
        file.write('// -- Symmetry\n')
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d},{5:d},{6:d},{7:d},{8:d}}};\n'.format(self.surN[0][0], self.wing.linaN[0][0], self.wing.linaN[0][1], self.wing.linaN[0][2], self.wake.linN[1][-1], -self.linxN[0][2], -self.linxN[0][1], -self.linxN[0][0], -self.wake.linN[1][0]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d},{5:d},{6:d},{7:d},{8:d}}};\n'.format(self.surN[0][1], self.wing.linaN[0][3], self.wing.linaN[0][4], self.wing.linaN[0][5], self.wake.linN[1][0], -self.linxN[0][5], -self.linxN[0][4], -self.linxN[0][3], -self.wake.linN[1][-1]))
        file.write('// -- Downstream\n')
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},'.format(self.surN[0][2], self.linyN[0][0], -self.linxN[1][0]))
        for i in range(0, self.wing.n):
            file.write('{0:d},'.format(-self.wake.linN[0][self.wing.n-i-1]))
        file.write('{0:d}}};\n'.format(self.linxN[0][0]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},'.format(self.surN[0][3], -self.linyN[0][3], self.linxN[0][5]))
        for i in range(0, self.wing.n):
            file.write('{0:d},'.format(self.wake.linN[0][i]))
        file.write('{0:d}}};\n'.format(-self.linxN[1][5]))
        file.write('// -- Farfield\n')
        # upstream
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},'.format(self.surN[0][4], -self.linyN[0][1], self.linxN[0][2]))
        for i in range(0, self.wing.n):
            file.write('{0:d},'.format(-self.wake.linN[0][-i-1]))
        file.write('{0:d}}};\n'.format(-self.linxN[1][2]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},'.format(self.surN[0][5], self.linyN[0][2], -self.linxN[1][3]))
        for i in range(0, self.wing.n):
            file.write('{0:d},'.format(self.wake.linN[0][-self.wing.n+i]))
        file.write('{0:d}}};\n'.format(self.linxN[0][3]))
        # back
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d},{5:d},{6:d},{7:d},{8:d}}};\n'.format(self.surN[0][6], self.linxN[1][0], self.linxN[1][1], self.linxN[1][2], -self.wake.linN[0][self.wing.n+4], -self.wake.linN[0][self.wing.n+3], -self.wake.linN[0][self.wing.n+2], -self.wake.linN[0][self.wing.n+1], -self.wake.linN[0][self.wing.n]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d},{5:d},{6:d},{7:d},{8:d}}};\n'.format(self.surN[0][7], self.linxN[1][3], self.linxN[1][4], self.linxN[1][5], self.wake.linN[0][self.wing.n], self.wake.linN[0][self.wing.n+1], self.wake.linN[0][self.wing.n+2], self.wake.linN[0][self.wing.n+3], self.wake.linN[0][self.wing.n+4]))
        # top and bottom
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][8], self.linxN[0][1], self.linyN[0][1], -self.linxN[1][1], -self.linyN[0][0]))
        file.write('Line Loop({0:d}) = {{{1:d},{2:d},{3:d},{4:d}}};\n'.format(self.surN[0][9], self.linxN[0][4], self.linyN[0][3], -self.linxN[1][4], -self.linyN[0][2]))
        # surfaces
        for i in range(0, self.surN[0].shape[0]):
            file.write('Plane Surface({0:d}) = {{{0:d}}};\n'.format(self.surN[0][i]))
        file.write('\n')
        file.close()

    def writeVolumes(self, fname):
        """Write computational volume
        """
        file = open(fname, 'a')
        file.write('// --- Computational volumes ---\n')
        file.write('// -- Upper\n')
        # surface loops
        file.write('Surface Loop({0:d}) = {{'.format(1))
        for i in range(0, self.wing.n-1):
            for j in range(0, 3):
                file.write('{0:d},'.format(self.wing.surN[i][j]))
        for j in range(0, 3):
            file.write('{0:d},'.format(self.tip.surN[0][j]))
        for j in range(0, self.wake.surN[0].shape[0]):
            file.write('{0:d},'.format(self.wake.surN[0][j]))
        for j in range(0, 4):
            file.write('{0:d},'.format(self.surN[0][2*j]))
        file.write('{0:d}}};\n'.format(self.surN[0][8]))
        file.write('// -- Lower\n')
        file.write('Surface Loop({0:d}) = {{'.format(2))
        for i in range(0, self.wing.n-1):
            for j in range(3, 6):
                file.write('{0:d},'.format(self.wing.surN[i][j]))
        for j in range(3, 6):
            file.write('{0:d},'.format(self.tip.surN[0][j]))
        for j in range(0, self.wake.surN[0].shape[0]):
            file.write('{0:d},'.format(self.wake.surN[0][j]))
        for j in range(0, 4):
            file.write('{0:d},'.format(self.surN[0][2*j+1]))
        file.write('{0:d}}};\n'.format(self.surN[0][9]))
        # volumes
        file.write('Volume({0:d}) = {{{0:d}}};\n'.format(1))
        file.write('Volume({0:d}) = {{{0:d}}};\n'.format(2))
        file.write('\n')
        file.close()

    def writePhysical(self, fname):
        """Write box physical groups
        """
        file = open(fname, 'a')
        file.write('// --- Box physical groups ---\n')
        file.write('Physical Surface("symmetry") = {{{0:d}}};\n'.format(self.surN[0][0]))
        file.write('Physical Surface("symmetry_") = {{{0:d}}};\n'.format(self.surN[0][1]))
        file.write('Physical Surface("downstream") = {{{0:d}}};\n'.format(self.surN[0][2]))
        file.write('Physical Surface("downstream_") = {{{0:d}}};\n'.format(self.surN[0][3]))
        file.write('Physical Surface("farfield") = {{{0:d},{1:d},{2:d},{3:d},{4:d},{5:d}}};\n'.format(self.surN[0][4],self.surN[0][5],self.surN[0][6],self.surN[0][7],self.surN[0][8],self.surN[0][9]))
        file.write('Physical Volume("field") = {{{0:d}}};\n'.format(1))
        file.write('Physical Volume("field_") = {{{0:d}}};\n'.format(2))
        file.write('\n')
        file.close()
