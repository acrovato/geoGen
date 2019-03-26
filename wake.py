#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package FPE grid creator
#
# Create a rectangular unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow Full Potential solver
# Adrien Crovato

import numpy as np

## Manage wake data
#
# Adrien Crovato
class Wake:
    def __init__(self, xO, xF, yF, nSlope, _wing, _tip):
        self.wing = _wing
        self.tip = _tip

        self.midplane(xO, xF, yF, nSlope)

    def midplane(self, xO, xF, yF, nSlope):
        """Compute wake points and define numbering
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