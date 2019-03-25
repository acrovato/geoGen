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
    def __init__(self, xO, xF, yF, nSlope, spanPos):
        self.vertices(xO, xF, yF, nSlope, spanPos)

    def vertices(self, xO, xF, yF, nSlope, spanPos):
        """Compute wake points
        """
        n = len(spanPos)
        pts = np.zeros([n*2+6,3])
        for i in range(0, n):
            pts[i,:] = np.array([xF, spanPos[i], 0.])
        pts[n,:] = np.array([xF, yF, 0.])
        pts[n+1,:] = np.array([1.0, yF, 0.])
        pts[n+2,:] = np.array([1.2, yF, 0.])
        pts[n+3,:] = np.array([1.4, yF, 0.])
        pts[n+4,:] = np.array([1.8, yF, 0.])
        pts[n+5,:] = np.array([xO, yF, 0.])
        for i in range(0, n):
            pts[n+6+i,:] = np.array([xO, spanPos[n-1-i], 0.])
        self.pts = [pts]
        self.ptsN = [np.arange(5351, 5351+n*2+6)]

    def writePoints(self, fname):
        """Write wake points
        """
        file = open(fname, 'a')
        file.write('// --- Wake points ---\n')
        for i in range(0, self.pts[0].shape[0]):
            file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msF}};\n'.format(self.ptsN[0][i], self.pts[0][i,0], self.pts[0][i,1], self.pts[0][i,2]))
        file.write('\n')
        file.close()