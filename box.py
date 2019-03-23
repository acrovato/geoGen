#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package FPE grid creator
#
# Create a rectangular unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow Full Potential solver
# Adrien Crovato

import numpy as np

## Handle box data
#
# Adrien Crovato
class Box:
    def __init__(self, xO, xF, yF, zO, zF):
        self.corners(xO, xF, yF, zO, zF)

    def corners(self, xO, xF, yF, zO, zF):
        """Initialize default corner points
        """
        self.pts = [np.array([[xF, 0., zO],
                              [xO, 0., zO],
                              [xO, 0., zF],
                              [xF, 0., zF]]),
                    np.array([[xF, yF, zO],
                              [xO, yF, zO],
                              [xO, yF, zF],
                              [xF, yF, zF]])]
        self.ptsN = [np.arange(5000, 5004), np.arange(5004, 5008)]

    def writePoints(self, fname):
        """Write box points
        """
        file = open(fname, 'a')
        file.write('// --- Box points ---\n')
        for i in range(0, 2):
            for j in range(0,4):
                file.write('Point({0:d}) = {{{1:f},{2:f},{3:f},msF}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
        file.write('\n')