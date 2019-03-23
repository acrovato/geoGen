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
    def __init__(self, rootChord, span):
        # Default size and origin
        self.xO = -3.5*rootChord
        self.xF = (3.5+1.)*rootChord
        self.yO = 0.
        self.yF = 2*span
        self.zO = -3.5*rootChord
        self.zF = 3.5*rootChord
        self.mshSize = 0.5*rootChord

        #self.PtsN = np.arange()

    def writeOpts(self, fname):
        """Write box gmsh options
        """
        file = open(fname, 'a')
        file.write('// --- Box options ---\n')
        file.write('DefineConstant[ xO = {{ {0:f}, Name "Box x-origin" }} ];\n'.format(self.xO))
        file.write('DefineConstant[ xF = {{ {0:f}, Name "Box x-end" }} ];\n'.format(self.xF))
        file.write('DefineConstant[ yO = {{ {0:f}, Name "Box y-end" }} ];\n'.format(self.yO))
        file.write('DefineConstant[ zO = {{ {0:f}, Name "Box z-origin" }} ];\n'.format(self.zO))
        file.write('DefineConstant[ zF = {{ {0:f}, Name "Box z-end" }} ];\n'.format(self.zF))
        file.write('DefineConstant[ msF = {{ {0:f}, Name "Farfield mesh size" }} ];\n'.format(self.mshSize))
        file.write('\n')
        file.close()

    def writePoints(self, fname):
        """Write box points
        """
        #file.write('Point({0:d}) = {{{1:f},{2:f},{3:f}}};\n'.format(self.ptsN[i][j], self.pts[i][j,0], self.pts[i][j,1], self.pts[i][j,2]))
        