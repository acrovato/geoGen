#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## Onera M6 configuration module for fpeGrid
#
# Adrien Crovato

def getParams():
    p = {}
    # Wing parameters (nP planforms for half-wing)
    p['airfPath'] = '../airfoils' # path pointing the airfoils directory (relative to this config file)
    p['airfName'] = ['oneraM6.dat', 'oneraM6.dat'] # names of file containing airfoil (Selig formatted) data (size: nP+1)
    p['span'] = [1.196] # span of each planform (size: nP)
    p['taper'] = [0.562] # taper of each planform (size: nP)
    p['sweep'] = [30] # leading edge sweep of each planform (size: nP) 
    p['dihedral'] = [0] # dihedral angle of each planform (size: nP)
    p['twist'] = [0, 0] # twist angle of each airfoil (size: nP+1)
    p['rootChord'] = 0.8059
    # Box
    # Wake

    return p
