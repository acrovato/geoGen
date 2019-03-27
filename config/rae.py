#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## Custom rae wing configuration module for geoGen
#
# Adrien Crovato

def getParams():
    p = {}
    # Wing parameters (nP planforms for half-wing)
    p['airfPath'] = '../airfoils' # path pointing the airfoils directory (relative to this config file)
    p['airfName'] = ['rae2822.dat', 'rae2822.dat', 'rae2822.dat'] # names of file containing airfoil (Selig formatted) data (size: nP+1)
    p['span'] = [.5, 2.] # span of each planform (size: nP)
    p['taper'] = [.8, .4] # taper of each planform (size: nP)
    p['sweep'] = [20., 20.] # leading edge sweep of each planform (size: nP) 
    p['dihedral'] = [2., 1.] # dihedral angle of each planform (size: nP)
    p['twist'] = [1., 0, -1.] # twist angle of each airfoil (size: nP+1)
    p['rootChord'] = 1.0 # root chord
    p['coWingtip'] = True # cut-off wingtip (not supported yet)
    # Box
    p['xoBox'] = -3.5*p['rootChord']
    p['xfBox'] = 4.5*p['rootChord']
    p['yfBox'] = 2*sum(p['span'])
    p['zoBox'] = -3.5*p['rootChord']
    p['zfBox'] = 3.5*p['rootChord']
    # Wake
    p['nSlope'] = 10 # number of airfoil TE points to compute wake slope

    return p
