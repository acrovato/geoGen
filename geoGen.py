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

import wing as w
import tip as t
import wake as wk
import domain as d

def main(_module, _output):
    # Get config
    p = getConfig(_module)

    # Create wing, wingtip, wake and bounding domain
    wing = w.Wing(p['airfName'], p['span'], p['taper'], p['sweep'], p['dihedral'], p['twist'], p['rootChord'], p['offset'])
    if p['coWingtip']:
        tip = t.CTip(wing)
    else:
        tip = t.RTip(wing)
    if p['domType'] == 'box':
        wake = wk.Wake(p['xoBox'], p['xfBox'], p['yfBox'], p['nSlope'], wing, tip)
        dom = d.Box(p['xoBox'], p['xfBox'], p['yfBox'], p['zoBox'], p['zfBox'], wing, tip, wake)
    elif p['domType'] == 'sphere':
        wake = wk.GWake()
        dom = d.Sphere(p['rSphere'], wing, tip)
    else:
        raise Exception('"domType" parameter can be either "box" or "sphere", but', p['domType'], 'was given!\n')

    # Switch to workspace and write
    createWdir()
    outFile = _output
    # misc
    writeHeader(outFile, _module)
    wing.writeInfo(outFile)
    tip.writeInfo(outFile)
    dom.writeInfo(outFile)
    wing.writeOpts(outFile)
    dom.writeOpts(outFile)
    # points
    wing.writePoints(outFile)
    tip.writePoints(outFile)
    wake.writePoints(outFile)
    dom.writePoints(outFile)
    # lines
    wing.writeLines(outFile)
    tip.writeLines(outFile)
    wake.writeLines(outFile)
    dom.writeLines(outFile)
    # surfaces
    wing.writeSurfaces(outFile)
    tip.writeSurfaces(outFile)
    wake.writeSurfaces(outFile)
    dom.writeSurfaces(outFile)
    # volumes
    dom.writeVolumes(outFile)
    # physical
    wing.writePhysical(outFile)
    tip.writePhysical(outFile)
    wake.writePhysical(outFile)
    dom.writePhysical(outFile)
    # mesh options
    writeOpts(outFile, tip.surN)

    # Printout
    printInfo(outFile)   

    # eof
    print('')

def getConfig(_module):
    # Get prarmeters from config file
    import os, sys, ntpath
    sys.path.append(os.path.abspath(os.path.dirname(_module))) # tmp append module directory to pythonpath
    module = __import__(ntpath.basename(_module)) # import config as module
    p = module.getParams()
    # Fix path
    for i in range(0, len(p['airfName'])):
        p['airfName'][i] = os.path.join(os.path.abspath(os.path.dirname(_module)), p['airfPath'],p['airfName'][i])
    return p

def createWdir():
    import os
    wdir = os.path.join(os.getcwd(), 'workspace')
    if not os.path.isdir(wdir):
        print("creating", wdir)
        os.makedirs(wdir)
    os.chdir(wdir)

def writeHeader(fname, _module):
    import ntpath
    file = open(fname, 'w')
    file.write('/******************************************/\n')
    file.write('/* Gmsh geometry for {0:>20s} */\n'.format(ntpath.basename(_module)))
    file.write('/* Generated by      {0:>20s} */\n'.format(ntpath.basename(__file__)))
    file.write('/* Adrien Crovato                         */\n')
    file.write('/* ULiege, 2018-2019                      */\n')
    file.write('/******************************************/\n\n')
    file.close()

def writeOpts(fname, tipSur):
    """Write misc options
    """
    import os
    file = open(fname, 'a')
    file.write('// --- Misc Meshing options ---\n')
    file.write('Mesh.Algorithm = 5; // Delaunay\n')
    file.write('MeshAlgorithm Surface {{{0:d},{1:d}}} = 1; // Mesh-adapt\n'.format(tipSur[0][2], tipSur[0][3]))
    file.write('Mesh.Algorithm3D = 2; // New Delaunay\n')
    file.write('Mesh.OptimizeNetgen = 1;\n')
    file.write('Mesh.Smoothing = 10;\n')
    file.write('Mesh.SmoothNormals = 1;\n')
    file.write('\n')
    file.close()

def printInfo(fname):
    """Print info
    """
    import os
    print('*' * 79)
    print('* geoGen')
    print('* Adrien Crovato')
    print('* ULiege, 2018-2020')
    print('* Distributed under Apache license 2.0')
    print('*' * 79)
    print(os.path.abspath(os.path.join(os.getcwd(), fname)), 'has been successfully written!')
    print('Visual file check in gmsh recommended before further use!')
    print('*' * 79)

if __name__ == "__main__":
    # Arguments parser
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='input config .py file')
    parser.add_argument('-o', dest='out', help='output .geo file', default='grid.geo')
    args = parser.parse_args()

    main(args.file[:-3], args.out)
