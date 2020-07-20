# GeoGen
Adrien Crovato  
ULiege, 2018-2020  

[![Run Status](https://api.shippable.com/projects/5c98c1405142dd0007ecd6c0/badge?branch=master)]()

## Features and limitations

### What is GeoGen?
GeoGen is a python code that generates an aerodynamic geometry which can then be meshed using gmsh and used in CFD solvers.  
GeoGen v2 and up is compatible with python 3, while previous versions are compatible with python 2.

### What can GeoGen do?
GeoGen currently supports the following configurations:
  - [x] Arbitrary isolated wing
    - [x] Sharp trailing edge
    - [x] Cutoff wingtip
    - [ ] Rounded wingtip
  - [ ] Generic fuselage
  - [ ] Horizontal tail
  - [ ] Multibody (e.g. several isolated wings)


GeoGen was primarly designed to be used with [SU2](https://github.com/su2code/SU2) and [waves](https://gitlab.uliege.be/am-dept/waves), but any solver interfaced with [gmsh](http://gmsh.info/) can be used!

## Usage
The script is run through the command line:
```sh
python geoGen.py path/to/config/file.py <-o path/to/output/file.geo>
```
If no output file is provided, a workspace directory will be created and the geometry will be stored inside as `grid.geo`.

The geometry is generated from a python file containing a dictionary of parameters. Examples are given in [config](config/) and the main options are summurized hereunder.

**Parameters**

Wing definition:
 - `airfPath`: relative path to the airfoils directory
 - `airfName`: array (size: nP+1) of names of file containing airfoil (Selig formatted) coordinates
 - `span`: array (size: nP) of span for each planform of the wing
 - `taper`: array (size: nP)) of taper of each planform of the wing
 - `sweep`: array (size: nP) of leading edge sweep of each planform of the wing 
 - `dihedral`: array (size: nP) of dihedral angle of each planform of the wing
 - `twist`: array (size: nP+1) of twist angle of each airfoil of the wing
 - `rootChord`: root chord (scalar) of the wing
 - `offset`: array of x and z offset (size: 2) applied to the leading edge of the root section
 - `coWingtip`: boolean, True for cutoof wingtip, Fasle for rounded wingtip (not supported yet)


Domain definition:
 - `domType`: string, box for box-shaped domain or shpere for shperical-shaped domain


if domain type is shpere, typical for Euler equations:
 - `rSphere`: radius of the sphere (scalar)


if domain type is a box (a wake will then be defined), typically for potential equations:
 - `xoBox`: x-coordinate (scalar) of the origin of the box 
 - `xfBox`: x-coordinate (scalar) of the end of the box
 - `yfBox`: y-coordinate (scalar) of the end of the box
 - `zoBox`: z-coordinate (scalar) of the origin of the box
 - `zfBox`: z-coordinate (scalar) of the end of the box
 - `nSlope`: number (scalar) of airfoil geometrical points counted from TE used to compute wake slope

