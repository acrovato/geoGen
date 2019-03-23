#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
## @package FPE grid creator
#
# Create a rectangular unstructured tetrahedral grid around a wing
# to be meshed with gmsh for Flow Full Potential solver
# Adrien Crovato

import numpy as np
import wing as w
import box as b

## Manage interface between modules
#
# Adrien Crovato
class Interface:
    # Interface is the manager that will
    # handle data from modules (wing, box)
    # create the .geo logic
    # write the .geo