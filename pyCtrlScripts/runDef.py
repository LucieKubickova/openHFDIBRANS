#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python script to set up and run bread baking simulations according to 
# custom experiments

# IMPORTS===============================================================
import sys
from OF_caseClass import OpenFOAMCase
import numpy as np
# from blockMeshDictClassV8 import *
# from meshGeneration import *
# from myAddFcs import *
import re
import matplotlib.pyplot as plt
import os

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/openHFDIBRANS_noGap/' # -- base case for simulation
outFolder = '../ZZ_cases/openHFDIBRANS_noGap_def_IB/'
# expDir = os.path.join('..', 'Experiments2026') # -- when comparing experiments

# SIMULATION INPUTS=====================================================
# -- cluster settings
slurmRun = True
slurmRun = False

# -- decomposition
nProc = 1

# -- solvers
solverDef = 'solids4Foam' 

# -- numerics
nDefTIters = 10000  # -- number of solid/temperature iterations at each time level

# SCRIPT ITSELF (DO NOT EDIT)===========================================                       
# -- create OpenFOAMCase object to change values in dictionaries
baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase(baseCaseDir)
baseCase.changeOFCaseDir(outFolder)
baseCase.copyBaseCase()

# 1) constant/solidProperties
baseCase.setParameters(
    [
        ["constant/solidProperties", 'nCorrectors', str(nDefTIters), ""],
    ]
)

# 2) system/decomposePar
baseCase.setParameters(
    [
        ["system/decomposeParDict", "numberOfSubdomains", str(nProc), ""],
    ]
)

# -- run deformation simulation
if nProc > 1:
    if slurmRun:
        pass
    else:
        baseCase.runCommands(
            [
                'decomposePar > log.decomposePar',
                'foamJob -parallel -screen %s > log.%s' % (solverDef, solverDef)
            ]
        )
else:
    baseCase.runCommands(
        [
            '%s > log.%s' % (solverDef, solverDef)
        ]
    )