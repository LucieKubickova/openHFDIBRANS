#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python script to set up and run bread baking simulations according to 
# custom experiments

# IMPORTS===============================================================
import sys
from OF_caseClass import OpenFOAMCase
import numpy as np
# from meshGeneration import *
# from myAddFcs import *
import re
import matplotlib.pyplot as plt
import os
from mesher import *
from rewriteLambdaBoundary import *

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/openHFDIBRANS_flowTwoChannels/' # -- base case for simulation
# baseCaseDir = '../ZZ_cases/97_V5_onlyRect_par_DefIncl/'
# baseCaseDir = '../ZZ_cases/94_VIB_flow/'
outFolder = '../ZZ_cases/90_VIB_def_SIMP_par/'
# expDir = os.path.join('..', 'Experiments2026') # -- when comparing experiments

# WHAT TO RUN===========================================================
crtMesh = True
# crtMesh = False

flowRun = True  # -- simulation of flow
# flowRun = False

defRun = True   # -- simulation of deformation
# defRun = False

# SIMULATION INPUTS=====================================================
# -- cluster settings
slurmRun = True
slurmRun = False

# -- boundary and initial conditions
uIn = 0.001 # -- velocity
lScale = 0.005 # -- for computation of turbulence variables
Iper = 2.0
Iper /= 100
Cmu = 0.09

# -- decomposition
nProc = 2

# -- fields to copy for postProcessing
fldsToCopy = ["T", "U"]

# -- solvers
solverFlow = 'buoyantBoussinesqSimpleHFDIBRANSFoam' 
solverDef = 'solids4Foam' 

# -- numerics
nDefTIters = 5000  # -- number of solid/temperature iterations at each time level
flowEndTime = 1000

# SCRIPT ITSELF (DO NOT EDIT)===========================================                       
# -- create OpenFOAMCase object to change values in dictionaries
baseCase = OpenFOAMCase()
baseCase.loadOFCaseFromBaseCase(baseCaseDir)
baseCase.changeOFCaseDir(outFolder)
baseCase.copyBaseCase()

# -- calculate initial values for turbulence
k0 = 3/2*uIn**2*Iper**2
e0 = Cmu*k0**(3/2)/lScale
o0 = e0/Cmu/k0
n0 = k0/o0
values = [k0, e0, o0, n0]
names = ["k", "epsilon", "omega", "nut"]

# 1) boundary conditions
# -- turbulence variables
ls0 = os.listdir('%s/0.org/'%(baseCase.dir))
for name in names:
    if name in ls0:
        baseCase.setParameters(
            [
                ['0.org/%s'%(name), 'internalField', 'uniform %g' % (values[names.index(name)]),'']
            ]
        )
# -- velocity
baseCase.setParameters(
    [
        ['0.org/U', 'internalField', 'uniform (%s 0 0)'%(uIn), '']
    ]
)

# 2) constant/solidProperties
baseCase.setParameters(
    [
        ["constant/solidProperties", 'nCorrectors', str(nDefTIters), ""],
    ]
)

# 3) system/decomposePar
baseCase.setParameters(
    [
        ["system/decomposeParDict", "numberOfSubdomains", str(nProc), ""],
    ]
)

# 4) system/controlDict
baseCase.setParameters(
    [
        ["system/controlDict", "endTime", str(flowEndTime), ""],
    ]
)

# -- create mesh and lambda field
if crtMesh:
    createMeshLamBet(baseCase)
    baseCase.runCommands(
        [
            'blockMesh > log.blockMesh',
            'mkdir -p 0',
            'cp -rf 0.org/* 0/',
            'generateLambda > log.generateLambda',
            'createLambdaBoundary > log.createLambdaBoundary',
            # 'createPatch -overwrite -dict system/createPatchDictNewPatch > log.createPatch',
        ]
    )
    rewriteBoundary('%s/constant/polyMesh/boundary' % (baseCase.dir), '%s/0/lambda' % (baseCase.dir))
    baseCase.runCommands(
        [
            'correctLambdaBoundary > log.correctLambdaBoundary',
            'paraFoam -touch',
        ]
    )

# -- run flow simulation
if flowRun:
    if nProc > 1:
        if slurmRun:
            pass
        else:
            baseCase.runCommands(
                [
                    'decomposePar > log.decomposePar',
                    'foamJob -parallel -screen %s > log.%s' % (solverFlow, solverFlow)
                ]
            )
    else:
        baseCase.runCommands(
            [
                'chmod 755 ./Allclean && chmod 755 ./Allrun-serialAll',
                '%s > log.%s'% (solverFlow, solverFlow),
            ]
        )

# -- run deformation simulation
if defRun:
    baseCase.setParameters(
        [
            ['system/controlDict', 'endTime', str(flowEndTime + 1),''],
        ]
    )
    if nProc > 1:
        baseCase.updateTimesParallel()
        latestTimeOld = baseCase.latestParTime
        if slurmRun:
            pass
        else:
            baseCase.runCommands(
                [
                    'mkdir %g' % (baseCase.latestParTime),
                    'cp 0.org/D %g/'% (baseCase.latestParTime),
                    'decomposePar -fields > log.decomposeParFields',
                    'foamJob -parallel -screen %s > log.%s' % (solverDef, solverDef)
                ]
            )
        for procI in range(nProc):
            baseCase.updateTimesParallel()
            for fldToCopy in fldsToCopy:
                baseCase.runCommands(
                    [
                        'cp processor%d/%g/%s processor%d/%g' % (procI, latestTimeOld, fldToCopy, procI, baseCase.latestParTime)
                    ]
                )
    else:
        baseCase.updateTimes()
        latestTimeOld = baseCase.latestTime
        baseCase.runCommands(
            [
                'cp 0.org/D %g/'% (baseCase.latestTime),
                '%s > log.%s' % (solverDef, solverDef)
            ]
        )
        baseCase.updateTimes()
        for fldToCopy in fldsToCopy:
            baseCase.runCommands(
                [
                    'cp %g/%s %g' % (latestTimeOld, fldToCopy, baseCase.latestTime)
                ]
            )
