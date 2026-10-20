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
from writeChangeDictDict import *

# CASE FOLDERS==========================================================
baseCaseDir = '../tutorials/multiRegionBBFoam_flowTwoChannels/' # -- base case for simulation
outFolder = '../ZZ_cases/94_V5_multiRegion/'
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

# -- regions
regions = ["fluidCold", "fluidHot", "solid"]
fldInReg = [
    ["U", "T", "omega", "nut", "k", "alphat", "p_rgh"], # fluidCol
    ["U", "T", "omega", "nut", "k", "alphat", "p_rgh"], # fluidHot
    ["T"],  # solid
]

# -- decomposition
nProc = 1

# -- fields to copy for postProcessing
fldsToCopy = ["T", "U"]

# -- solvers
solverFlow = 'chtMultiRegionBuoyantBoussinesqSimpleFoam' 
solverDef = 'solids4CustomMultiRegionFoam' 

# -- numerics
nDefTIters = 10000  # -- number of solid/temperature iterations at each time level
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
values = [k0, o0, n0]
names = ["k", "omega", "nut"]

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
    baseCase.setParameters(
        [
            ['0.org/fluidCold/%s'%(name), 'internalField', 'uniform %g' % (values[names.index(name)]),'']
        ]
    )
    baseCase.setParameters(
        [
            ['0.org/fluidHot/%s'%(name), 'internalField', 'uniform %g' % (values[names.index(name)]),'']
        ]
    )
# -- velocity
baseCase.setParameters(
    [
        ['0.org/U', 'internalField', 'uniform (%s 0 0)'%(uIn), ''],
        ['0.org/fluidCold/U', 'internalField', 'uniform (%s 0 0)'%(uIn), ''],
        ['0.org/fluidHot/U', 'internalField', 'uniform (%s 0 0)'%(uIn), ''],
    ]
)

# 2) constant/solidProperties
baseCase.setParameters(
    [
        ["constant/solid/solidProperties", 'nCorrectors', str(nDefTIters), ""],
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

# -- create mesh 
if crtMesh:
    createMeshMulRegOnlyChannels(baseCase)
    baseCase.runCommands(
        [
            'blockMesh > log.blockMesh',
            'mkdir -p 0',
            'cp -rf 0.org/* 0/',
            'paraFoam -touch',
        ]
    )
    writeChangeDict(baseCase, regions, fldInReg)
    baseCase.runCommands(
        [
            'splitMeshRegions -cellZones -overwrite > log.splitMeshRegions',
        ]
    )
    for region in regions:
        baseCase.runCommands(
            [
                'changeDictionary -region %s > log.changeDictionary' %(region),
                # './Allrun-serialAll',
                # 'runApplication -s %s changeDictionary -region %s' %(region, region)
            ]
        )
    baseCase.runCommands(
        [
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
                # 'chmod 755 ./Allclean && chmod 755 ./Allrun-serialAll',
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
        if slurmRun:
            pass
        else:
            baseCase.runCommands(
                [
                    'mkdir %g' % (baseCase.latestParTime),
                    'cp 0.org/D %g/'% (baseCase.latestParTime),
                    'decomposePar > log.decomposePar',
                    'foamJob -parallel -screen %s > log.%s' % (solverDef, solverDef)
                ]
            )
    else:
        baseCase.updateTimes()
        latestTimeOld = baseCase.latestTime
        baseCase.runCommands(
            [
                'cp 0.org/solid/D %g/solid/'% (baseCase.latestTime),
                '%s > log.%s' % (solverDef, solverDef)
            ]
        )
        baseCase.updateTimes()
        for fldToCopy in fldsToCopy:
            baseCase.runCommands(
                [
                    'cp %g/fluidCold/%s %g' % (latestTimeOld, fldToCopy, baseCase.latestTime),
                    'cp %g/fluidHot/%s %g' % (latestTimeOld, fldToCopy, baseCase.latestTime)
                ]
            )
        baseCase.runCommands(
            [
                'cp %g/solid/T %g' % (latestTimeOld, baseCase.latestTime),
                'cp %g/solid/T %g' % (latestTimeOld, baseCase.latestTime)
            ]
        )