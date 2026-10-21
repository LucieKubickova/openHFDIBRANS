
from blockMeshDictClassV8 import *

def createMesh(baseCase):
    # parameters
    lPipe = 0.3
    wPipe = 0.05

    wSkip = 0.1
    wShift = wSkip + wPipe

    cStart = 0.09
    cEnd = 0.11
    cLen = cEnd - cStart

    dX = dY = dZ = 1.0e-3

    nCZ = 1

    x0 = y0 = z0 = 0.0
    grX = grY = grZ = "1.0"

    # create mesh
    fvMesh = mesh()

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, y0, z0
    xE, yE, zE = xC+cStart, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- top second block
    xC, yC, zC = xE, y0, z0
    xE, yE, zE = xC+cLen, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- top third block
    xC, yC, zC = xE, y0, z0
    xE, yE, zE = x0+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topThird = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom first block
    xC, yC, zC = x0, y0-wShift, z0
    xE, yE, zE = xC+cStart, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom second block
    xC, yC, zC = xE, y0-wShift, z0
    xE, yE, zE = xC+cLen, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom third block
    xC, yC, zC = xE, y0-wShift, z0
    xE, yE, zE = x0+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botThird = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- lambda block
    xC, yC, zC = x0+cStart, y0-wSkip, z0
    xE, yE, zE = xC+cLen, y0, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst, topSecond, topThird, botFirst, botSecond, botThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    lambdaBlock = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    ### PATCHES ###
    # -- empty
    empty = list()
    for block in fvMesh.blocks:
        empty.append(block.retFXY0())
        empty.append(block.retFXYE())

    fvMesh.addPatch("frontAndBack", "empty", empty)

    # -- inlets
    inletTop = list()
    inletTop.append(topFirst.retFYZ0())

    fvMesh.addPatch("inletTop", "patch", inletTop)

    inletBot = list()
    inletBot.append(botFirst.retFYZ0())

    fvMesh.addPatch("inletBot", "patch", inletBot)

    # -- outlets
    outletTop = list()
    outletTop.append(topThird.retFYZE())

    fvMesh.addPatch("outletTop", "patch", outletTop)

    outletBot = list()
    outletBot.append(botThird.retFYZE())

    fvMesh.addPatch("outletBot", "patch", outletBot)

    # -- walls
    walls = list()
    walls.append(topFirst.retFXZ0())
    walls.append(topFirst.retFXZE())

    walls.append(topSecond.retFXZE())

    walls.append(topThird.retFXZ0())
    walls.append(topThird.retFXZE())

    walls.append(botFirst.retFXZ0())
    walls.append(botFirst.retFXZE())

    walls.append(botSecond.retFXZ0())

    walls.append(botThird.retFXZ0())
    walls.append(botThird.retFXZE())

    walls.append(lambdaBlock.retFYZ0())
    walls.append(lambdaBlock.retFYZE())

    fvMesh.addPatch("walls", "wall", walls)

    ### WRITE ###
    # fvMesh.writeBMD("./system/")
    fvMesh.writeBMD("%s/system/" % baseCase.dir)

def createMeshLamBet(baseCase, dA = 1e-3):
    # parameters
    lPipe = 0.3
    wPipe = 0.05

    wSkip = 0.1
    wShift = wSkip + wPipe

    cStart = 0.09
    cEnd = 0.11
    cLen = cEnd - cStart

    dX = dY = dZ = dA

    nCZ = 1

    x0 = y0 = z0 = 0.0
    grX = grY = grZ = "1.0"

    # create mesh
    fvMesh = mesh()

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, y0, z0
    xE, yE, zE = xC+cStart, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- top second block
    xC, yC, zC = xE, y0, z0
    xE, yE, zE = xC+cLen, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- top third block
    xC, yC, zC = xE, y0, z0
    xE, yE, zE = x0+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topThird = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom first block
    xC, yC, zC = x0, y0-wShift, z0
    xE, yE, zE = xC+cStart, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom second block
    xC, yC, zC = xE, y0-wShift, z0
    xE, yE, zE = xC+cLen, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- bottom third block
    xC, yC, zC = xE, y0-wShift, z0
    xE, yE, zE = x0+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botThird = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- mid lambda block
    xC, yC, zC = x0+cStart, y0-wSkip, z0
    xE, yE, zE = xC+cLen, y0, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst, topSecond, topThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    midSolid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- left solid block
    xC, yC, zC = x0, y0-wSkip, z0
    xE, yE, zE = x0+cStart, y0, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [midSolid, topFirst, topSecond, botFirst, botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    leftSolid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    # -- right solid block
    xC, yC, zC = x0+cEnd, y0-wSkip, z0
    xE, yE, zE = x0+lPipe, y0, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [midSolid, topFirst, topSecond, botFirst, botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    rightSolid = fvMesh.addBlock(vertices, neighbours, nCells, grading)

    ### PATCHES ###
    # -- empty
    empty = list()
    for block in fvMesh.blocks:
        empty.append(block.retFXY0())
        empty.append(block.retFXYE())

    fvMesh.addPatch("frontAndBack", "empty", empty)

    # -- inlets
    inletTop = list()
    inletTop.append(topFirst.retFYZ0())

    fvMesh.addPatch("inletTop", "patch", inletTop)

    inletBot = list()
    inletBot.append(botFirst.retFYZ0())

    fvMesh.addPatch("inletBot", "patch", inletBot)

    # -- outlets
    outletTop = list()
    outletTop.append(topThird.retFYZE())

    fvMesh.addPatch("outletTop", "patch", outletTop)

    outletBot = list()
    outletBot.append(botThird.retFYZE())

    fvMesh.addPatch("outletBot", "patch", outletBot)

    # -- walls
    walls = list()
    walls.append(topFirst.retFXZE())
    walls.append(topSecond.retFXZE())
    walls.append(topThird.retFXZE())

    walls.append(botFirst.retFXZ0())
    walls.append(botSecond.retFXZ0())
    walls.append(botThird.retFXZ0())

    fvMesh.addPatch("walls", "wall", walls)

    # -- solid walls
    solidLeft = list()
    solidLeft.append(leftSolid.retFYZ0())

    fvMesh.addPatch("solidLeft", "wall", solidLeft)

    solidRight = list()
    solidRight.append(rightSolid.retFYZE())

    fvMesh.addPatch("solidRight", "wall", solidRight)

    fvMesh.writeBMD("%s/system/" % baseCase.dir)

def createMeshMulReg(baseCase):
    # parameters
    lPipe = 0.3
    wPipe = 0.05
    wCube = 0.015

    wSkip = 0.1
    wShift = wSkip + wCube

    cStart = 0.085
    cEnd = 0.115

    bStart = 0.09
    bEnd = 0.11

    cLen = cEnd - cStart

    dX = dY = dZ = 1.0e-3

    nCZ = 1

    x0 = y0 = z0 = 0.0
    grX = grY = grZ = "1.0"

    # create mesh
    fvMesh = mesh()

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, y0, z0
    xE, yE, zE = xC+cStart, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- top second block
    xC, yC, zC = x0+cStart, y0, z0
    xE, yE, zE = x0+bStart, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- top third block
    xC, yC, zC = x0+bStart, y0, z0
    xE, yE, zE = x0+bEnd, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topThird = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- top fourth block
    xC, yC, zC = x0+bEnd, y0, z0
    xE, yE, zE = x0+cEnd, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topFourth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- top fifth block
    xC, yC, zC = xE, y0, z0
    xE, yE, zE = x0+lPipe, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFourth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    topFifth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- above first block
    xC, yC, zC = x0, y0+wCube, z0
    xE, yE, zE = xC+cStart, y0+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topFirst, topSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    aboveFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- above second block
    xC, yC, zC = x0+cStart, y0+wCube, z0
    xE, yE, zE = x0+bStart, y0+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [aboveFirst, topFirst, topSecond, topThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    aboveSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- above third block
    xC, yC, zC = x0+bStart, y0+wCube, z0
    xE, yE, zE = x0+bEnd, y0+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [aboveSecond, topSecond, topThird, topFourth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    aboveThird = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- above fourth block
    xC, yC, zC = x0+bEnd, y0+wCube, z0
    xE, yE, zE = x0+cEnd, y0+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [aboveThird, topThird, topFourth, topFifth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    aboveFourth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- above fifth block
    xC, yC, zC = xE, y0+wCube, z0
    xE, yE, zE = x0+lPipe, y0+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [aboveFourth, topFourth, topFifth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    aboveFifth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    # -- bottom first block
    xC, yC, zC = x0, y0-wShift, z0
    xE, yE, zE = xC+cStart, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- bottom second block
    xC, yC, zC = x0+cStart, y0-wShift, z0
    xE, yE, zE = x0+bStart, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botFirst]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- bottom third block
    xC, yC, zC = x0+bStart, y0-wShift, z0
    xE, yE, zE = x0+bEnd, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botThird = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- bottom fourth block
    xC, yC, zC = x0+bEnd, y0-wShift, z0
    xE, yE, zE = x0+cEnd, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botFourth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    # -- bottom fifth block
    xC, yC, zC = xE, y0-wShift, z0
    xE, yE, zE = x0+lPipe, yC+wCube, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botFourth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    botFifth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- below first block
    xC, yC, zC = x0, y0-wShift+wCube-wPipe, z0
    xE, yE, zE = xC+cStart, y0-wShift, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [botFirst, botSecond]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    belFirst = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- below second block
    xC, yC, zC = x0+cStart, y0-wShift+wCube-wPipe, z0
    xE, yE, zE = x0+bStart, y0-wShift, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [belFirst, botFirst, botSecond, botThird]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    belSecond = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- below third block
    xC, yC, zC = x0+bStart, y0-wShift+wCube-wPipe, z0
    xE, yE, zE = x0+bEnd, y0-wShift, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [belSecond, botSecond, botThird, botFourth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    belThird = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- below fourth block
    xC, yC, zC = x0+bEnd, y0-wShift+wCube-wPipe, z0
    xE, yE, zE = x0+cEnd, y0-wShift, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [belThird, botThird, botFourth, botFifth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    belFourth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- below fifth block
    xC, yC, zC = xE, y0-wShift+wCube-wPipe, z0
    xE, yE, zE = x0+lPipe, y0-wShift, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [belFourth, botFourth, botFifth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    belFifth = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    # -- lambda block
    xC, yC, zC = x0+bStart, y0-wSkip, z0
    xE, yE, zE = x0+bEnd, y0, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = [topSecond, topThird, topFourth, botSecond, botThird, botFourth]

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    lambdaBlock = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    ### PATCHES ###
    # -- empty
    empty = list()
    for block in fvMesh.blocks:
        empty.append(block.retFXY0())
        empty.append(block.retFXYE())

    fvMesh.addPatch("frontAndBack", "empty", empty)

    # -- inlets
    inletTop = list()
    inletTop.append(topFirst.retFYZ0())
    inletTop.append(aboveFirst.retFYZ0())

    fvMesh.addPatch("inletTop", "patch", inletTop)

    inletBot = list()
    inletBot.append(botFirst.retFYZ0())
    inletBot.append(belFirst.retFYZ0())

    fvMesh.addPatch("inletBot", "patch", inletBot)

    # -- outlets
    outletTop = list()
    outletTop.append(topFifth.retFYZE())
    outletTop.append(aboveFifth.retFYZE())

    fvMesh.addPatch("outletTop", "patch", outletTop)

    outletBot = list()
    outletBot.append(botFifth.retFYZE())
    outletBot.append(belFifth.retFYZE())

    fvMesh.addPatch("outletBot", "patch", outletBot)

    # -- walls
    walls = list()
    walls.append(topFirst.retFXZ0())
    walls.append(topSecond.retFXZ0())
    walls.append(topFourth.retFXZ0())
    walls.append(topFifth.retFXZ0())

    walls.append(aboveFirst.retFXZE())
    walls.append(aboveSecond.retFXZE())
    walls.append(aboveThird.retFXZE())
    walls.append(aboveFourth.retFXZE())
    walls.append(aboveFifth.retFXZE())

    walls.append(botFirst.retFXZE())
    walls.append(botSecond.retFXZE())
    walls.append(botFourth.retFXZE())
    walls.append(botFifth.retFXZE())

    walls.append(belFirst.retFXZ0())
    walls.append(belSecond.retFXZ0())
    walls.append(belThird.retFXZ0())
    walls.append(belFourth.retFXZ0())
    walls.append(belFifth.retFXZ0())

    walls.append(lambdaBlock.retFYZ0())
    walls.append(lambdaBlock.retFYZE())

    fvMesh.addPatch("walls", "wall", walls)

    ### WRITE ###
    fvMesh.writeBMD("%s/system/" % baseCase.dir)

def createMeshMulRegOnlyChannels(baseCase):
    # parameters
    lPipe = 0.3
    wPipe = 0.05
    wCube = 0.015
    wSkip = 0.1

    dX = dY = dZ = 1.0e-3

    nCZ = 1

    x0 = y0 = z0 = 0.0
    grX = grY = grZ = "1.0"

    # create mesh
    fvMesh = mesh()

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, y0, z0
    xE, yE, zE = xC+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    fluidCold = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidCold")

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, yC+wPipe, z0
    xE, yE, zE = xC+lPipe, yC+wSkip, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    solid = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "solid")

    ### BLOCKS ###
    # -- top first block
    xC, yC, zC = x0, yC+wSkip, z0
    xE, yE, zE = xC+lPipe, yC+wPipe, zC+dZ

    # vertices
    vertices = [
            [xC, yC, zC],
            [xE, yC, zC],
            [xE, yE, zC],
            [xC, yE, zC],
            [xC, yC, zE],
            [xE, yC, zE],
            [xE, yE, zE],
            [xC, yE, zE],
        ]

    # neighbouring blocks
    neighbours = []

    # number of cells
    nCX = int(round(abs(xE-xC)/dX))
    nCY = int(round(abs(yE-yC)/dY))
    nCells = [nCX, nCY, nCZ]

    # grading
    grading = [grX, grY, grZ]

    # create the block
    fluidHot = fvMesh.addBlock(vertices, neighbours, nCells, grading, name = "fluidHot")

    ### PATCHES ###
    # -- empty
    empty = list()
    for block in fvMesh.blocks:
        empty.append(block.retFXY0())
        empty.append(block.retFXYE())

    fvMesh.addPatch("frontAndBack", "empty", empty)

    # -- inlets
    inletTop = list()
    inletTop.append(fluidHot.retFYZ0())
    # inletTop.append(aboveFirst.retFYZ0())

    fvMesh.addPatch("inletTop", "patch", inletTop)

    inletBot = list()
    inletBot.append(fluidCold.retFYZ0())
    # inletBot.append(belFirst.retFYZ0())

    fvMesh.addPatch("inletBot", "patch", inletBot)

    # -- outlets
    outletTop = list()
    outletTop.append(fluidHot.retFYZE())
    # outletTop.append(aboveFifth.retFYZE())

    fvMesh.addPatch("outletTop", "patch", outletTop)

    outletBot = list()
    outletBot.append(fluidCold.retFYZE())
    # outletBot.append(belFifth.retFYZE())

    fvMesh.addPatch("outletBot", "patch", outletBot)

    # -- walls
    walls = list()
    walls.append(fluidCold.retFXZ0())
    # walls.append(topSecond.retFXZ0())
    # walls.append(topFourth.retFXZ0())
    # walls.append(topFifth.retFXZ0())

    walls.append(fluidHot.retFXZE())

    # -- solid walls
    solidLeft = list()
    solidLeft.append(solid.retFYZ0())

    fvMesh.addPatch("solidLeft", "wall", solidLeft)

    solidRight = list()
    solidRight.append(solid.retFYZE())

    fvMesh.addPatch("solidRight", "wall", solidRight)
    # walls.append(aboveSecond.retFXZE())
    # walls.append(aboveThird.retFXZE())
    # walls.append(aboveFourth.retFXZE())
    # walls.append(aboveFifth.retFXZE())

    # walls.append(botFirst.retFXZE())
    # walls.append(botSecond.retFXZE())
    # walls.append(botFourth.retFXZE())
    # walls.append(botFifth.retFXZE())

    # walls.append(belFirst.retFXZ0())
    # walls.append(belSecond.retFXZ0())
    # walls.append(belThird.retFXZ0())
    # walls.append(belFourth.retFXZ0())
    # walls.append(belFifth.retFXZ0())

    # walls.append(lambdaBlock.retFYZ0())
    # walls.append(lambdaBlock.retFYZE())

    fvMesh.addPatch("walls", "wall", walls)

    ### WRITE ###
    fvMesh.writeBMD("%s/system/" % baseCase.dir)