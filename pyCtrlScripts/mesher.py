
def mesh():
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

    fvMesh.addPatch("solidLeft", "patch", solidLeft)

    solidRight = list()
    solidRight.append(rightSolid.retFYZE())

    fvMesh.addPatch("solidRight", "patch", solidRight)

    ### WRITE ###
    fvMesh.writeBMD("./system/")
