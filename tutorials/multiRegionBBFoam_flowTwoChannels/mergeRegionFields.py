
import os

# pars
regions = ["fluidCold", "solid", "fluidHot"]
mapDict = "./0/cellToRegion"
fields = ["T", "U", "kappa"]
types = ["Scalar", "Vector", "Scalar"]
dimensions = ["[0 0 0 1 0 0 0]", "[0 1 -1 0 0 0 0]", "[1 1 -3 -1 0 0 0]"]

# read end time
with open("./system/controlDict", 'r') as file:
    lines = file.readlines()

for line in lines:
    if line[:7] == "endTime":
        time = line.split()[-1][:-1]
        break

timeDir = "./" + time + "/"

dictStart = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2412                                  |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
"""

dictMid = """
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

dictEnd = """
boundaryField
{
    inletTop
    {
        type            zeroGradient;
    }

    inletBot
    {
        type            zeroGradient;
    }

    "(outletTop|outletBot)"
    {
        type            zeroGradient;
    }

    walls
    {
        type            zeroGradient;
    }

    wallInsideLambda
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
"""

# read mapping
with open(mapDict, 'r') as file:
    lines = file.readlines()

for l in range(len(lines)):
    line = lines[l]
    if "internalField" in line:
        lS = l+3
        break

for l in range(len(lines)):
    line = lines[l]
    if "boundaryField" in line:
        lE = l-3
        break

cellMap = list()
for line in lines[lS:lE]:
    cellMap.append(int(line))

# loop over fields
for f in range(len(fields)):
    # read data
    data = list()
    for region in regions:
        fileName = timeDir + region + "/" + fields[f]
        with open(fileName, 'r') as file:
            lines = file.readlines()
    
        uniform = False
        for l in range(len(lines)):
            line = lines[l]
            if "internalField" in line:
                if " uniform" in line:
                    uniform = True
                    uniValue = line.split(" uniform ")[-1][:-2] + "\n"

                lS = l+3
                break
    
        for l in range(len(lines)):
            line = lines[l]
            if "boundaryField" in line:
                lE = l-3
                break
    
        regData = list()
        for line in lines[lS:lE]:
            regData.append(line)
    
        if uniform:
            count = cellMap.count(regions.index(region))
            regData = [uniValue]*count

        data.append(regData)
    
    # write merged file
    inds = [0]*len(regions)
    
    with open(timeDir + fields[f], 'w') as file:
        file.write(dictStart)

        file.write("    class       vol" + types[f] + "Field;\n")
        file.write("    location    " + time + ";\n")
        file.write("    object      " + fields[f] + ";\n")

        file.write(dictMid)
        file.write("dimensions      " + dimensions[f] + ";\n\n")
        file.write("internalField   nonuniform List<" + types[f].lower() + ">\n")
        file.write(str(len(cellMap)) + "\n(\n")
    
        for r in cellMap:
            value = data[r][inds[r]]
            inds[r] += 1
    
            file.write(value)
    
        file.write(")\n;\n\n")
        file.write(dictEnd)
