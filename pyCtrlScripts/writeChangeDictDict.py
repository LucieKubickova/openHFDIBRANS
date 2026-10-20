
import os

# pars

# dict parts
dictStart = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

dictEnd = """
// ************************************************************************* //
"""

def writeChangeDict(baseCase, regions, fldInReg):
    # loop over regions
    for r in range(len(regions)):
        region = regions[r]
        changeFile = open("%s/system/"%(baseCase.dir) + region + "/changeDictionaryDict", 'w')

        changeFile.write(dictStart)

        for a in range(len(fldInReg[r])):
            changeFile.write(fldInReg[r][a] + "\n{\n")
            with open("%s/0.org/" %(baseCase.dir) + region + "/" + fldInReg[r][a], 'r') as file:
                lines = file.readlines()

            for l in range(len(lines)):
                if "internalField" in lines[l]:
                    lS = l
                    break

            changeFile.writelines(lines[lS:-1])
            changeFile.write("}\n\n")

        changeFile.write(dictEnd)
