
import os

# pars
regions = ["fluidCold", "fluidHot", "solid"]

toAdd = [
    ["U", "T", "omega", "nut", "k", "alphat", "p_rgh"],
    ["U", "T", "omega", "nut", "k", "alphat", "p_rgh"],
    ["T"],
    ]

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

# loop over regions
for r in range(len(regions)):
    region = regions[r]
    changeFile = open("./system/" + region + "/changeDictionaryDict", 'w')

    changeFile.write(dictStart)

    for a in range(len(toAdd[r])):
        changeFile.write(toAdd[r][a] + "\n{\n")
        with open("./0.org/" + region + "/" + toAdd[r][a], 'r') as file:
            lines = file.readlines()

        for l in range(len(lines)):
            if "internalField" in lines[l]:
                lS = l
                break

        changeFile.writelines(lines[lS:-1])
        changeFile.write("}\n\n")

    changeFile.write(dictEnd)
