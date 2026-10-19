
import os

# functions
def writeValue(value, fileName):
    with open(fileName, 'r') as file:
        lines = file.readlines()

    for l in range(len(lines)):
        line = lines[l]
        if "internalField" in line:
            newLine = line.split("uniform")[0]
            newLine += "uniform " + str(value) + ";\n"
            break

    with open(fileName, 'w') as file:
        file.writelines(lines[:l])
        file.write(newLine)
        file.writelines(lines[l+1:])

def writeVelocity(value, fileName):
    with open(fileName, 'r') as file:
        lines = file.readlines()

    for l in range(len(lines)):
        line = lines[l]
        if "internalField" in line:
            newLine = line.split("uniform")[0]
            newLine += "uniform (" + str(value) + " 0 0);\n"
            break

    with open(fileName, 'w') as file:
        file.writelines(lines[:l])
        file.write(newLine)
        file.writelines(lines[l+1:])

# pars
lScale = 0.005
uIn = 0.001
Iper = 2.0
Iper /= 100
Cmu = 0.09

names = ["k", "epsilon", "omega", "nut"]

# calculate initial values
k0 = 3/2*uIn**2*Iper**2
e0 = Cmu*k0**(3/2)/lScale
o0 = e0/Cmu/k0
n0 = k0/o0

values = [k0, e0, o0, n0]

# read 0.org
fileDir = "./0.org/"
ls0 = os.listdir(fileDir)

# write to 0.org
for name in names:
    if name in ls0:
        writeValue(values[names.index(name)], fileDir + ls0[ls0.index(name)])

if "U" in ls0:
    writeVelocity(uIn, fileDir + ls0[ls0.index("U")])

# read col region
fileDir = "./0.org/fluidCold/"
ls0 = os.listdir(fileDir)

# write to 0.org
for name in names:
    if name in ls0:
        writeValue(values[names.index(name)], fileDir + ls0[ls0.index(name)])

if "U" in ls0:
    writeVelocity(uIn, fileDir + ls0[ls0.index("U")])

# read hot region
fileDir = "./0.org/fluidHot/"
ls0 = os.listdir(fileDir)

# write to 0.org
for name in names:
    if name in ls0:
        writeValue(values[names.index(name)], fileDir + ls0[ls0.index(name)])

if "U" in ls0:
    writeVelocity(uIn, fileDir + ls0[ls0.index("U")])
