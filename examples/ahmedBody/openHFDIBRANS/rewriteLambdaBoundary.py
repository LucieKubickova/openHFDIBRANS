
fileName = './constant/polyMesh/boundary'
with open(fileName, 'r') as file:
    data = file.readlines()

patches = list()
for i in range(len(data)):
    if 'type' in data[i]:
        name = data[i-2]
        if 'wedge' in data[i] or 'empty' in data[i]:
            info = name + data[i-1] + data[i]

        else:
            info = name + data[i-1] + "\t\ttype\t\tcalculated;\n"

        patches.append(info)

fileName = './0/lambda'
with open(fileName, 'r') as file:
    data = file.readlines()

for i in range(len(data)):
    if 'boundaryField' in data[i]:
        bIndex = i
        break

del data[bIndex:]

with open(fileName, 'w') as file:
    for line in data:
        file.write(line)

    file.write("boundaryField\n{\n")

    for patch in patches:
        file.write(patch)
        file.write("\t\tvalue\t\tuniform 0;\n")
        file.write("\t}\n")

    file.write("}\n\n\n")
    file.write("// ************************************************************************* //")
