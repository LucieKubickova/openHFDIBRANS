
import csv
import math
import random
import numpy as np
import matplotlib.pyplot as plt

field = 'k'
pos = 0
time = 1000

fileName = "./ZZ_python/" + str(time) + "/interpolation_" + field + ".dat"
with open(fileName, 'r') as file:
    reader = csv.reader(file)

    cols = next(reader)

    data = list()
    for line in reader:
        data.append(line)

ids = cols.index("ds")
iphi = cols.index("phi")
icfs = cols.index("coeffs")

n = 9
total = len(data)

lines = [random.randrange(total-1) for i in range(n)]

fig = plt.figure(figsize = (20,12))

k = 1
for i in lines:
    ax = fig.add_subplot(int(math.sqrt(n)),int(math.sqrt(n)),k)

    line = data[i]

    x = [float(i) for i in line[1][2:-1].split()]

    if field == "grad(U)" or field == "U":
        y = [float(i.split()[pos]) for i in line[2][3:-2].split(") (")]
        coeffs = [float(i.split()[pos]) for i in line[3][3:-2].split(") (")]

    else:
        y = [float(i) for i in line[2][2:-1].split()]
        coeffs = [float(i) for i in line[3][2:-1].split()]

    ax.scatter(x[0], y[0], color = 'black', marker = 'x')
    ax.scatter(x, y)

    xx = np.linspace(min(x), max(x), 100)

    if len(coeffs) == 2:
        if field == "grad(U)":
            yy = [coeffs[0]*i**2 + coeffs[1]*i + y[1] for i in xx]

        else:
            yy = [coeffs[0]*i**2 + coeffs[1]*i + y[0] for i in xx]

        ax.plot(xx, yy, color = 'black')

    elif len(coeffs) == 1:
        yy = [coeffs[0]*i + y[0] for i in xx]
        ax.plot(xx, yy, color = 'black')

    k += 1

plt.savefig("./" + field + "_" + str(time) + ".png")
plt.show()
