
import csv
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

# geometric parameters
grid = False

# show box
xMin = 0.0
xMax = 0.1
yMin = 0.0
yMax = 0.015

# read data
fileName = "ZZ_python/interpolationInfo.dat"
with open(fileName, 'r') as file:
    reader = csv.reader(file)

    cols = next(reader)

    data = list()
    for line in reader:
        data.append(line)

# find indexes of relevant data
iCC = cols.index("cellCenter")
iIP = cols.index("intPoints")

# save relevant data
xCC = list()
yCC = list()
xIP = list()
yIP = list()
xSF = list()
ySF = list()

for line in data:
    CC = [float(i) for i in line[iCC][1:-1].split()]

    if CC[0] > xMin and CC[0] < xMax and CC[1] > yMin and CC[1] < yMax:
        xCC.append(CC[0])
        yCC.append(CC[1])

        IP = line[iIP][3:-2].split(") (")
        xIPs = [float(ip.split()[0]) for ip in IP]
        yIPs = [float(ip.split()[1]) for ip in IP]

        xSF.append(float(IP[0].split()[0]))
        ySF.append(float(IP[0].split()[1]))

        xIP.append(xIPs)
        yIP.append(yIPs)

# prepare things for plots
if grid:
    vs = np.linspace(LMxT, LMxT+LDiff, nCX+1)
    hs = np.linspace(WMxT, WDiff, nCY+1)
    
    vlines = list()
    hlines = list()
    
    for i in range(len(vs)):
        if vs[i] > xMin and vs[i] < xMax:
            vlines.append(vs[i])
    
    for i in range(len(hs)):
        if hs[i] > yMin and hs[i] < yMax:
            hlines.append(hs[i])

nCols = 10
colors = cm.rainbow(np.linspace(0, 1, nCols))

# plot data
fig = plt.figure(figsize = (20,12))
ax = fig.add_subplot(111)

# plot cell centers and interpolation points
for i in range(len(xCC)):
    x = [xCC[i]]
    y = [yCC[i]]

    for j in range(len(xIP[i])):
        x.append(xIP[i][j])
        y.append(yIP[i][j])

    ax.plot(x, y, color = colors[i%nCols])
    ax.scatter(x, y, color = colors[i%nCols])
    ax.scatter(x[0], y[0], color = "black", marker = "x")

# plot surface line
ax.plot(xSF, ySF, color = "black")

# plot help grid
if grid:
    ax.hlines(hlines, min(vlines), max(vlines), color = "grey")
    ax.vlines(vlines, min(hlines), max(hlines), color = "grey")

plt.savefig("./interpolationPoints.png")
plt.show()
