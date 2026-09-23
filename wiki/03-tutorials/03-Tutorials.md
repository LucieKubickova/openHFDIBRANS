Here, we showcase the features of the library on select tutorial cases included in the project's repository.

| Tutorial | Location | Description | Link |
| :------- | :------- | :---------- | :--- |
| Backward-facing step | tutorials/backwardFacingStep | 2D verification case, lambda field generation, turbulence models, wall shear stress | [View](#2D-flow-over-a-backward-facing-step) |

## 2D flow over a backward-facing step
### Case description
![The backward-facing step domain diagram with measurements and boundaries](03-tutorials/bfs.png)

This tutorial introduces the case setup for a 2D separated flow over a backward-facing step. The diagram above details the domain measurements and patch labels.

For detailed domain settings see `system/blockMeshDict`. For boundary conditions of every field see the `0.org/` folder and files therein.

### Prerequisites
To run this tutorial case, you need to have the OpenHFDIBRANS library installed and sourced. Verify this with:
```bash
source path/to/OpenFOAM-8/etc/bashrc
echo "$WM_PROJECT $WM_PROJECT_VERSION"
```

Next, copy the tutorial from the library source to a desired destination:
```bash
cp -r path/to/openHFDIBRANS/tutorials/backwardFacingStep path/to/my/directory
cd path/to/my/directory/backwardFacingStep
```

### Lambda field setup
The lambda field defines the presence of the immersed boundary in the computational domain. The field can be generated from an STL surface geometry.

To do so, place the desired STL file into `constant/triSurface`, define the `stlName` in `constant/HFDIBDEMDict` and run the lambda generation utility:
```bash
generateLambda
```
This utility generates a `lambda` field file in `0.org` from the provided STL file. Alternatively, set *generateLambda* to *true* in `constant/HFDIBDEMDict` to have the solver generate the lambda field automatically at startup.

### Important inputs
Below is the summary of important dictionary entries introduced by OpenHFDIB-RANS. For detailed explanations of individual parameters, please refer to the [user guide](#02-User-guide).

| Input | Description |
| :---- | :---------- |
| `constant/HFDIBDEMDict` | Configures the general OpenHFDIB-RANS behaviour. |
| `constant/turbulenceProperties` | Configures the flow regime and turbulence model. |
| `system/fvSchemes` | Configures the discretization schemes. Mind the terms in *gradSchemes* and *divSchemes*. |
| `system/fvSolution` | Configures the solver settings. |
| `system/controlDict` | General simulation settings. Set the application to *simpleHFDIBRANSFoam*. |
| `0.org/` | Configure the initial and boundary conditions for various fields. Add a `lambda` field file here to represent the immersed body. |

### Running the case
Run the tutorial case by calling:
```bash
./Allrun
```

The `Allrun` script contains the following:
```bash
#!/bin/bash

# Set OpenFOAM environment
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Clear old initial conditions
rm -rf 0

# Create base mesh
runApplication blockMesh

# Copy initial conditions
mkdir -p 0
cp -rf 0.org/* 0

# Create files for postprocessing
paraFoam -touch

# Get the solver application
application=`getApplication`

# Run the application
runApplication $application

# Calculate the wall shear stress
runApplication calculateWallShearStress
```

Once the simulation finishes, the resulting fields, along with the `tauw` field produced by `calculateWallShearStress`, are written to the latest time directory. Open the case in ParaView to further inspect the results.
