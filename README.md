# openHFDIBRANS
openHFDIBRANS is an open-source CFD library based on OpenFOAM (https://openfoam.org) that is linked to openHFDIB-DEM (https://github.com/MartinKotoucSourek/openHFDIB-DEM) and it is prepared to run Reynolds-averaged simulations (RAS) with Hybrid fictitious domain-immersed boundary method (HFDIB). The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB), but the code was heavily modified and the connection to RAS is new.

![introImage](https://user-images.githubusercontent.com/114754867/223365550-765b8994-b5b6-45f7-8797-ac02a854197c.png)

## Cite this work as
L. Kubíčková and M. Isoz.: On Reynolds-Averaged Turbulence Modeling with Immersed Boundary Method. In Proceedings of Topical Problems of Fluid Mechanics 2023, Prague, 2023, Edited by David Šimurda and Tomáš Bodnár, pp. 104–111., DOI: 10.14311/TPFM.2023.015

## Compatibility
The code is prepared for compilation with OpenFOAMv8 (https://openfoam.org/version/8/).

## Compilation
./Allwclean && ./Allwmake

## Tutorials
cd tutorials && ./Allrun
