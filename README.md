# openHFDIBRANS
Hybrid fictitious domain-immersed boundary (HFDIB) method extended for Reynolds-averaged simulation (RAS) with wall functions at the immersed boundaries. The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB), but the code was heavily modified.

## Initial verification result on a bent pipe
![introImage](https://user-images.githubusercontent.com/114754867/223365550-765b8994-b5b6-45f7-8797-ac02a854197c.png)

## Cite this work as
L. Kubíčková and M. Isoz.: On Reynolds-Averaged Turbulence Modeling with Immersed Boundary Method. In Proceedings of Topical Problems of Fluid Mechanics 2023, Prague, 2023, Edited by David Šimurda and Tomáš Bodnár, pp. 104–111., DOI: https://doi.org/10.14311/TPFM.2020.020 

## Compatibility
The code is prepared for compilation with OpenFOAMv8 (https://openfoam.org/version/8/).

## Compilation
./Allwclean && ./Allwmake

## Tutorials
cd tutorials && ./Allrun
Note: the results for simpleFoam and HFDIB-RAS (simpleHFDIBFoam) will be stored in respective folders
