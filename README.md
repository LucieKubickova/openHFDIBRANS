# openHFDIBRANS with cut cells
Hybrid fictitious domain-immersed boundary (HFDIB) method extended for Reynolds-averaged simulation (RAS) with wall functions at the immersed boundaries. The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB), but the code was heavily modified.

## Initial verification result on a bent pipe
<p align="center">
  <img src="https://github.com/LucieKubickova/openHFDIBRANS/assets/114754867/a76d5efd-734c-4487-ba38-6e41232f7ddd.png">
</p>

## Cite this work as
L. Kubíčková and M. Isoz.: On Reynolds-Averaged Turbulence Modeling with Immersed Boundary Method. In Proceedings of Topical Problems of Fluid Mechanics 2023, Prague, 2023, Edited by David Šimurda and Tomáš Bodnár, pp. 104–111., DOI:  https://doi.org/10.14311/TPFM.2023.015

## Compatibility
The code is prepared for compilation with OpenFOAM v8 (https://openfoam.org/version/8/).

## Compilation
./Allwclean && ./Allwmake

## Tutorials
cd tutorials && ./Allrun
Note: the results for simpleFoam and HFDIB-RAS (simpleHFDIBRANSFoam) will be stored in respective folders
