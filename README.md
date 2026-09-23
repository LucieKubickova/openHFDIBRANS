# openHFDIBRANS
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22917049-blue.svg)](https://doi.org/10.5281/zenodo.22917049)

Hybrid fictitious domain-immersed boundary (HFDIB) method extended for Reynolds-averaged simulation (RAS) with wall functions at the immersed boundaries. The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB), but the code was heavily modified.

## Solver results on the backward facing step benchmark
<p align="center">
  <img src="https://github.com/LucieKubickova/openHFDIBRANS/blob/main/Images/backwardFacingStepBenchmark.png">
</p>

## Cite this work as
article in preparation

## Compatibility
The code is prepared for compilation with OpenFOAM v8 (https://openfoam.org/version/8/).

## Compilation
./Allwclean && ./Allwmake

## Tutorials
cd tutorials && ./Allrun

## License
openHFDIBRANS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.
