# openHFDIBRANS
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22917049-blue.svg)](https://doi.org/10.5281/zenodo.22917049)

openHFDIBRANS is an open-source library for OpenFOAM that features the extension of the hybrid fictitious domain-immersed boundary (HFDIB) method for steady-state Reynolds-averaged simulation (RAS). The library includes custom implementation of wall functions at the immersed boundary, modified turbulence models, and a steady-state solver.

The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB), but the code was heavily modified. Its variant for CFD-DEM simulations can be found in (https://github.com/techMathGroup/openHFDIB-DEM).

## Solver results on the backward facing step benchmark
<p align="center">
  <img src="https://github.com/LucieKubickova/openHFDIBRANS/blob/main/Images/backwardFacingStepBenchmark.png">
</p>

## Compatibility
The code is prepared for compilation with OpenFOAM v8 (https://openfoam.org/version/8/).

## For users
For information regarding the compilation and usage of this solver, please refer to the [wiki](https://github.com/LucieKubickova/openHFDIBRANS/wiki).

## Cite this work as
* L. Kubíčková and M. Isoz.: Extending the hybrid fictitious domain-immersed boundary method for reynolds-averaged turbulence modeling, 2026. URL: https://arxiv.org/abs/2606.06135. arXiv:2606.06135
* L. Kubíčková and M. Isoz.: On Reynolds-Averaged Turbulence Modeling with Immersed Boundary Method. In Proceedings of Topical Problems of Fluid Mechanics 2023, Prague, 2023, Edited by David Šimurda and Tomáš Bodnár, pp. 104–111., DOI: https://doi.org/10.14311/TPFM.2023.015

## License
openHFDIBRANS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.
