To install the OpenHFDIB-RANS library, you need to meet the prerequisites. Then, the installation consists of two steps: cloning and compiling the library.

## Prerequisites
* OpenFOAM 8 built or installed according to the OpenFOAM Foundation [instructions](https://openfoam.org/download/8-source/)
* (Optional) ParaView to visualize OpenFOAM fields

## Cloning the repository
First, source your OpenFOAM 8 installation and verify it:
```bash
source path/to/OpenFOAM-8/etc/bashrc
echo "$WM_PROJECT $WM_PROJECT_VERSION"
```

Navigate to your folder of choice and run:
```bash
git clone https://github.com/LucieKubickova/openHFDIBRANS.git
```
The library source code should now be located in the `./openHFDIBRANS` folder.

## Compiling the library
To compile the library, ensure your OpenFOAM 8 installation is still sourced:
```bash
echo "$WM_PROJECT $WM_PROJECT_VERSION"
```
Proceed into the source code folder:
```bash
cd openHFDIBRANS
```
Compile the library and solvers by running:
```bash
./Allwclean && ./Allwmake
```

If you run into any problems, feel free to open an [issue](https://github.com/LucieKubickova/openHFDIBRANS/issues).
