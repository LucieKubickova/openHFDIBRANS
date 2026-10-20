Here, the individual solver options are introduced and described. Below is the reference index for all custom dictionaries, subdictionaries, and utilities introduced by this library.

| Dictionary | Location | Description | Link |
| :--------- | :------- | :---------- | :--- |
| HFDIBDEMDict | `constant/` | Configures the interpolation and immersed boundary settings.  | [View](#HFDIBDEMDict) |
| HFDIBRAS subdictionary | `constant/turbulenceProperties` | Configures the flow turbulence models used. | [View](#HFDIBRAS-in-turbulenceProperties)  |
| HFDIBSchemes subdictionary | `system/fvSchemes` | Configures the discretization schemes used. | [View](#HFDIBSchemes-in-fvSchemes) |
| HFDIB subsubdictionary | `system/fvSolution` | Configures the SIMPLE algorithm used by the solver.  | [View](#HFDIB-in-SIMPLE-in-fvSolution) |

| Utility | Description | Link |
| :------ | :---------- | :--- |
| calculateWallShearStress | Calculates the wall shear stress exerted on the immersed boundary. | [View](#calculateWallShearStress) |
| calculateForceCoeffsLaminar | Calculates the aerodynamic force coefficients acting on the immersed boundary for laminar flow. | [View](#calculateForceCoeffs)  |
| calculateForceCoeffs | Calculates the aerodynamic force coefficients acting on the immersed boundary for turbulent flow. | [View](#calculateForceCoeffs)  |

## HFDIBDEMDict

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      HFDIBDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceThreshold    1e-4;
interfaceSpan       1.0;
sdBasedLambda       false;
stlName             cylinder;
saveIntInfo         false;
surfAreaType        cutEdges;

wallFunctions
{
	nut		        nutkWallFunction;
	k		        kLowReWallFunction;
	omega	        omegaWallFunction;
	epsilon	        epsilonWallFunction;
}

// ************************************************************************* //
```

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction).

**interfaceSpan** - *required* > Used for signed distance based reconstruction of the interface. Best results with 1.0.

**sdBasedLambda** - *optional* > Whether to calculate lambda based on the signed distance from an STL body surface.
> Default: *false*\
> Possible values: {*true*, *false*}

**stlName** - *required if sdBasedLambda* > Name (without the extension) of an STL file to read surface normals from.

**generateLambda** - *optional* > Whether to generate the lambda field based on an STL file.
> Default: *false*\
> Possible values: {*true*, *false*}

**cutCellType** - *optional* > ???
> Default: *cutCell*\
> Possible values: {*cutCell*, *cutEdges*}

**boundarySearch** - *optional* > Boundary cells found as face or vertex neighbors.
> Default: *face*\
> Possible values: {*face*, *vertex*}

**saveIntInfo** - *optional* > Whether to output info about interpolation.
> Default: *false*\
> Possible values: {*true*, *false*}

**correctIntPoints** - *optional* > ???
> Default: *false*\
> Possible values: {*true*, *false*}

**excludePatch** - *optional* > Selection of a patch name the cells of which are excluded from boundary cells.
> Default: *none*

**readSurfaceNormal** - *optional* > Whether to read surface normal field from file.
> Default: *false*\
> Possible values: {*true*, *false*}

**averageYOrtho** - *optional* > Whether to average orthogonal distance over face neighbors.
> Default: *false*\
> Possible values: {*true*, *false*}

**totalYOrthoAverage** - *optional* > ???
> Default: *false*\
> Possible values: {*true*, *false*}

**nAveragingYOrtho** - *optional* > For averageYOrtho. The number of averaging cycles.
> Default: *1*

**averagingCoeff** - *optional* > ???
> Default: *1.0*

**averageVolume** - *optional* > Whether to use average volume in reconstruction of distance from lambda field.
> Default: *false*\
> Possible values: {*true*, *false*}

**copyDisToInner** - *optional* > Whether to copy distance to boundary to inner cells.
> Default: *false*\
> Possible values: {*true*, *false*}

**scaleG** - *optional* > ???
> Default: *true*\
> Possible values: {*true*, *false*}

**scaleDisG** - *optional* > Whether to scale value of omega/epsilon and G in boundary cells.
> Default: *false*\
> Possible values: {*true*, *false*}

**scaleCoeff** - *optional* > For scaleDisG. Value of the scaling coefficient.
> Default: *1.0*

**scaleCoeffG** - *optional* > ???
> Default: *scaleCoeff*

**useEffectiveDist** - *optional* > Whether to use yEff instead of yOrtho for yPlus calculation.
> Default: *true*\
> Possible values: {*true*, *false*}

**assignNut** - *optional* > ???
> Default: *false*\
> Possible values: {*true*, *false*}

**uTauType** - *optional* > Method to calculate friction velocity.
> Default: *freeStreamCell*\
> Possible values: {*freeStreamCell*, *boundaryCell*, *effectiveDistance*, *cellDistance*, *coeffDistance*, *interpPoint*}

**uTauCoeff** - *optional* > Coefficient used in calculation of friction velocity. Method specific.
> Default: *1.0*

**wallFunctions** - *required if not laminar* > Subdict with specified wall functions.
> **nut** - *required* > wall function for nut.
>> Possible values: {*nutkWallFunction*}

> **k** - *required* > wall function for k.
>> Possible values: {*kLowReWallFunction*}

> **omega** - *required* > wall function for omega.
>> Possible values: {*omegaWallFunction*}

> **epsilon** - *required* > wall function for epsilon.
>> Possible values: {*epsilonWallFunction*}

## HFDIBRAS in turbulenceProperties

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  HFDIBRAS;


HFDIBRAS
{
    HFDIBRASModel	    HFDIBKOmega;
    turbulence          on;
    printCoeffs         on;
    disGSurfaceType     setValue;
    disGBoundaryValue   1.0;
    kSurfaceType        switched;
    kBoundaryValue      1.0;
}

// ************************************************************************* //
```

**HFDIBRASModel** - *required if HFDIBRAS* > Turbulence model used.
> Possible values: {*HFDIBKOmega*, *HFDIBKOmegaSST*, *HFDIBKEpsilon*, *HFDIBRealizableKE*}

**disGSurfaceType** - *required* > Type of surface field used for omega/epsilon and G.
> Possible values: {*setValue*, *switched*, *lambdaBased*}

**disGBoundaryValue** - *required* > For disGSurfaceType. Value defining the surface field based on the surface type.

**kSurfaceType** - *required* > Type of surface field used for turbulence kinetic energy.
> Possible values: {*setValue*, *switched*, *lambdaBased*}

**kBoundaryValue** - *required* > For kSurfaceType. Value defining the surface field based on the surface type.

**useKSource** - *optional* > Whether to use immersed boundary induced source term for k or matrix manipulation.
> Default: *true*\
> Possible values: {*true*, *false*}

## HFDIBSchemes in fvSchemes

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

.
.
.

gradSchemes
{
    .
    .
    .
    grad(lambda)        pointCellsLeastSquares;
}

divSchemes
{
    .
    .
    .
    div(((lambda*(1|A(U)))*f))  Gauss linear;
}

HFDIBSchemes
{
    outerSchemes
    {
        U               outer quadratic;
        k               switched quadratic logarithmic;
    }

    innerSchemes
    {
        U               cellPointFace;
        nut             cellPointFace;
        k               cellPointFace;
        omega           cellPointFace;
        epsilon         cellPointFace;
    }
}

.
.
.

// ************************************************************************* //
```

**grad(lambda)** - *required* > Gradient discretization scheme for the gradient of lambda.
> Possible values: native OpenFOAM gradient schemes

**div(((lambda\*(1|A(U)))\*f))** - *required* > Divergence scheme for the lambda term.
> Possible values: native OpenFOAM divergence schemes - Gauss linear is highly recommended

**outerSchemes** - *required* > Interpolation schemes for reconstruction of the surface value at the immersed boundary. Each field has a separate entry consisting of two tokens.
> **U** - *required* > Scheme for velocity.
>> Possible first tokens: {*outer*, *lambdaBased*, *switched*, *outerInner*, *inner*}\
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}. Two for *switched* and *outerInner*

> **k** - *required if present* > Schemes for turbulence kinetic energy.
>> Possible first tokens: {*outer*, *switched*, *outerInner*, *inner*}\
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}. Two for *switched* and *outerInner*

> **T** - *required if present* > Schemes for scalar T.
>> Possible first tokens: {*outer*}\
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}

**innerSchemes** - *required* > Interpolation schemes for interpolating values from cell centers to points. Separate entry for each field.
>> Possible values: native OpenFOAM interpolationSchemes

## HFDIB in SIMPLE in fvSolution

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

.
.
.

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;

    residualControl
    {
        p           1e-5;
    }

    HFDIB
    {
        U
        {
              surfaceType     setValue;
              boundaryValue   0.0;
        }
    }
}

.
.
.

// ************************************************************************* //
```

**U** - *required if present* > Subdict used by simpleHFDIBRANSFoam.
> **surfaceType** - *required* > Type of surface field used for velocity.
>> Possible values: {*setValue*, *switched*, *lambdaBased*}

> **boundaryValue** - *required* > For surfaceType. Value defining the surface field based on the surface type.

> **useNormSurface** - *optional* > Whether to correct surface field based on the surface normal.
>> Default: *false*\
>> Possible values: {*true*, *false*}

> **normalCorrectionLimit** - *optional* > Limiting value of lambda field for the correction of surface field by the surface normal.
>> Default: *0.5*

**T** - *required if present* > Subdict used by scalarHFDIBTransportFoam.
> **surfaceType** - *required* > Type of surface field used for scalar T.
>> Possible values: {*setValue*, *switched*, *lambdaBased*}

> **boundaryValue** - *required* > For surfaceType. Value defining the surface field based on the surface type.

## calculateWallShearStress
This utility function calculates the wall shear stress exerted by the fluid on the immersed boundary.

**Usage:**\
Run the utility from your case directory:
```bash
calculateWallShearStress
```
Alternatively, add the following line to your `Allrun` file:
```bash
runApplication calculateWallShearStress
```

**Output:**\
The results are written as a `tauw` file inside the respective time directories.

## calculateForceCoeffs
This utility function calculates the force coefficients acting on the immersed boundary.

Before running, ensure that your reference values are correctly defined in the `constant/HFDIBDEMDict` dictionary.\
Here is an example:
```cpp
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      HFDIBDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

.
.
.

forceCoeffs
{
    rhoInf      1000.0;
    magUInf     10.0;
    dragDir     (1 0 0);
    liftDir     (0 1 0);
    ARef        0.004;
}

// ************************************************************************* //
```

**rhoInf** - *required* > Freestream fluid density.

**magUInf** - *required* > Freestream velocity magnitude.

**dragDir** - *required* > Drag direction.

**liftDir** - *required* > Lift direction.

**ARef** - *required* > Reference area of the immersed body.

**Usage:**\
For laminar flows:
```bash
calculateForceCoeffsLaminar
```
For turbulent flows:
```bash
calculateForceCoeffs
```
Alternatively, add the following line to your `Allrun` file:
```bash
runApplication calculateForceCoeffsLaminar
# OR
runApplication calculateForceCoeffs
```

**Output:**\
The calculated coefficients are written to the `log.calculateForceCoeffs` file.
