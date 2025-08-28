# Documentation for HFDIBDEMdict

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**interfaceSpan** - *required* > For signed distance based reconstruction of interface. Best results with 1.0

**wallFunctions** - *required if not laminar* > Subdict with specified wall functions
> **nut** - *required* > Wall function for nut
>> Possible values: {*nutkWallFunction*}

> **k** - *required* > Wall function for k
>> Possible values: {*kLowReWallFunction*}

> **omega** - *required* > Wall function for omega
>> Possible values: {*omegaWallFunction*}

> **epsilon** - *required* > Wall function for epsilon
>> Possible values: {*epsilonWallFunction*}

**boundarySearch** - *optional* > Boundary cells found as face or vertex neighbors. Default face
> Possible values: {*face*, *vertex*}

**saveIntInfo** - *optional* > If to output info about interpolation. Default false
> Possible values: {*true*, *false*}

**excludeWalls** - *optional* > If to exclude cells adjacent to walls from boundary cells. Default false
> Possible values: {*true*, *false*}

**readSurfaceNormal** - *optional* > If to read surface normal field from file. Default false
> Possible values: {*true*, *false*}

**averageYOrtho** - *optional* > If to average orthogonal distance over face neighbors. Default false
> Possible values: {*true*, *false*}

**nAveragingYOrtho** - *optional* > For averageYOrtho. The number of averaging cycles. Default 1.0

**averageVolue** - *optional* > If to use average volume in reconstruction of distance from lambda field. Default false
> Possible values: {*true*, *false*}

**copyDisToInner** - *optional* > If to copy distance to boundary to inner cells. Default false
> Possible values: {*true*, *false*}

**scaleDisG** - *optional* > If to scale value of omega/epsilon and G in boundary cells. Default false
> Possible values: {*true*, *false*}

**scaleCoeff** - *optional* > For scaleDisG. Value of the scaling coefficient. Default 1.0

**useEffectiveDistance** - *optional* > If to yEff instead of yOrtho for yPlus calculation. Default true
> Possible values: {*true*, *false*}

# Documentation for HFDIBSchemes subdict in fvSchemes
**outerSchemes** - *required* > Interpolation schemes for reconstruction of the surface value at the immerse boundary. Separate entry for each field consisting of two tokens
> **U** - *required* > Scheme for velocity
>> Possible first tokens: {*unifunctional*, *lambdaBased*, *switched*, *outerInner*, *inner*}
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}. Two for *switched* and *outerInner*

> **k** - *required if present* > Schemes for turbulence kinetic energy
>> Possible first tokens: {*unifunctional*, *switched*, *outerInner*, *inner*}
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}. Two for *switched* and *outerInner*

> **T** - *required if present* > Schemes for scalar T
>> Possible first tokens: {*unifunctional*}
>> Possible second tokens: {*constant*, *linear*, *quadratic*, *logarithmic*, *fixedGradient*, *zeroGradient*}

**innerSchemes** - *required* > Interpolation schemes for interpolating values from cell centers to points. Separate entry for each field
>> Possible values: native openfoam interpolationSchemes

# Documentation for HFDIB subdict in SIMPLE dict in fvSolution
**U** - *required if present* > Subdict used by simpleHFDIBRANSFoam
> **surfaceType** - *required* > Type of surface field used for velocity
>> Possible values: {*setValue*, *switched*, *lambdaBased*}

> **boundaryValue** - *required* > For surfaceType. Value defining the surface field based on the surface type

> **tolEqn** - *required* > Tolerance for convergence of U to UIB in boundary cells

> **maxEqnIters** - *required* > For tolEqn. Maximum number of iterations to drop below required tolerance

> **cutForce** - *optional* > If to cut the induced source term to surface normal direction. Used for testing. Default false
>> Possible values: {*true*, *false*}

> **cutVelocity** - *optional* > If to cut velocity when it goes in oposity direction to surface normal. Used for testing. Default false
>> Possible values: {*true*, *false*}

> **cutPhi** - *optional* > If to cut fluxes inside immersed boundary. Used for testing. Default false
>> Possible values: {*true*, *false*}

> **enforceVelocity** - *optional* > If to enforce velocity values in boundary cell directly. Used for testing. Default false
>> Possible values: {*true*, *false*}

**T** - *required if present* > Subdict used by scalarHFDIBTransportFoam
> **surfaceType** - *required* > Type of surface field used for scalar T
>> Possible values: {*setValue*, *switched*, *lambdaBased*}

> **boundaryValue** - *required* > For surfaceType. Value defining the surface field based on the surface type

> **tolEqn** - *required* > Tolerance for convergence of T to TIB in boundary cells

> **maxEqnIters** - *required* > For tolEqn. Maximum number of iterations to drop below required tolerance

# Documentation for HFDIBRAS subdict in turbulenceProperties -- all entries from native RAS subdict present
**kSurfaceType** - *required* > Type of surface field used for turbulence kinetic energy
> Possible values: {*setValue*, *switched*, *lambdaBased*}

**kBoundaryValue** - *required* > For kSurfaceType. Value defining the surface field based on the surface type

**disGSurfaceType** - *required* > Type of surface field used for omega/epsilon and G
> Possible values: {*setValue*, *switched*, *lambdaBased*}

**disGBoundaryValue** - *requied* > Fot disGSurfaceType. Value defining the surface field based on the surface type

**tolKEqn** - *required* > Tolerance for convergence of k to kIB in boundary cells

**maxKEqnIters** - *required* > For tolKEqn. Maximum number of iterations to drop below required tolerance

**useKSource** - *optional* > If to use immersed boundary induced source term for k. Default true
> Possible values: {*true*, *false*}
