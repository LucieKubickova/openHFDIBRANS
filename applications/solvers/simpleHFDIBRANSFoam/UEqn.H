// Momentum predictor

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvOptions(U)
  + g
);

UEqn.relax();

fvOptions.constrain(UEqn);

// HFDIBRANS: set up imposed velocity
HFDIBRANS.computeUi(U, Ui, surfaceType);
Ui.correctBoundaryConditions();

// HFDIBRANS: set up surfaces
HFDIBRANS.updateSurface(surface, surfaceType);
surface.correctBoundaryConditions();
HFDIBRANS.correctSurfaceByNormal(normSurface, surface, normCorrLimit);
normSurface.correctBoundaryConditions();

// HFDIBRANS: immersed boundary induced source term
if (useNormSurface)
{
    f = normSurface & (UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
}

else
{
    f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
}

// HFDIBRANS: add to U equation
UEqn -= f;

// HFDIBRANS: calculate f from under-relaxation // Note (LK): does not work well
//~ scalar alpha = mesh.equationRelaxationFactor("U");
//~ f = surface*(1 - alpha)/alpha*UEqn.A()*(Ui - U);
//~ UEqn -= f;

// HFDIBRANS: matrix manipulate
//~ matrixManipulate(UEqn, Ui, surface);
//~ {
    //~ DynamicList<label> cells;
    //~ DynamicList<vector> phis;

    //~ forAll(lambda, cellI)
    //~ {
        //~ if (surface[cellI] == 1.0)
        //~ {
            //~ cells.append(cellI);
            //~ phis.append(Ui[cellI]);
        //~ }
    //~ }

    //~ labelUList UCells(cells);
    //~ UList<vector> UPhis(phis);

    //~ UEqn.setValues(UCells, UPhis);
//~ }

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}

// calculate error
scalar error(0.0);
if (useNormSurface)
{
    error = max(mag(normSurface & (Ui - U))).value();
}

else
{
    error = max(mag(surface*(Ui - U))).value();
}
Info << "HFDIBRANS: Max error in U -> Ui is " << error << endl;

// apply correction
if (useNormSurface)
{
    U += 1.0*normSurface & (Ui - U);
}

else
{
    U += 1.0*surface*(Ui - U);
}
