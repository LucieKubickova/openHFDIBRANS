// Momentum predictor

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

HFDIBRANS.computeUi(U, Ui, turbulence->k(), turbulence->nu()); // problem for laminar, or... is it?

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) + f);
    fvOptions.correct(U);
}
