{
    volScalarField alphat("alphat", turb.nut()/Prt);
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turb.nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rho = 1.0 - beta*(T - TRef);

    kappa = alphaEff*cpRef*rhoRef;
}
