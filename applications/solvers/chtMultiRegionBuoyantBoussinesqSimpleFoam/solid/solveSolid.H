{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix TEqn
        (
          - fvm::laplacian(alpha, T)
         ==
            fvOptions(T)
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        Info<< "\nSolving for solid region "
            << solidRegions[i].name() << endl;

        TEqn.solve();
        fvOptions.correct(T);
    }
}
