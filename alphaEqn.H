{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phia - phib);

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phic, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, (1.0-alpha), schemer), alpha, schemer)
        );

        alphaEqn.relax();
        alphaEqn.solve();
        alpha.correctBoundaryConditions();
        int i=0;
        while (i<500 && max(alpha).value()>=(alphaMax-0.001))
        {
             #include "packingLimiter.H"
             i++;
        }
        alpha.min(alphaMax - 0.0001);
        alpha.max(0.0);

        Info<< "Dispersed phase volume fraction = "
            << alpha.weightedAverage(mesh.V()).value()
            << "  Min(alpha) = " << min(alpha).value()
            << "  Max(alpha) = " << max(alpha).value()
            << endl;
    }
}

