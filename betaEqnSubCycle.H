if (nBetaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> betaSubCycle(alpha1, nBetaSubCycles);
        !(++betaSubCycle).end();
    )
    {
        #include "betaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "betaEqn.H"
}
