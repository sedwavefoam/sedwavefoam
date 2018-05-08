{
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField voidfractionf(scalar(1) - alphaf);
    surfaceScalarField rhof(fvc::interpolate(rho));

    rUaA=1.0/UaEqn.A();
    volScalarField rUbA(1.0/UbEqn.A());

    rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf(fvc::interpolate(rUbA));

    forAll (alpha, celll)
    {
     if(alpha[celll]<=1e-6)
       {
          Ua[celll] = Ub[celll];
       }
    }

    phia = fvc::interpolate(Ua) & mesh.Sf();
    phib = fvc::interpolate(Ub) & mesh.Sf();
			
    surfaceScalarField phiDraga
    (
      - (fvc::interpolate(scalar(1.0)/(((alpha+scalar(1e-6)))))*fvc::snGrad(kineticTheory.pa()))*rUaAf*mesh.magSf()
       -(fvc::interpolate(scalar(1.0)/(((alpha+scalar(1e-6)))))*fvc::snGrad(kineticTheory.pf_()))*rUaAf*mesh.magSf()
      +rUaAf*fvc::interpolate(rhoa-rho)*(g & mesh.Sf())
      +rUaAf*((gradPOSC) & mesh.Sf())
    );

    surfaceScalarField phiDragb
    (
       rUbAf*(gradPOSC & mesh.Sf())
       -rUbAf*fvc::interpolate(interface.sigmaK())*fvc::snGrad(gamma)*mesh.magSf()
       -(rUbAf)*ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    Ua = rUaA*UaEqn.H();
    Ub = rUbA*UbEqn.H();

    #include "UaUbLimiter.H"

    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (isA<symmetryFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (p.boundaryField().types()[patchi] == "symmetryPlane")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (p.boundaryField().types()[patchi] == "zeroGradient")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }
// for version later than OF2.3.0
    surfaceScalarField phiUa = (fvc::interpolate(Ua) & mesh.Sf()) + fvc::interpolate(rhoa*rUaA)*fvc::ddtCorr(Ua, phia);
    surfaceScalarField phiUb = (fvc::interpolate(Ub) & mesh.Sf()) + fvc::interpolate(rho*rUbA)*fvc::ddtCorr(Ub, phib);

    adjustPhi(phiUa, Ua, p_rgh);
    adjustPhi(phiUb, Ub, p_rgh);

    phia = phiUa + phiDraga;
    phib = phiUb + phiDragb;

    phi = alphaf*phia + voidfractionf*phib;

    surfaceScalarField Dp
    (
        "(rho*(1|A(U)))",
        alphaf*rUaAf + voidfractionf*rUbAf
    );

    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
        phi.boundaryField()
        -  (alphaf.boundaryField()*(mesh.Sf().boundaryField() & Ua.boundaryField())
            +voidfractionf.boundaryField()*(mesh.Sf().boundaryField() & Ub.boundaryField())
           )
        ) / (mesh.magSf().boundaryField()*Dp.boundaryField())

    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p_rgh) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        p_rgh.correctBoundaryConditions();

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phia -= rUaAf*SfGradp;
            phib -= rUbAf*SfGradp;
            phi -= pEqn.flux();
            p.relax();
        }
    }
    Ua += rUaA*fvc::reconstruct((phia - phiUa)/rUaAf);
    Ua.correctBoundaryConditions();
    Ub += rUbA*fvc::reconstruct((phib - phiUb)/rUbAf);
    Ub.correctBoundaryConditions();
    forAll (alpha, celll)
    {
     if(alpha[celll]<=1e-6)
       {
          Ua[celll] = Ub[celll];
       }
    }
    #include "UaUbLimiter.H"
    Ua.correctBoundaryConditions();
    U = Ub; //used in waves2Foam
    fvOptions.correct(U);
    #include "continuityErrs.H"
    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}

