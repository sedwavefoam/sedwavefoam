{
    word betaScheme("div(phi,beta)");
    word betarScheme("div(phirb,beta)");

    beta = 1 - gamma;
 
    tmp<fv::ddtScheme<scalar> > ddtBeta
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(beta)")
        )
    );

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    if
    (
        isType<fv::EulerDdtScheme<scalar> >(ddtBeta())
     || isType<fv::localEulerDdtScheme<scalar> >(ddtBeta())
    )
    {
        ocCoeff = 0;
    }
	// if you are using openfoam over 2.4.x you can use crank-nicolson 
    /*else if (isType<fv::CrankNicolsonDdtScheme<scalar> >(ddtBeta()))
    { 
	if (nBetaSubCycles > 1)
        {
            FatalErrorIn(args.executable())
                << "Sub-cycling is not supported "
                   "with the CrankNicolson ddt scheme"
                << exit(FatalError);
        }

        ocCoeff =
            refCast<fv::CrankNicolsonDdtScheme<scalar> >(ddtBeta()).ocCoeff();
    }
    */
    else
    {
        FatalErrorIn(args.executable())
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phico(interface.cAlpha1()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icBeta > 0)
    {
        phico *= (1.0 - icBeta);
        phico += (interface.cAlpha1()*icBeta)*fvc::interpolate(mag(U));
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phico.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phico.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }

    if (MULESCorr)
    {
        fvScalarMatrix betaEqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(beta)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(beta)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, beta)
        );

        betaEqn.solve();

        tmp<surfaceScalarField> tphiBetaUD(betaEqn.flux());
        phiBeta = tphiBetaUD();

        if (betaApplyPrevCorr && tphiBetaCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            #ifdef LTSSOLVE
            MULES::LTScorrect(beta, phiBeta, tphiBetaCorr0(), 1, 0);
            #else
            MULES::correct(beta, phiBeta, tphiBetaCorr0(), 1, 0);
            #endif

            phiBeta += tphiBetaCorr0();
        }

        // Cache the upwind-flux
        tphiBetaCorr0 = tphiBetaUD;

        gamma = 1.0 - beta;

        interface.correct();
    }

    for (int aCorr=0; aCorr<nBetaCorr; aCorr++)
    {
        surfaceScalarField phirc(phico*interface.nHatf());

        tmp<surfaceScalarField> tphiBetaUn
        (
            fvc::flux
            (
                phi,
                beta,
                betaScheme
            )
          + fvc::flux
            (
               //-fvc::flux(-phirc, gamma, betarScheme),
                -fvc::flux(-phirc/fvc::interpolate(1.0-alpha), gamma, betarScheme), //check Eq 18 in the first sedwavefoam paper
		beta,
                betarScheme
            )
        );

        // Calculate the Crank-Nicolson off-centred alpha flux
        if (ocCoeff > 0)
        {
            tphiBetaUn =
                cnCoeff*tphiBetaUn + (1.0 - cnCoeff)*phiBeta.oldTime();
        }

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiBetaCorr(tphiBetaUn() - phiBeta);
            volScalarField beta0("beta0", beta);

            #ifdef LTSSOLVE
            MULES::LTScorrect(beta, tphiBetaUn(), tphiBetaCorr(), 1, 0);
            #else
            MULES::correct(beta, tphiBetaUn(), tphiBetaCorr(), 1, 0);
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiBeta += tphiBetaCorr();
            }
            else
            {
                beta = 0.5*beta + 0.5*beta0;
                phiBeta += 0.5*tphiBetaCorr();
            }
        }
        else
        {
            phiBeta = tphiBetaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve(beta, phi, phiBeta, 1, 0);
            #else
            MULES::explicitSolve(beta, phiCN, phiBeta, 1, 0);
            #endif
        }

        gamma = 1.0 - beta;

        interface.correct();
    }

    if (betaApplyPrevCorr && MULESCorr)
    {
        tphiBetaCorr0 = phiBeta - tphiBetaCorr0;
    }

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
//        rhoPhi = phiBeta*(rhob - rhoc) + phiCN*rhoc;
	rhoPhi = phiBeta*(rhob - rhoc);
	phib = phiCN; //UEqn and check
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            phiBeta = (phiBeta - (1.0 - cnCoeff)*phiBeta.oldTime())/cnCoeff;
        }

     //   rhoPhi = phiBeta*(rhob - rhoc) + phib*rhoc;
   	    rhoPhi = phiBeta*(rhob-rhoc);
        }

    Info<< "Phase-1 volume fraction = "
        << beta.weightedAverage(mesh.Vsc()).value()
        << "  Min(beta) = " << min(beta).value()
        << "  Max(beta) = " << max(beta).value()
        << endl;
}
