scalar pi = constant::mathematical::pi;

// following is the 2nd order stokes waves
gradPOSC = gradPMEAN+gradPAMP1*Foam::cos(initTheta*(pi)/180.0 + 2.0*(pi)*runTime.value()/oscpT)+gradPAMP2*Foam::sin(initTheta*(pi)/90.0 + 4.0*(pi)*runTime.value()/oscpT);     // Add oscilitary pressure gradient)

fvVectorMatrix UaEqn(Ua, rhoa.dimensions()*Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, rhoa.dimensions()*Ub.dimensions()*dimVol/dimTime);

{
    rho == (gamma*rhoc + alpha1*rhob)/(scalar(1.0)-alpha);
    nugl=(alpha1*rhob*nub + gamma*rhoc*nuc)/(gamma*rhoc + alpha1*rhob);  //// YAR
		
    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions(); 
    {
        volTensorField gradUbT = fvc::grad(Ub)().T();
        volTensorField gradUb = fvc::grad(Ub);
        nuEffb = nugl + nutb;//sgsModelb->nuSgs();
		
        volTensorField Rcb
        (
            "Rcb",
            -nuEffb*dev(T(fvc::grad(Ub))) 
        );
        volVectorField URb = -(rho)*fvc::grad(alpha)/fvc::average(1.0-alpha);

        surfaceScalarField phiRb =
            fvc::interpolate(rho*nuEffb)*mesh.magSf()*(fvc::snGrad(alpha))
            /fvc::interpolate(1.0-alpha);
        volVectorField susb
        (
            "susb",
            // turbulence suspension in all directions
            SUS*(K)*nutb*fvc::grad(alpha)
        );

	rhoPhiU == rhoPhi/fvc::interpolate(1.0-alpha) + phib*rhoc; // this is what I explained



	UbEqn =
        (
	  fvm::ddt(rho,Ub)
	  + fvm::div(rhoPhiU,Ub,"div(phib,Ub)")
	  - fvm::Sp(rho*fvc::div(phib), Ub)
      - fvm::laplacian(rho*nuEffb, Ub)
      + fvc::div(rho*Rcb)
      // Stress part2
      + fvm::div(phiRb, Ub, "div(phiRb,Ub)")
      - fvm::Sp(fvc::div(phiRb), Ub)
      + (URb & Rcb)
         ==
            - fvm::Sp(alpha*K, Ub)   // explicit part of Ua has been obtained in p-equation
            + (alpha)*K*Ua               // Explicit drag transfered to p-equation
            + alpha*(liftCoeff + Cvm*rho*DDtUa)
            + susb                        // add the suspension term by C.Z(Mar,2011)
        );
        UbEqn.relax();
	fvOptions.constrain(UbEqn);
    }

    {
        volTensorField Rca
        (
            "Rca",
            -nuEffa*dev(T(fvc::grad(Ua)))
        );

        if (kineticTheory.on())
        {
            Rca -= ((lambdaUa/((alpha+1e-8)*rhoa))*tr(gradUaT))*tensor(I);
        }

        surfaceScalarField phiRa =
            -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
            /fvc::interpolate((alpha+ scalar(1e-8)));

        volVectorField susa
        (
            "susa",
              // turbulent suspension in all directions
              SUS*(K)*nutb*((1.0-alpha)/max(alpha,scalar(1e-10)))*(fvc::grad(alpha))
        );

        UaEqn = // rhoa is not changing anyway
        (
            (scalar(1.0) + Cvm*rhob*(1.0-alpha)/rhoa)*
            (
                fvm::ddt(rhoa,Ua)
              + fvm::div(rhoa*phia, Ua, "div(phia,Ua)")
              - fvm::Sp(rhoa*fvc::div(phia), Ua)
            )
          // The following 5 lines are treatment for granular shear stress term.
          + fvc::div(rhoa*Rca)
          // The following to include interfacial effect
          - fvm::laplacian(rhoa*nuEffa, Ua)
          + fvm::div(rhoa*phiRa, Ua, "div(phia,Ua)")
          - fvm::Sp(rhoa*fvc::div(phiRa), Ua)
          + rhoa*(fvc::grad(alpha)/(alpha+scalar(1e-8)) & Rca)
         ==
          - fvm::Sp((1.0-alpha)*K, Ua) //explit part of Ub is obtained in p-equation
          + ((1.0-alpha))*K*Ub             // Explicit drag transfered to p-equation
          - (1.0-alpha)*(liftCoeff - Cvm*rho*DDtUb)
          - susa
        );
          
        UaEqn.relax();
    }
}
