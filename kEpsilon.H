if(turbulence)
{
  if (kEpsilon)
    {
      if (mesh.changing())
	{
	  y.correct();
	}
      // to be consistent, let's define nutb = C_mu*k_f^2/epsilon_f
	  
	  rho == (gamma*rhoc + alpha1*rhob)/(scalar(1.0)-alpha); //// YAR
      nugl=(alpha1*rhob*nub + gamma*rhoc*nuc)/(gamma*rhoc + alpha1*rhob);  //// YAR
	
	
      tmp<volTensorField> tgradUb = fvc::grad(Ub);
      volTensorField TgradUb = fvc::grad(Ub)().T();
      // same as that in momentum equation
      volScalarField G = fvc::grad(Ub) && (-(2.0/3.0)*(k + nutb*tr(TgradUb))*I + nutb*(tgradUb+TgradUb));
      tgradUb.clear();
      dimensionedScalar Tpsmall
	(
	 "Tpsmall",
	 dimensionSet(1, -3, -1, 0, 0, 0, 0),
	 scalar(1e-40)
	 );                                  
      dimensionedScalar kSmall
	(
	 "ksmall",
	 dimensionSet(0, 2, -2, 0, 0, 0, 0),
	 scalar(1e-10)// this should be the same as your initial condition
	 );
      dimensionedScalar epsilonSmall
	(
	 "epsilonsmall",
	 dimensionSet(0, 2, -3, 0, 0, 0, 0),
	 scalar(1e-10) // this should be the same as your initial condition 
	 );
      volScalarField tmf = Foam::exp(-B*rhoa*6.0*epsilon/(max(k,kSmall)*max((1.0-alpha)*K,Tpsmall)));
      // just need the horizontal ones;
      volScalarField rtemp = (Ub-Ua) & g;
      volVectorField Urtemp = rtemp*g/magSqr(g);
      volVectorField Urh = (Ub-Ua) - Urtemp;
      volScalarField ESD1 = nutb*(K/rho)*(fvc::grad(alpha)) & Urh;
      volScalarField ESD2 = -2.0*K*alpha*(1.0-tmf)/rho; //k_f was included in the equation implicitly
      volScalarField ESD3 = K*alpha*magSqr(Ub-Ua)/((1.0-alpha)*rho);
      volScalarField ESD4 = -SUS*(nutb/(1.0-alpha))*(rhoa/rho-1.0)*g & fvc::grad(alpha);
      ESD4 = ESD4 + KE1*ESD1;
      // the following term is implemented in epsilon equation.
      volScalarField ESD5 = -SUS*(nutb/(1.0-alpha))*(rhoa/rho-1.0)*g & fvc::grad(alpha);
      // if ESD5 > 0, then it should be implemented, else, no
      ESD5 = max(ESD5,scalar(0.0)*ESD5);
      ESD5 = ESD5 + KE1*ESD1;
      volScalarField ESD = KE3*ESD3;   //need to check here!
      
      rhoPhiKE == rhoPhi/fvc::interpolate(1.0-alpha) + phib*rhoc;	
  
      #include "wallFunctions.H"
      
      // Dissipation equation
      //epsilon = epsilon_f
      fvScalarMatrix epsEqn
	(
          fvm::ddt(rho,epsilon)
	    + fvm::div(rhoPhiKE, epsilon)
	

	- fvm::laplacian
        (
          rho*(alphaEps*nutb+nugl), epsilon,
          "laplacian(DepsilonEff,epsilon)"
        )
	==
	fvm::Sp(C1*rho*G/max(k,kSmall) ,epsilon)
	// The following one line should be Sp for stability
	- fvm::SuSp(C2*rho*epsilon/max(k,kSmall), epsilon)
	+ fvm::Sp(C3ep*rho*ESD/max(k,kSmall),epsilon)
	+ fvm::Sp(C3ep*rho*KE2*ESD2,epsilon)
	+ fvm::Sp((C4ep*rho*KE4*ESD5/max(k,kSmall)),epsilon)
	);
      
      #include "wallDissipation.H"
      
      epsEqn.relax();
      epsEqn.solve();

      epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), epsilonSmall.value()));
      epsilon.correctBoundaryConditions();

      // Turbulent kinetic energy equation
      // k=k_f 
      fvScalarMatrix kEqn
	(
	 fvm::ddt(rho,k)
	 + fvm::div(rhoPhiKE, k)
	 - fvm::laplacian // RHS2
	 (
	  rho*(alphak*nutb+nugl), k,
	  "laplacian(DkEff,k)"
	  )
	 ==
	 fvm::Sp(rho*G/max(k,kSmall),k)
	 // The following one line should be Sp for stability
	 - fvm::Sp(rho*epsilon/max(k,kSmall),k) // 3rd RHS
	 + fvm::Sp(ESD*rho/max(k,kSmall),k)
	 + fvm::Sp(KE4*ESD4*rho/max(k,kSmall),k) // the lst
	 + fvm::Sp(KE2*ESD2*rho,k) //4th 
	 );
      kEqn.relax();
      kEqn.solve();

      k.max(dimensionedScalar("zero", k.dimensions(), kSmall.value()));
      k.correctBoundaryConditions();
      //- Re-calculate turbulence viscosity
      nutb = Cmu*sqr(k)/max(epsilon,epsilonSmall);
      #include "wallViscosity.H"
    }
  else if (MixingLength)
    {
//Efective fluid viscosity
        tmp<volTensorField> gradUb = fvc::grad(Ub);
        volSymmTensorField Db = symm(gradUb);
        volScalarField magDb = ::sqrt(2.0)*mag(Db);	
        gradUb.clear();

        volVectorField centres = Ub.mesh().C();
        volScalarField Y = centres.component(1);
        scalarField& alphaCells = alpha.internalField();
        scalarField& magDbCells = magDb.internalField();
        scalarField& YCells = Y.internalField();

        scalar Lm = 0.;
        scalar dY;
        scalar kappaLMs=kappaLM.value();
        scalar alphaMaxLMs = alphaMaxLM.value();
        scalar LmPhi = 0.;

        nutb.storePrevIter();

        forAll(YCells, cellI)
        {
          dY = YCells[cellI]-YCells[cellI-1];
          LmPhi = LmPhi + kappaLMs*max(alphaMaxLMs-alphaCells[cellI],0.)/alphaMaxLMs*dY;
          Lm = LmPhi;
          nutb[cellI] = pow(Lm,2)*magDbCells[cellI];
         }
        nutb.relax();
    }
}
Info<< "max(k) = " << max(k).value() <<" max(nutb) = "<<max(nutb).value()<< endl;

// limit nutb for stability reason, 
nutb.min(dimensionedScalar("nutlimit", nub.dimensions(), 1));
nutb.correctBoundaryConditions();

// NuEffb corresponds to the eddy viscosity + the molecular diffusion (for the granular rheology nub will be replace by the effective viscosity in callFrictionStress.H)
nuEffb = (nutb + nugl);
