    volVectorField Ur(-Ua + Ub);
    volScalarField magUr(mag(Ur));

    volScalarField Ka(draga->K(magUr));
    volScalarField K(Ka);

    if (dragPhase == "b")
    {
        volScalarField Kb(dragb->K(magUr));
        K = Kb;
    }
    else if (dragPhase == "blended")
    {
        volScalarField Kb(dragb->K(magUr));
        K = ((1-alpha)*Ka + alpha*Kb);
    }

    volVectorField liftCoeff(Cl*(gamma*rhoc+alpha1*rhob + alpha*rhoa)*(Ur ^ fvc::curl(U)));
