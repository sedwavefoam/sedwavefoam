    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<phaseModel> phasea = phaseModel::New
    (
        mesh,
        transportProperties,
        "a"
    );

    autoPtr<phaseModel> phaseb = phaseModel::New
    (
        mesh,
        transportProperties,
        "b"
    );

    autoPtr<phaseModel> phasec = phaseModel::New
    (
        mesh,
        transportProperties,
        "c"
    );

    volVectorField& Ua = phasea->U();
    surfaceScalarField& phia = phasea->phi();
    const dimensionedScalar& rhoa = phasea->rho();
    const dimensionedScalar& nua = phasea->nu();
    const dimensionedScalar& da = phasea->d();

    volVectorField& Ub = phaseb->U();
    surfaceScalarField& phib = phaseb->phi();
    const dimensionedScalar& rhob = phaseb->rho();
    const dimensionedScalar& nub = phaseb->nu();
    const dimensionedScalar& db = phaseb->d();

    volVectorField& Uc = phasec->U(); //dummy
    surfaceScalarField& phic = phasec->phi(); //dummy
    const dimensionedScalar& rhoc = phasec->rho();
    const dimensionedScalar& nuc = phasec->nu();
    const dimensionedScalar& dc = phasec->d(); //dummy
    
    //const dimensionedScalar& Yar = dimensionSet(1,-2,-1,0,0,0,0);
    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh 
    );

    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        scalar(1)-alpha-alpha1
    );

    gamma = max(gamma,scalar(0));

    Info<<"after reading alpha, alpha1, gamma, now recalculate alpha1"<<endl;

    alpha1 == scalar(1)-alpha-gamma;

    volScalarField beta 
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1)-gamma
       ,alpha.boundaryField().types()
    );

////// YAROOO
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
////////////////////////
//
    volVectorField U //used in waves2Foam
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ub
    );

    dimensionedScalar Cvm
    (
        transportProperties.lookup("Cvm")
    );

    dimensionedScalar Cl
    (
        transportProperties.lookup("Cl")
    );

    dimensionedScalar Ct
    (
        transportProperties.lookup("Ct")
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh
        ),
        fvc::interpolate(alpha)*phia + fvc::interpolate(1.0-alpha)*phib
    );

    volScalarField nugl
    (
        IOobject
        (
            "nugl",
            runTime.timeName(),
            mesh
        ),
        (alpha1*rhob*nub + gamma*rhoc*nuc)/(alpha1*rhob+gamma*rhoc)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        (alpha1*rhob + gamma*rhoc)/(alpha1+gamma)
    );

    
	surfaceScalarField rhoPhiU
    (
        IOobject
        (
            "rhoPhiU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phib
	//fvc::interpolate(rho)*phib
    );
	
	surfaceScalarField rhoPhiKE
    (
        IOobject
        (
            "rhoPhiKE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phib
	//fvc::interpolate(rho)*phib
    );
	
	surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phib
    );

    #include "createRASTurbulence.H"

    Info<< "Calculating field DDtUa and DDtUb\n" << endl;

    volVectorField DDtUa =
        fvc::ddt(Ua)
      + fvc::div(phia, Ua)
      - fvc::div(phia)*Ua;

    volVectorField DDtUb =
        fvc::ddt(Ub)
      + fvc::div(phib, Ub)
      - fvc::div(phib)*Ub;


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    Info<< "Reading field p_rgh\n" << endl;


    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    };

    fv::IOoptionList fvOptions(mesh);

    IOdictionary interfacialProperties
    (
        IOobject
        (
            "interfacialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<dragModel> draga = dragModel::New
    (
        interfacialProperties,
        alpha,
        alpha1,
        phasea,
        phaseb,
        phasec
    );

    autoPtr<dragModel> dragb = dragModel::New
    (
        interfacialProperties,
        alpha1,
        alpha,
        phaseb,
        phasea,
        phasec
    );

    word dragPhase("blended");
    if (interfacialProperties.found("dragPhase"))
    {
        dragPhase = word(interfacialProperties.lookup("dragPhase"));

        bool validDrag =
            dragPhase == "a" || dragPhase == "b" || dragPhase == "blended";

        if (!validDrag)
        {
            FatalErrorIn(args.executable())
                << "invalid dragPhase " << dragPhase
                << exit(FatalError);
        }
    }

    Info << "dragPhase is " << dragPhase << endl;
    kineticTheoryModel kineticTheory
    (
        phasea,
        Ub,
        alpha,
        draga
    );

    surfaceScalarField rUaAf
    (
        IOobject
        (
            "rUaAf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), 0.0)
    );

    volScalarField rUaA
    (
        IOobject
        (
            "rUaA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), 0.0)
    );

    surfaceScalarField ppMagf
    (
        IOobject
        (
            "ppMagf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    );

    volScalarField lambdaUa
    (
        IOobject
        (
            "lambdaUa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    );
    volScalarField K
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        draga->K(mag(Ub-Ua))
    );

    volVectorField gradPs
    (
        IOobject
        (
            "gradPs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
       dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
    );

    volVectorField gradPf
    (
        IOobject
        (
            "gradPf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
    );

    // Construct interface from gamma distribution
    interfaceProperties interface(beta,gamma, Ub, transportProperties);

   
     // MULES flux from previous time-step
    surfaceScalarField phiBeta
    (
        IOobject
        (
            "phiBeta",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(beta)
    );
 
    tmp<surfaceScalarField> tphiBetaCorr0;
    relaxationZone relaxing(mesh, Ub, beta); // we are defining the relaxation of Ub as relaxing
    relaxationZone relaxing2(mesh, Ua, beta); // we are defining the relaxation of Ua as relaxing2
